//! Air Quality and Pollen Module
//!
//! Fetches Air Quality Index (European AQI) and pollen data from the
//! Open-Meteo air-quality API and classifies pollen concentrations into
//! medical/aerobiology severity bands.
//!
//! The HTTP response is cached on disk (default `$XDG_RUNTIME_DIR/suntimes/<lat>_<lon>.json`)
//! with a 1-hour TTL. The cache is invalidated and re-fetched if the file
//! is missing, older than the TTL, empty, or unparseable as JSON.
//! Writes to the cache are atomic (write to temp file + rename) so a
//! concurrently-running invocation never sees a half-written file.
//!
//! The `<lat>_<lon>` portion of the filename is governed by
//! `--location-privacy` (and the user's input precision) — see
//! `cli::Coordinate::for_api`. This module takes pre-formatted
//! coordinate strings and uses them verbatim for both the URL query
//! and the cache filename, so URL and cache key are always byte-equal.
//!
//! Every fetch requests *all* fields (AQI + all pollen species)
//! regardless of which `--aqi` / `--pollen` flag is set, so a cached
//! response is always complete; the display layer is what filters down
//! to what the user asked for. This avoids surprises where running
//! `--pollen` first and `--aqi` later within the cache window would
//! otherwise leave AQI fields missing for up to an hour.
//!
//! Proxy detection is automatic via `ureq::Proxy::try_from_env()` which
//! consults `ALL_PROXY`, `HTTPS_PROXY`, `HTTP_PROXY` (and `NO_PROXY`).
//!
//! Pollen severity thresholds follow the European Aeroallergen Network
//! (EAN) / RNSA conventions widely used in clinical aerobiology:
//! values below the species' "Low" minimum (1 grain/m³ for all six
//! species we query) are treated as trace and hidden from output.

use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};
use tempfile::NamedTempFile;

// ===================== CONSTANTS =====================

const API_BASE: &str = "https://air-quality-api.open-meteo.com/v1/air-quality";
/// Subdirectory under `dirs::runtime_dir()` where per-location caches
/// live. Created on demand by `write_cache_atomic`.
const CACHE_SUBDIR: &str = "suntimes";
const CACHE_TTL: Duration = Duration::from_secs(3600);
const HTTP_TIMEOUT: Duration = Duration::from_secs(5);

/// AQI fields (queried on every fetch).
///
/// `european_aqi` is the overall European AQI index.
/// `sulphur_dioxide`, `pm10`, `pm2_5`, `nitrogen_dioxide`, `carbon_monoxide`,
/// `ozone`, and `dust` are raw concentrations in μg/m³.
/// `aerosol_optical_depth` is dimensionless (column extinction at
/// 550 nm); we map it to a descriptive category for display.
const AQI_FIELDS: &[&str] = &[
    "european_aqi",
    "pm10",
    "pm2_5",
    "nitrogen_dioxide",
    "carbon_monoxide",
    "ozone",
    "sulphur_dioxide",
    "dust",
    "aerosol_optical_depth",
];

// ===================== POLLEN TYPES =====================

/// Severity tier for a pollen reading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PollenSeverity {
    Low,
    Moderate,
    High,
    VeryHigh,
}

impl PollenSeverity {
    pub fn label(&self) -> &'static str {
        match self {
            PollenSeverity::Low => "Low",
            PollenSeverity::Moderate => "Moderate",
            PollenSeverity::High => "High",
            PollenSeverity::VeryHigh => "Very High",
        }
    }
}

/// A pollen species we query, plus its display name and severity thresholds.
///
/// `thresholds[i]` is the *minimum* concentration (grains/m³) required to
/// reach severity tier `i`:
/// * 0: Low
/// * 1: Moderate
/// * 2: High
/// * 3: Very High
///
/// Values below `thresholds[0]` are treated as trace and hidden.
///
/// Reference scales: European Aeroallergen Network (EAN) / RNSA
/// (Réseau National de Surveillance Aérobiologique). These match the
/// breakpoints commonly used by national pollen services across Europe.
pub struct PollenSpecies {
    /// Open-Meteo API field name (e.g. `birch_pollen`)
    pub api_key: &'static str,
    /// Human-readable name (e.g. `Birch`)
    pub display_name: &'static str,
    /// Severity tier minima in grains/m³
    pub thresholds: [f64; 4],
}

/// All six pollen species we report on.
pub const POLLEN_SPECIES: &[PollenSpecies] = &[
    // Tree pollens — alder and birch share the same scale (high counts common).
    PollenSpecies {
        api_key: "alder_pollen",
        display_name: "Alder",
        thresholds: [1.0, 11.0, 101.0, 1001.0],
    },
    PollenSpecies {
        api_key: "birch_pollen",
        display_name: "Birch",
        thresholds: [1.0, 11.0, 101.0, 1001.0],
    },
    // Grass — much lower clinical thresholds than tree pollen.
    PollenSpecies {
        api_key: "grass_pollen",
        display_name: "Grass",
        thresholds: [1.0, 6.0, 31.0, 201.0],
    },
    // Weeds — mugwort and ragweed are highly allergenic at low counts.
    PollenSpecies {
        api_key: "mugwort_pollen",
        display_name: "Mugwort",
        thresholds: [1.0, 6.0, 26.0, 51.0],
    },
    PollenSpecies {
        api_key: "olive_pollen",
        display_name: "Olive",
        thresholds: [1.0, 16.0, 51.0, 201.0],
    },
    PollenSpecies {
        api_key: "ragweed_pollen",
        display_name: "Ragweed",
        thresholds: [1.0, 6.0, 26.0, 101.0],
    },
];

impl PollenSpecies {
    /// Classify a concentration. Returns `None` for trace values that
    /// should be hidden from output.
    pub fn classify(&self, value: f64) -> Option<PollenSeverity> {
        if !value.is_finite() || value < self.thresholds[0] {
            return None;
        }
        if value >= self.thresholds[3] {
            Some(PollenSeverity::VeryHigh)
        } else if value >= self.thresholds[2] {
            Some(PollenSeverity::High)
        } else if value >= self.thresholds[1] {
            Some(PollenSeverity::Moderate)
        } else {
            Some(PollenSeverity::Low)
        }
    }
}

/// One pollen species worth showing to the user, after trace filtering.
pub struct PollenReading {
    pub display_name: &'static str,
    pub severity: PollenSeverity,
    pub value: f64,
}

// ===================== AQI TYPES =====================

/// European AQI category, per https://airindex.eea.europa.eu/.
pub fn european_aqi_category(value: f64) -> &'static str {
    // EEA bands: 0-20 Good, 20-40 Fair, 40-60 Moderate,
    // 60-80 Poor, 80-100 Very Poor, >100 Extremely Poor.
    match value as i32 {
        i32::MIN..=19 => "Good",
        20..=39 => "Fair",
        40..=59 => "Moderate",
        60..=79 => "Poor",
        80..=99 => "Very Poor",
        _ => "Extremely Poor",
    }
}

/// Descriptive category for Aerosol Optical Depth (550 nm, dimensionless).
///
/// Bands derived from NASA Earth Observatory and NOAA SURFRAD guidance:
/// AOD ≤ 0.1 = crystal-clear sky; 0.1–0.15 ≈ continental U.S. average;
/// ≥ 0.4 = very hazy; ≥ 1.0 = wildfire smoke / dust-storm territory.
/// We use the same five-tier vocabulary as the EAQI bands so the display
/// reads consistently.
pub fn aerosol_optical_depth_category(value: f64) -> &'static str {
    if !value.is_finite() || value < 0.1 {
        "Excellent"
    } else if value < 0.2 {
        "Good"
    } else if value < 0.3 {
        "Moderate"
    } else if value < 0.5 {
        "Poor"
    } else {
        "Very Poor"
    }
}

/// Dust below this concentration (μg/m³) is essentially background and
/// not worth flagging — we hide the line in that case to keep the
/// Argos panel from filling with noise on normal days.
pub const DUST_DISPLAY_THRESHOLD: f64 = 1.0;

// ===================== PER-POLLUTANT CATEGORIES =====================

/// Five-tier severity bands for an individual pollutant.
///
/// `thresholds[i]` is the *minimum* concentration (in the pollutant's
/// native unit) required to enter band `i+1`. Values strictly below
/// `thresholds[0]` map to "Good"; values ≥ `thresholds[3]` map to
/// "Very Poor". Non-finite inputs (NaN, ±∞) collapse to "Good".
struct PollutantBands {
    thresholds: [f64; 4],
}

impl PollutantBands {
    fn category(&self, value: f64) -> &'static str {
        if !value.is_finite() || value < self.thresholds[0] {
            "Good"
        } else if value < self.thresholds[1] {
            "Fair"
        } else if value < self.thresholds[2] {
            "Moderate"
        } else if value < self.thresholds[3] {
            "Poor"
        } else {
            "Very Poor"
        }
    }
}

// EAQI 1-hour sub-index breakpoints (μg/m³), per
// https://airindex.eea.europa.eu/AQI/index.html. The EAQI defines six
// bands (Good, Fair, Moderate, Poor, Very Poor, Extremely Poor); we
// collapse the top two into a single "Very Poor" tier so all per-pollutant
// categories share the same five-tier vocabulary as
// `aerosol_optical_depth_category`.
const PM25_BANDS: PollutantBands = PollutantBands { thresholds: [10.0, 20.0, 25.0, 50.0] };
const PM10_BANDS: PollutantBands = PollutantBands { thresholds: [20.0, 40.0, 50.0, 100.0] };
const NO2_BANDS: PollutantBands = PollutantBands { thresholds: [40.0, 90.0, 120.0, 230.0] };
const O3_BANDS: PollutantBands = PollutantBands { thresholds: [50.0, 100.0, 130.0, 240.0] };
const SO2_BANDS: PollutantBands = PollutantBands { thresholds: [100.0, 200.0, 350.0, 500.0] };

// CO has no EAQI sub-index. Bands derived from US EPA 8-hour AQI
// breakpoints (Good / Moderate / Unhealthy-for-Sensitive-Groups /
// Unhealthy / Very-Unhealthy+) converted to μg/m³ at standard conditions
// (1 ppm CO ≈ 1145 μg/m³ at 25 °C / 1013 hPa).
const CO_BANDS: PollutantBands = PollutantBands { thresholds: [4400.0, 9400.0, 12400.0, 15400.0] };

// Mineral dust has no EAQI sub-index. Bands reflect operational
// thresholds used in European Saharan-dust monitoring and WMO/MENA
// dust-storm classifications: tens of μg/m³ during routine intrusions,
// > 500 μg/m³ in severe storm events.
const DUST_BANDS: PollutantBands = PollutantBands { thresholds: [50.0, 150.0, 300.0, 500.0] };

/// Five-band category for PM2.5 (μg/m³), EAQI sub-index thresholds.
pub fn pm25_category(value: f64) -> &'static str {
    PM25_BANDS.category(value)
}

/// Five-band category for PM10 (μg/m³), EAQI sub-index thresholds.
pub fn pm10_category(value: f64) -> &'static str {
    PM10_BANDS.category(value)
}

/// Five-band category for NO2 (μg/m³), EAQI sub-index thresholds.
pub fn nitrogen_dioxide_category(value: f64) -> &'static str {
    NO2_BANDS.category(value)
}

/// Five-band category for O3 (μg/m³), EAQI sub-index thresholds.
pub fn ozone_category(value: f64) -> &'static str {
    O3_BANDS.category(value)
}

/// Five-band category for SO2 (μg/m³), EAQI sub-index thresholds.
pub fn sulphur_dioxide_category(value: f64) -> &'static str {
    SO2_BANDS.category(value)
}

/// Five-band category for CO (μg/m³), US EPA AQI 8-hour breakpoints.
pub fn carbon_monoxide_category(value: f64) -> &'static str {
    CO_BANDS.category(value)
}

/// Five-band category for mineral dust (μg/m³), operational dust-storm
/// thresholds.
pub fn dust_category(value: f64) -> &'static str {
    DUST_BANDS.category(value)
}

// ===================== API RESPONSE =====================

/// Subset of the Open-Meteo air-quality response we care about.
///
/// `elevation` is returned at the top level of every open-meteo
/// response (in metres above mean sea level) — we use it to back the
/// `--altitude auto` mode.
///
/// `current_units` carries the unit string that Open-Meteo associates
/// with each measurement. Display code consumes those verbatim instead
/// of hard-coding "μg/m³" so the output stays correct if the API ever
/// returns a different unit (e.g. ppb) for a given field.
///
/// `Option<f64>` everywhere on `current` — Open-Meteo omits fields the
/// caller didn't request, and we make this resilient to missing/null
/// values. `current_units` defaults to all-None for backward
/// compatibility with cache files written before this field existed.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AirQualityResponse {
    pub elevation: Option<f64>,
    #[serde(default)]
    pub current_units: CurrentUnits,
    pub current: CurrentValues,
}

/// Per-field unit strings reported by Open-Meteo (e.g. "μg/m³",
/// "grains/m³", "EAQI", or empty for dimensionless quantities like
/// aerosol optical depth).
#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct CurrentUnits {
    pub european_aqi: Option<String>,
    pub pm10: Option<String>,
    pub pm2_5: Option<String>,
    pub nitrogen_dioxide: Option<String>,
    pub carbon_monoxide: Option<String>,
    pub ozone: Option<String>,
    pub sulphur_dioxide: Option<String>,
    pub dust: Option<String>,
    pub aerosol_optical_depth: Option<String>,
    pub alder_pollen: Option<String>,
    pub birch_pollen: Option<String>,
    pub grass_pollen: Option<String>,
    pub mugwort_pollen: Option<String>,
    pub olive_pollen: Option<String>,
    pub ragweed_pollen: Option<String>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct CurrentValues {
    pub european_aqi: Option<f64>,
    pub pm10: Option<f64>,
    pub pm2_5: Option<f64>,
    pub nitrogen_dioxide: Option<f64>,
    pub carbon_monoxide: Option<f64>,
    pub ozone: Option<f64>,
    pub sulphur_dioxide: Option<f64>,
    pub dust: Option<f64>,
    pub aerosol_optical_depth: Option<f64>,
    pub alder_pollen: Option<f64>,
    pub birch_pollen: Option<f64>,
    pub grass_pollen: Option<f64>,
    pub mugwort_pollen: Option<f64>,
    pub olive_pollen: Option<f64>,
    pub ragweed_pollen: Option<f64>,
}

impl CurrentValues {
    /// Look up a species' concentration by its API key.
    fn pollen_value(&self, api_key: &str) -> Option<f64> {
        match api_key {
            "alder_pollen" => self.alder_pollen,
            "birch_pollen" => self.birch_pollen,
            "grass_pollen" => self.grass_pollen,
            "mugwort_pollen" => self.mugwort_pollen,
            "olive_pollen" => self.olive_pollen,
            "ragweed_pollen" => self.ragweed_pollen,
            _ => None,
        }
    }

    /// Build the list of pollen readings worth displaying (skipping
    /// trace/zero values per species' Low threshold).
    pub fn pollen_readings(&self) -> Vec<PollenReading> {
        POLLEN_SPECIES
            .iter()
            .filter_map(|sp| {
                let v = self.pollen_value(sp.api_key)?;
                let severity = sp.classify(v)?;
                Some(PollenReading { display_name: sp.display_name, severity, value: v })
            })
            .collect()
    }
}

// ===================== PUBLIC FETCH ENTRY =====================

/// Fetch air-quality / pollen / elevation data, using the on-disk
/// cache when fresh.
///
/// Returns `None` if both the cache is stale/absent *and* the network
/// fetch fails — callers should silently skip the section in that case
/// (see `print_argos` / `print_air_quality_terminal`). This keeps the
/// Argos panel clean when offline.
///
/// We always request the full AQI + pollen field set so the cache is
/// usable for any later flag combination within the same TTL window.
/// `elevation` comes back at the top level on every open-meteo
/// response, regardless of which `current` fields are requested.
///
/// `lat_str` and `lon_str` must already have been rounded to the
/// caller's chosen privacy level (see `cli::Coordinate::for_api`) —
/// this function takes pre-formatted strings rather than raw floats so
/// that the cache key and the URL stay byte-identical and we don't
/// silently re-introduce extra precision.
///
/// The caller is responsible for deciding whether a fetch is wanted at
/// all (i.e. only when `--aqi`, `--pollen`, or `--altitude auto` is
/// active) — this function unconditionally fetches when called.
pub fn fetch(lat_str: &str, lon_str: &str) -> Option<AirQualityResponse> {
    let cache = cache_path(lat_str, lon_str);

    // Try cache first — read_fresh_cache returns None for missing/stale/empty/invalid.
    if let Some(path) = cache.as_ref()
        && let Some(parsed) = read_fresh_cache(path)
    {
        return Some(parsed);
    }

    // Cache miss — fetch fresh.
    let url = build_url(lat_str, lon_str);
    let body = match http_get(&url) {
        Ok(b) => b,
        Err(_) => return None,
    };

    // Parse to validate before caching — refuse to write garbage.
    let parsed: AirQualityResponse = serde_json::from_str(&body).ok()?;

    // Best-effort cache write; failure is non-fatal (we still have the data).
    if let Some(path) = cache.as_ref() {
        let _ = write_cache_atomic(path, body.as_bytes());
    }

    Some(parsed)
}

// ===================== CACHE =====================

/// Resolve the cache file path for a given pre-rounded location.
///
/// Layout: `<XDG_RUNTIME_DIR>/suntimes/<lat>_<lon>.json`. The
/// subdirectory is created on demand by `write_cache_atomic`. Each
/// rounded location gets its own file, so changing `--latitude` /
/// `--longitude` (or `--location-privacy`) cannot return stale data
/// from a different place / privacy level.
///
/// Returns `None` on platforms without an XDG runtime dir (typical on
/// macOS/Windows) — caching is skipped silently in that case and every
/// invocation will re-fetch.
fn cache_path(lat_str: &str, lon_str: &str) -> Option<PathBuf> {
    let filename = format!("{}_{}.json", lat_str, lon_str);
    dirs::runtime_dir().map(|d| d.join(CACHE_SUBDIR).join(filename))
}

/// Read and parse the cache if it exists and is fresh.
///
/// Returns `None` (triggering a re-fetch) when:
/// * the file does not exist,
/// * its mtime is older than `CACHE_TTL`,
/// * it is empty,
/// * its contents fail to parse as `AirQualityResponse`.
fn read_fresh_cache(path: &Path) -> Option<AirQualityResponse> {
    let meta = std::fs::metadata(path).ok()?;
    let modified = meta.modified().ok()?;
    let age = SystemTime::now().duration_since(modified).ok()?;
    if age > CACHE_TTL {
        return None;
    }

    let s = std::fs::read_to_string(path).ok()?;
    if s.trim().is_empty() {
        return None;
    }

    serde_json::from_str(&s).ok()
}

/// Atomically write `bytes` to `target` via a same-directory temp file rename.
/// The rename is atomic on POSIX provided source and target share a filesystem;
/// placing the temp file alongside the target guarantees that.
fn write_cache_atomic(target: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let parent = target.parent().ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::InvalidInput, "cache path has no parent directory")
    })?;
    std::fs::create_dir_all(parent)?;

    let mut tmp = NamedTempFile::new_in(parent)?;
    tmp.write_all(bytes)?;
    tmp.flush()?;
    tmp.persist(target).map_err(|e| e.error)?;
    Ok(())
}

// ===================== HTTP =====================

/// Build the open-meteo URL. Coordinates are expected pre-rounded to
/// Build the open-meteo URL. Coordinate strings come pre-formatted to
/// the caller's chosen precision (`cli::Coordinate::for_api` →
/// privacy-rounded, no padding). We pass them through verbatim.
///
/// All AQI and pollen fields are always requested, regardless of which
/// CLI flag the user set. This guarantees the cache file is complete
/// so that a later invocation with a different flag combination still
/// has the data it needs.
///
/// We pin `timezone=UTC` regardless of the user's `--timezone` setting:
/// we only consume the `current` block (numeric measurements that don't
/// vary with timezone) and never display the API's `time` string, so
/// the parameter is decorative for our use case. Pinning it stabilises
/// the URL — switching `--timezone` between invocations does not cause
/// extra fetches and the cache stays usable.
fn build_url(lat_str: &str, lon_str: &str) -> String {
    let mut params: Vec<&str> = Vec::with_capacity(AQI_FIELDS.len() + POLLEN_SPECIES.len());
    params.extend_from_slice(AQI_FIELDS);
    for sp in POLLEN_SPECIES {
        params.push(sp.api_key);
    }

    format!(
        "{}?latitude={}&longitude={}&current={}&timezone=UTC",
        API_BASE,
        lat_str,
        lon_str,
        params.join(","),
    )
}

/// Issue the GET, transparently honouring `HTTP_PROXY` / `HTTPS_PROXY`
/// / `ALL_PROXY` from the environment.
fn http_get(url: &str) -> Result<String, Box<dyn std::error::Error>> {
    let user_agent = concat!("suntimes-rs/", env!("CARGO_PKG_VERSION"));

    let mut builder =
        ureq::Agent::config_builder().timeout_global(Some(HTTP_TIMEOUT)).user_agent(user_agent);

    // Pick up proxy from environment if set (ALL_PROXY / HTTPS_PROXY / HTTP_PROXY).
    if let Some(proxy) = ureq::Proxy::try_from_env() {
        builder = builder.proxy(Some(proxy));
    }

    let agent: ureq::Agent = builder.build().into();

    let body = agent.get(url).call()?.body_mut().read_to_string()?;
    Ok(body)
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_european_aqi_category_bands() {
        assert_eq!(european_aqi_category(0.0), "Good");
        assert_eq!(european_aqi_category(19.9), "Good");
        assert_eq!(european_aqi_category(20.0), "Fair");
        assert_eq!(european_aqi_category(39.0), "Fair");
        assert_eq!(european_aqi_category(40.0), "Moderate");
        assert_eq!(european_aqi_category(60.0), "Poor");
        assert_eq!(european_aqi_category(80.0), "Very Poor");
        assert_eq!(european_aqi_category(100.0), "Extremely Poor");
        assert_eq!(european_aqi_category(250.0), "Extremely Poor");
    }

    #[test]
    fn test_aerosol_optical_depth_category_bands() {
        // NASA Earthdata / NOAA SURFRAD reference points
        assert_eq!(aerosol_optical_depth_category(0.0), "Excellent");
        assert_eq!(aerosol_optical_depth_category(0.05), "Excellent"); // very clean
        assert_eq!(aerosol_optical_depth_category(0.099), "Excellent");
        assert_eq!(aerosol_optical_depth_category(0.1), "Good"); // continental average
        assert_eq!(aerosol_optical_depth_category(0.15), "Good");
        assert_eq!(aerosol_optical_depth_category(0.2), "Moderate");
        assert_eq!(aerosol_optical_depth_category(0.3), "Poor");
        assert_eq!(aerosol_optical_depth_category(0.4), "Poor"); // very hazy per NOAA
        assert_eq!(aerosol_optical_depth_category(0.5), "Very Poor");
        assert_eq!(aerosol_optical_depth_category(1.0), "Very Poor"); // wildfire/dust storm
        assert_eq!(aerosol_optical_depth_category(4.0), "Very Poor");

        // NaN / weird inputs treated as Excellent rather than panicking.
        assert_eq!(aerosol_optical_depth_category(f64::NAN), "Excellent");
        assert_eq!(aerosol_optical_depth_category(-0.1), "Excellent");
    }

    #[test]
    fn test_pm25_category_bands() {
        // EAQI 1-hour PM2.5: <10 Good, <20 Fair, <25 Moderate, <50 Poor, ≥50 Very Poor.
        assert_eq!(pm25_category(0.0), "Good");
        assert_eq!(pm25_category(9.9), "Good");
        assert_eq!(pm25_category(10.0), "Fair");
        assert_eq!(pm25_category(19.9), "Fair");
        assert_eq!(pm25_category(20.0), "Moderate");
        assert_eq!(pm25_category(24.9), "Moderate");
        assert_eq!(pm25_category(25.0), "Poor");
        assert_eq!(pm25_category(49.9), "Poor");
        assert_eq!(pm25_category(50.0), "Very Poor");
        assert_eq!(pm25_category(500.0), "Very Poor");
        // Non-finite / negative degrade to "Good".
        assert_eq!(pm25_category(f64::NAN), "Good");
        assert_eq!(pm25_category(-1.0), "Good");
    }

    #[test]
    fn test_pm10_category_bands() {
        // EAQI 1-hour PM10: <20 Good, <40 Fair, <50 Moderate, <100 Poor, ≥100 Very Poor.
        assert_eq!(pm10_category(0.0), "Good");
        assert_eq!(pm10_category(20.0), "Fair");
        assert_eq!(pm10_category(40.0), "Moderate");
        assert_eq!(pm10_category(50.0), "Poor");
        assert_eq!(pm10_category(100.0), "Very Poor");
        assert_eq!(pm10_category(500.0), "Very Poor");
    }

    #[test]
    fn test_nitrogen_dioxide_category_bands() {
        // EAQI 1-hour NO2: <40, <90, <120, <230.
        assert_eq!(nitrogen_dioxide_category(0.0), "Good");
        assert_eq!(nitrogen_dioxide_category(40.0), "Fair");
        assert_eq!(nitrogen_dioxide_category(90.0), "Moderate");
        assert_eq!(nitrogen_dioxide_category(120.0), "Poor");
        assert_eq!(nitrogen_dioxide_category(230.0), "Very Poor");
    }

    #[test]
    fn test_ozone_category_bands() {
        // EAQI 1-hour O3: <50, <100, <130, <240.
        assert_eq!(ozone_category(0.0), "Good");
        assert_eq!(ozone_category(50.0), "Fair");
        assert_eq!(ozone_category(100.0), "Moderate");
        assert_eq!(ozone_category(130.0), "Poor");
        assert_eq!(ozone_category(240.0), "Very Poor");
        assert_eq!(ozone_category(67.0), "Fair");
    }

    #[test]
    fn test_sulphur_dioxide_category_bands() {
        // EAQI 1-hour SO2: <100, <200, <350, <500.
        assert_eq!(sulphur_dioxide_category(0.0), "Good");
        assert_eq!(sulphur_dioxide_category(100.0), "Fair");
        assert_eq!(sulphur_dioxide_category(200.0), "Moderate");
        assert_eq!(sulphur_dioxide_category(350.0), "Poor");
        assert_eq!(sulphur_dioxide_category(500.0), "Very Poor");
    }

    #[test]
    fn test_carbon_monoxide_category_bands() {
        // US EPA AQI 8-hour CO breakpoints, μg/m³.
        assert_eq!(carbon_monoxide_category(0.0), "Good");
        assert_eq!(carbon_monoxide_category(147.0), "Good");
        assert_eq!(carbon_monoxide_category(4400.0), "Fair");
        assert_eq!(carbon_monoxide_category(9400.0), "Moderate");
        assert_eq!(carbon_monoxide_category(12400.0), "Poor");
        assert_eq!(carbon_monoxide_category(15400.0), "Very Poor");
    }

    #[test]
    fn test_dust_category_bands() {
        // Operational dust-storm thresholds.
        assert_eq!(dust_category(0.0), "Good");
        assert_eq!(dust_category(50.0), "Fair");
        assert_eq!(dust_category(150.0), "Moderate");
        assert_eq!(dust_category(300.0), "Poor");
        assert_eq!(dust_category(500.0), "Very Poor");
    }

    #[test]
    fn test_current_units_default_when_missing() {
        // Backward compat: cache files written before current_units was
        // added should still deserialize cleanly.
        let json = r#"{"elevation": 18.0, "current": {"european_aqi": 27}}"#;
        let parsed: AirQualityResponse = serde_json::from_str(json).unwrap();
        assert_eq!(parsed.elevation, Some(18.0));
        assert!(parsed.current_units.pm2_5.is_none());
        assert!(parsed.current_units.european_aqi.is_none());
    }

    #[test]
    fn test_current_units_parsed_when_present() {
        let json = r#"{
            "elevation": 11.0,
            "current_units": {
                "european_aqi": "EAQI",
                "pm2_5": "μg/m³",
                "aerosol_optical_depth": ""
            },
            "current": {"european_aqi": 27, "pm2_5": 2.6}
        }"#;
        let parsed: AirQualityResponse = serde_json::from_str(json).unwrap();
        assert_eq!(parsed.current_units.european_aqi.as_deref(), Some("EAQI"));
        assert_eq!(parsed.current_units.pm2_5.as_deref(), Some("μg/m³"));
        // Empty unit string is preserved (AOD is dimensionless in the API).
        assert_eq!(parsed.current_units.aerosol_optical_depth.as_deref(), Some(""));
    }

    fn species(name: &str) -> &'static PollenSpecies {
        POLLEN_SPECIES.iter().find(|s| s.api_key == name).expect("species in table")
    }

    #[test]
    fn test_pollen_trace_hidden() {
        // Anything below the species' Low threshold is hidden.
        assert!(species("birch_pollen").classify(0.0).is_none());
        assert!(species("birch_pollen").classify(0.5).is_none());
        assert!(species("birch_pollen").classify(0.999).is_none());

        // Even strictly negative or NaN inputs should not crash.
        assert!(species("birch_pollen").classify(-1.0).is_none());
        assert!(species("birch_pollen").classify(f64::NAN).is_none());
    }

    #[test]
    fn test_pollen_birch_bands() {
        let s = species("birch_pollen");
        assert_eq!(s.classify(1.0), Some(PollenSeverity::Low));
        assert_eq!(s.classify(10.0), Some(PollenSeverity::Low));
        assert_eq!(s.classify(11.0), Some(PollenSeverity::Moderate));
        assert_eq!(s.classify(100.0), Some(PollenSeverity::Moderate));
        assert_eq!(s.classify(101.0), Some(PollenSeverity::High));
        assert_eq!(s.classify(1000.0), Some(PollenSeverity::High));
        assert_eq!(s.classify(1001.0), Some(PollenSeverity::VeryHigh));
        assert_eq!(s.classify(2533.2), Some(PollenSeverity::VeryHigh));
    }

    #[test]
    fn test_pollen_grass_bands_lower_thresholds() {
        // Grass has much tighter thresholds than birch — a value of 6 is
        // already Moderate for grass but only Low for birch.
        let grass = species("grass_pollen");
        let birch = species("birch_pollen");
        assert_eq!(grass.classify(6.0), Some(PollenSeverity::Moderate));
        assert_eq!(birch.classify(6.0), Some(PollenSeverity::Low));
    }

    #[test]
    fn test_pollen_readings_filters_zero_and_trace() {
        let v = CurrentValues {
            birch_pollen: Some(2533.2),
            grass_pollen: Some(0.0),
            mugwort_pollen: Some(0.0),
            olive_pollen: Some(0.0),
            ragweed_pollen: Some(0.0),
            alder_pollen: Some(0.7), // trace, hidden
            ..Default::default()
        };
        let readings = v.pollen_readings();
        assert_eq!(readings.len(), 1);
        assert_eq!(readings[0].display_name, "Birch");
        assert_eq!(readings[0].severity, PollenSeverity::VeryHigh);
    }

    #[test]
    fn test_url_always_includes_all_fields() {
        // build_url takes pre-formatted coordinate strings now —
        // privacy rounding lives in cli::Coordinate::for_api.
        let url = build_url("60.4", "25.3");
        assert!(url.starts_with("https://air-quality-api.open-meteo.com/v1/air-quality?"));
        assert!(url.contains("latitude=60.4"));
        assert!(url.contains("longitude=25.3"));
        // build_url passes its inputs through unchanged — no padding here either.
        assert!(!url.contains("latitude=60.40"));

        // Caller-supplied precision is honoured verbatim:
        let url_low = build_url("60.382", "25.273");
        assert!(url_low.contains("latitude=60.382"));
        assert!(url_low.contains("longitude=25.273"));

        // All AQI fields present (concentrations + AOD + overall index)
        for f in AQI_FIELDS {
            assert!(url.contains(f), "missing AQI field {} in URL: {}", f, url);
        }
        assert!(url.contains("sulphur_dioxide"));
        assert!(url.contains("dust"));
        assert!(url.contains("aerosol_optical_depth"));
        // The SO2 *sub-index* is intentionally NOT requested anymore —
        // we display SO2 as raw concentration like the other pollutants.
        assert!(!url.contains("european_aqi_sulphur_dioxide"));

        // All pollen species present
        for sp in POLLEN_SPECIES {
            assert!(url.contains(sp.api_key), "missing {} in URL: {}", sp.api_key, url);
        }

        assert!(url.contains("timezone=UTC"));
    }

    #[test]
    fn test_cache_path_filename_per_location() {
        // We can't assert the parent directory (depends on the test
        // runner's environment), but we can check the filename format.
        // Skip the assertion when XDG_RUNTIME_DIR is not set on the host.
        if let Some(path) = cache_path("60.4", "25.3") {
            let name = path.file_name().unwrap().to_string_lossy();
            assert_eq!(name, "60.4_25.3.json");
            // Subdirectory should be `suntimes`.
            let parent = path.parent().unwrap().file_name().unwrap().to_string_lossy();
            assert_eq!(parent, CACHE_SUBDIR);
        }

        if let Some(path) = cache_path("-60.4", "-25.3") {
            let name = path.file_name().unwrap().to_string_lossy();
            assert_eq!(name, "-60.4_-25.3.json");
        }

        // Different precisions = different cache files. This is by
        // design: switching --location-privacy after a fetch leaves the
        // old cache file behind and starts a new one for the new precision.
        if let Some(path) = cache_path("60.38", "25.27") {
            let name = path.file_name().unwrap().to_string_lossy();
            assert_eq!(name, "60.38_25.27.json");
        }
    }

    #[test]
    fn test_atomic_write_replaces_existing() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("cache.json");
        std::fs::write(&target, b"old").unwrap();
        write_cache_atomic(&target, b"{\"current\":{}}").unwrap();
        let s = std::fs::read_to_string(&target).unwrap();
        assert_eq!(s, "{\"current\":{}}");
    }

    #[test]
    fn test_cache_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("cache.json");
        let payload = serde_json::json!({
            "elevation": 18.0,
            "current": {
                "european_aqi": 34,
                "pm2_5": 3.9,
                "birch_pollen": 2533.2
            }
        });
        let bytes = serde_json::to_vec(&payload).unwrap();
        write_cache_atomic(&target, &bytes).unwrap();

        let parsed = read_fresh_cache(&target).expect("cache should be fresh and parseable");
        assert_eq!(parsed.elevation, Some(18.0));
        assert_eq!(parsed.current.european_aqi, Some(34.0));
        assert_eq!(parsed.current.pm2_5, Some(3.9));
        assert_eq!(parsed.current.birch_pollen, Some(2533.2));
    }

    #[test]
    fn test_cache_invalid_json_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("cache.json");
        std::fs::write(&target, b"not json").unwrap();
        assert!(read_fresh_cache(&target).is_none());
    }

    #[test]
    fn test_cache_empty_file_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("cache.json");
        std::fs::write(&target, b"").unwrap();
        assert!(read_fresh_cache(&target).is_none());
    }
}
