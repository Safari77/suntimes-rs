//! Command-Line Interface Module
//!
//! Handles argument parsing and validation for the suntimes-rs application.

use clap::{Parser, ValueEnum};
use serde::Deserialize;

// ===================== TYPES =====================

/// Altitude argument: either a literal number (metres above MSL) or
/// the keyword `auto`, which makes main.rs pull the value from
/// open-meteo's `elevation` field at runtime.
///
/// `auto` is silently downgraded to 0.0 if the open-meteo fetch fails
/// — by design, per CLI contract: a transient network outage should
/// not turn a valid invocation into an error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AltitudeArg {
    Auto,
    Fixed(f64),
}

impl std::fmt::Display for AltitudeArg {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AltitudeArg::Auto => f.write_str("auto"),
            AltitudeArg::Fixed(v) => write!(f, "{}", v),
        }
    }
}

/// Privacy level controlling how coarsely we round latitude/longitude
/// before they leave the machine (sent to open-meteo and stored as
/// part of the cache filename).
///
/// Rounding precision (per coordinate axis):
/// * `high`   → 0.1°  ≈ 11 km
/// * `medium` → 0.01° ≈ 1.1 km
/// * `low`    → 0.001° ≈ 110 m  (useful for mountain UV / `--altitude auto`)
///
/// Solar calculations always use the user's full input precision —
/// this only affects what's sent over the network and the cache key.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocationPrivacy {
    High,
    Medium,
    Low,
}

impl LocationPrivacy {
    /// Maximum number of decimal places the privacy level permits in
    /// the outbound coordinate.
    pub fn max_decimal_places(&self) -> u8 {
        match self {
            LocationPrivacy::High => 1,
            LocationPrivacy::Medium => 2,
            LocationPrivacy::Low => 3,
        }
    }
}

impl std::fmt::Display for LocationPrivacy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Match clap's lowercased ValueEnum names so default_value_t
        // round-trips correctly through the parser.
        let s = match self {
            LocationPrivacy::High => "high",
            LocationPrivacy::Medium => "medium",
            LocationPrivacy::Low => "low",
        };
        f.write_str(s)
    }
}

/// A coordinate captured from the CLI, retaining the user's original
/// input precision
#[derive(Clone, Copy, Debug)]
pub struct Coordinate {
    /// Numeric value as parsed.
    pub value: f64,
    /// Number of decimal places the user typed. `0` if they typed an
    /// integer or used scientific notation without a fractional part.
    pub user_decimal_places: u8,
}

impl Coordinate {
    /// Format the coordinate for outbound use (URL + cache filename),
    /// applying the *coarser* of (privacy ceiling, user precision).
    pub fn for_api(&self, privacy: LocationPrivacy) -> String {
        let dp = privacy.max_decimal_places().min(self.user_decimal_places) as usize;
        // Round the f64 to that precision before formatting so that
        // float jitter never bleeds extra (false-precision) digits in.
        let factor = 10f64.powi(dp as i32);
        let rounded = ((self.value * factor).round() / factor) + 0.0;
        format!("{:.*}", dp, rounded)
    }
}

/// Count the decimal places in a numeric string. Counts only the
/// leading run of digits after the decimal point — anything after
/// (e.g. an `e`-style exponent) is ignored, so we don't mis-count
/// "55.12e1" as 4 dp.
fn count_decimal_places(s: &str) -> u8 {
    let frac = match s.split_once('.') {
        Some((_, f)) => f,
        None => return 0,
    };
    // Cap at 10 so this never overflows even for pathological inputs.
    frac.chars().take_while(|c| c.is_ascii_digit()).count().min(10) as u8
}

/// Sun position model.
///
/// Clap renders each variant lowercased in `--help`, preserving the original
/// string interface (`--model noaa | horizons | physical`).
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum SunModel {
    /// NOAA-style calculation (refraction + horizon dip for altitude)
    Noaa,
    /// Geometric horizon (refraction on, no altitude dip — matches JPL HORIZONS)
    Horizons,
    /// Physical model (altitude-adjusted pressure + temperature refraction)
    Physical,
}

impl std::fmt::Display for SunModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Must match clap's lowercased variant names so `default_value_t`
        // round-trips correctly through the ValueEnum parser.
        let s = match self {
            SunModel::Noaa => "noaa",
            SunModel::Horizons => "horizons",
            SunModel::Physical => "physical",
        };
        f.write_str(s)
    }
}

/// Individual pollutants that the `--aqi-display` filter can select.
///
/// Variant `value` names match the Open-Meteo air-quality API field
/// names verbatim (`pm2_5`, `nitrogen_dioxide`, …) so the CLI vocabulary
/// is the same as the API documentation. `european_aqi` and
/// `aerosol_optical_depth` are intentionally NOT in this enum: they are
/// always shown when `--aqi` is set, regardless of `--aqi-display`.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum AqiPollutant {
    #[value(name = "pm10")]
    Pm10,
    #[value(name = "pm2_5")]
    Pm25,
    #[value(name = "nitrogen_dioxide")]
    NitrogenDioxide,
    #[value(name = "carbon_monoxide")]
    CarbonMonoxide,
    #[value(name = "ozone")]
    Ozone,
    #[value(name = "sulphur_dioxide")]
    SulphurDioxide,
    #[value(name = "dust")]
    Dust,
}

// ===================== CLI =====================

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Args {
    /// Observer latitude in decimal degrees (-90 to 90)
    #[arg(
        long,
        allow_hyphen_values = true,
        value_parser = parse_latitude,
        env = "ARGOS_SUNTIMES_LATITUDE",
    )]
    pub latitude: Coordinate,
    /// Observer longitude in decimal degrees (-180 to 180)
    #[arg(
        long,
        allow_hyphen_values = true,
        value_parser = parse_longitude,
        env = "ARGOS_SUNTIMES_LONGITUDE",
    )]
    pub longitude: Coordinate,
    /// Time zone to use ("system", "location", or IANA time zone name)
    #[arg(long, default_value = "system", env = "ARGOS_SUNTIMES_TIMEZONE")]
    pub timezone: String,
    /// How coarsely to round coordinates before sending to external APIs
    /// (and using as the cache key). Solar/UV calculations always use
    /// the full input precision regardless of this setting.
    #[arg(long, default_value_t = LocationPrivacy::High, value_enum, env = "ARGOS_SUNTIMES_LOCATION_PRIVACY")]
    pub location_privacy: LocationPrivacy,
    /// Observer altitude above mean sea level (metres, may be negative)
    /// or the keyword "auto" to pull elevation from open-meteo.
    /// Numeric range: -500m (Dead Sea) to 11000m (Troposphere limit for ISA formula).
    /// "auto" forces an open-meteo fetch even without --aqi/--pollen,
    /// and silently falls back to 0.0 if the fetch fails.
    #[arg(long, default_value_t = AltitudeArg::Fixed(0.0), allow_hyphen_values = true, value_parser = parse_altitude, env = "ARGOS_SUNTIMES_ALTITUDE")]
    pub altitude: AltitudeArg,

    /// Sun position model to use
    #[arg(long, default_value_t = SunModel::Noaa, value_enum, env = "ARGOS_SUNTIMES_MODEL")]
    pub model: SunModel,
    /// Ambient temperature in °C for refraction correction (physical model only)
    #[arg(long, default_value_t = 15.0, allow_hyphen_values = true)]
    pub temperature: f64,
    /// Atmospheric pressure in hPa for refraction correction (physical model only)
    #[arg(long, default_value_t = 1013.25)]
    pub pressure: f64,

    /// Output in Argos (GNOME Shell) format
    #[arg(long, env = "ARGOS_SUNTIMES_ENABLE")]
    pub argos: bool,
    // This hidden field defines the negation flag
    #[arg(long, hide = true, action = clap::ArgAction::SetTrue)]
    pub no_argos: bool,

    /// Use civil, nautical, or astronomical twilight instead of sunrise/sunset
    #[arg(long)]
    pub civil: bool,
    #[arg(long)]
    pub nautical: bool,
    #[arg(long)]
    pub astro: bool,

    /// Date for calculations (e.g., "2024-12-25" or "today"); defaults to today
    #[arg(long)]
    pub date: Option<String>,
    /// Show Sun position at a specific time (HH:MM[:SS[.fffffffff]] or "now")
    #[arg(long)]
    pub at: Option<String>,
    /// Use UTC time zone
    #[arg(long)]
    pub utc: bool,
    /// Use the old empirical formula-based model for UV Index calculations
    #[arg(long)]
    pub formula_calcs: bool,

    /// Show build info from Cargo.lock at time of building
    #[arg(long)]
    pub show_build_info: bool,

    // ===================== AIR QUALITY & POLLEN =====================
    /// Enable Air Quality Index (European AQI) and core pollutants display
    #[arg(long, env = "ARGOS_SUNTIMES_AQI")]
    pub aqi: bool,

    /// Include regional pollen forecast (e.g., birch, grass) in output
    #[arg(long, env = "ARGOS_SUNTIMES_POLLEN")]
    pub pollen: bool,

    /// Pollutants to include in the air-quality display. Comma-separated subset
    /// of: pm10, pm2_5, nitrogen_dioxide, carbon_monoxide, ozone,
    /// sulphur_dioxide, dust. If omitted (the default), all are displayed.
    /// `european_aqi` and `aerosol_optical_depth` are always shown when --aqi
    /// is set. The display order is fixed by the program; this flag only
    /// filters. All values are always fetched and cached regardless of this
    /// flag, so flipping it between invocations does not trigger extra API
    /// calls.
    #[arg(long, value_delimiter = ',', env = "ARGOS_SUNTIMES_AQI_DISPLAY")]
    pub aqi_display: Vec<AqiPollutant>,

    // ===================== SOLAR PANEL OPTIONS =====================
    /// Solar panel area in square meters (enables solar output calculation)
    #[arg(long, value_parser = parse_positive_f64, env = "ARGOS_SUNTIMES_SOLARPANEL_SIZE")]
    pub solarpanel_size: Option<f64>,

    /// Solar panel tilt angle in degrees (0 = flat/horizontal, 90 = vertical)
    #[arg(long, default_value_t = 35.0, value_parser = parse_tilt, env = "ARGOS_SUNTIMES_SOLARPANEL_TILT")]
    pub solarpanel_tilt: f64,

    /// Solar panel azimuth in degrees (180 = facing south in northern hemisphere)
    #[arg(long, default_value_t = 180.0, value_parser = parse_azimuth, env = "ARGOS_SUNTIMES_SOLARPANEL_AZIMUTH")]
    pub solarpanel_azimuth: f64,

    /// Solar panel efficiency (0.0-1.0, typical ~0.18-0.22 for silicon)
    #[arg(long, default_value_t = 0.20, value_parser = parse_efficiency, env = "ARGOS_SUNTIMES_SOLARPANEL_EFFICIENCY")]
    pub solarpanel_efficiency: f64,

    /// Linke turbidity factor for clear-sky model (2-7 typical, 3 = clear)
    #[arg(long, default_value_t = 3.0, value_parser = parse_turbidity, env = "ARGOS_SUNTIMES_LINKE_TURBIDITY")]
    pub linke_turbidity: f64,

    /// Ground albedo for reflected radiation (0.0-1.0, 0.2 = grass, 0.8 = snow)
    #[arg(long, default_value_t = 0.2, value_parser = parse_albedo, env = "ARGOS_SUNTIMES_ALBEDO")]
    pub albedo: f64,

    /// Enable dual-axis tracking (panel always faces the sun)
    /// Provides +30-40% energy vs fixed, highest cost
    #[arg(long, conflicts_with_all = ["solarpanel_horizontal_tracking", "solarpanel_vertical_tracking"])]
    pub solarpanel_dual_axis: bool,

    /// Enable horizontal single-axis tracking (HSAT)
    /// Axis runs N-S, panel rotates E-W to track sun's azimuth
    /// Tilt stays fixed (use --solarpanel-tilt or --solarpanel-find-optimum)
    #[arg(long, conflicts_with_all = ["solarpanel_dual_axis", "solarpanel_vertical_tracking"])]
    pub solarpanel_horizontal_tracking: bool,

    /// Enable vertical single-axis tracking (VSAT)
    /// Axis is vertical, panel tilts to track sun's altitude
    /// Azimuth stays fixed (use --solarpanel-azimuth or --solarpanel-find-optimum)
    #[arg(long, conflicts_with_all = ["solarpanel_dual_axis", "solarpanel_horizontal_tracking"])]
    pub solarpanel_vertical_tracking: bool,

    /// Find optimum panel configuration for maximum daily energy
    /// - Fixed panel: finds optimal tilt and azimuth
    /// - HSAT: finds optimal tilt (azimuth is tracked)
    /// - VSAT: finds optimal azimuth (tilt is tracked)
    ///   Not compatible with dual-axis tracking
    #[arg(long, conflicts_with = "solarpanel_dual_axis")]
    pub solarpanel_find_optimum: bool,

    /// Find optimal yearly adjustment schedule for N adjustments (1-366)
    /// Outputs dates and tilt settings to maximize yearly kWh
    /// Valid for fixed panels and HSAT (not compatible with VSAT or dual-axis)
    #[arg(long, value_parser = parse_adjustments, env = "ARGOS_SUNTIMES_SOLARPANEL_YEARLY_ADJUSTMENTS",
          conflicts_with_all = ["solarpanel_dual_axis", "solarpanel_vertical_tracking"])]
    pub solarpanel_yearly_adjustments: Option<u32>,

    /// Tilt range constraint: "MIN-MAX" (e.g., "20-60" limits tilt to 20°-60°)
    /// If not specified, full range 0-90 is used
    #[arg(long, value_parser = parse_range, env = "ARGOS_SUNTIMES_SOLARPANEL_TILT_RANGE")]
    pub solarpanel_tilt_range: Option<(f64, f64)>,

    /// Azimuth range constraint: "MIN-MAX" (e.g., "150-210" limits azimuth to 150°-210°)
    /// If not specified, full range 0-360 is used
    #[arg(long, value_parser = parse_range, env = "ARGOS_SUNTIMES_SOLARPANEL_AZIMUTH_RANGE")]
    pub solarpanel_azimuth_range: Option<(f64, f64)>,
}

// Define the structure to match what we serialized in build.rs
#[derive(Debug, Deserialize)]
pub struct DepInfo {
    pub name: String,
    pub version: String,
    pub checksum: Option<String>,
    pub source: Option<String>,
}

// ===================== CLI VALUE PARSERS =====================

fn parse_latitude(s: &str) -> Result<Coordinate, String> {
    let value: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-90.0..=90.0).contains(&value) {
        return Err(format!("Latitude must be between -90 and 90, got {}", value));
    }
    Ok(Coordinate { value, user_decimal_places: count_decimal_places(s) })
}

fn parse_longitude(s: &str) -> Result<Coordinate, String> {
    let value: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-180.0..=180.0).contains(&value) {
        return Err(format!("Longitude must be between -180 and 180, got {}", value));
    }
    Ok(Coordinate { value, user_decimal_places: count_decimal_places(s) })
}

fn parse_altitude(s: &str) -> Result<AltitudeArg, String> {
    if s.eq_ignore_ascii_case("auto") {
        return Ok(AltitudeArg::Auto);
    }
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-500.0..=11000.0).contains(&v) {
        return Err(format!("Altitude must be between -500 and 11000 meters, got {}", v));
    }
    Ok(AltitudeArg::Fixed(v))
}

fn parse_positive_f64(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if v <= 0.0 {
        return Err(format!("Value must be positive, got {}", v));
    }
    Ok(v)
}

fn parse_tilt(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(0.0..=90.0).contains(&v) {
        return Err(format!("Tilt must be between 0 and 90 degrees, got {}", v));
    }
    Ok(v)
}

fn parse_azimuth(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(0.0..=360.0).contains(&v) {
        return Err(format!("Azimuth must be between 0 and 360 degrees, got {}", v));
    }
    Ok(v)
}

fn parse_efficiency(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(0.0..=1.0).contains(&v) {
        return Err(format!("Efficiency must be between 0.0 and 1.0, got {}", v));
    }
    Ok(v)
}

fn parse_turbidity(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(1.0..=10.0).contains(&v) {
        return Err(format!("Linke turbidity must be between 1.0 and 10.0, got {}", v));
    }
    Ok(v)
}

fn parse_albedo(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(0.0..=1.0).contains(&v) {
        return Err(format!("Albedo must be between 0.0 and 1.0, got {}", v));
    }
    Ok(v)
}

fn parse_adjustments(s: &str) -> Result<u32, String> {
    let v: u32 = s.parse().map_err(|_| format!("Invalid integer: {}", s))?;
    if !(1..=366).contains(&v) {
        return Err(format!("Number of adjustments must be between 1 and 366, got {}", v));
    }
    Ok(v)
}

fn parse_range(s: &str) -> Result<(f64, f64), String> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 2 {
        return Err(format!("Range must be in format MIN-MAX (e.g., '20-60'), got '{}'", s));
    }
    let min: f64 = parts[0].parse().map_err(|_| format!("Invalid minimum value: {}", parts[0]))?;
    let max: f64 = parts[1].parse().map_err(|_| format!("Invalid maximum value: {}", parts[1]))?;
    if min > max {
        return Err(format!("Minimum ({}) cannot be greater than maximum ({})", min, max));
    }
    Ok((min, max))
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_decimal_places() {
        assert_eq!(count_decimal_places("55"), 0);
        assert_eq!(count_decimal_places("55."), 0);
        assert_eq!(count_decimal_places("55.1"), 1);
        assert_eq!(count_decimal_places("55.12"), 2);
        assert_eq!(count_decimal_places("55.123"), 3);
        assert_eq!(count_decimal_places("55.123456"), 6);
        assert_eq!(count_decimal_places("-55.12"), 2);
        // Trailing zeros are preserved (literal precision).
        assert_eq!(count_decimal_places("55.10"), 2);
        assert_eq!(count_decimal_places("55.100"), 3);
        // Scientific notation: only the leading run of digits after
        // the dot counts; the exponent portion is ignored.
        assert_eq!(count_decimal_places("55.12e1"), 2);
        // Cap at 10 to avoid surprises with absurdly long inputs.
        assert_eq!(count_decimal_places("0.12345678901234567890"), 10);
    }

    #[test]
    fn test_privacy_decimal_places() {
        assert_eq!(LocationPrivacy::High.max_decimal_places(), 1);
        assert_eq!(LocationPrivacy::Medium.max_decimal_places(), 2);
        assert_eq!(LocationPrivacy::Low.max_decimal_places(), 3);
    }

    /// Build a Coordinate as if parsed from the given string.
    fn coord(s: &str) -> Coordinate {
        Coordinate { value: s.parse().unwrap(), user_decimal_places: count_decimal_places(s) }
    }

    #[test]
    fn test_for_api_user_coarser_than_privacy() {
        // The example from the spec: --location-privacy low --latitude 55.12
        // → user gave 2 dp, privacy permits 3, output stays at 2 dp.
        assert_eq!(coord("55.12").for_api(LocationPrivacy::Low), "55.12");
        assert_eq!(coord("55.1").for_api(LocationPrivacy::Low), "55.1");
        assert_eq!(coord("55").for_api(LocationPrivacy::Low), "55");
        // Same input, the other privacy levels — user precision still wins
        // because it's the coarser side.
        assert_eq!(coord("55.12").for_api(LocationPrivacy::Medium), "55.12");
        assert_eq!(coord("55.1").for_api(LocationPrivacy::Medium), "55.1");
        assert_eq!(coord("55.1").for_api(LocationPrivacy::High), "55.1");
    }

    #[test]
    fn test_for_api_privacy_coarser_than_user() {
        // User typed more precision than the privacy level allows —
        // the privacy ceiling kicks in and we round.
        assert_eq!(coord("60.382517").for_api(LocationPrivacy::High), "60.4");
        assert_eq!(coord("60.382517").for_api(LocationPrivacy::Medium), "60.38");
        assert_eq!(coord("60.382517").for_api(LocationPrivacy::Low), "60.383");

        // Negative coordinates work the same way.
        assert_eq!(coord("-25.273881").for_api(LocationPrivacy::High), "-25.3");
        assert_eq!(coord("-25.273881").for_api(LocationPrivacy::Medium), "-25.27");
        assert_eq!(coord("-25.273881").for_api(LocationPrivacy::Low), "-25.274");
    }

    #[test]
    fn test_for_api_no_padding() {
        // Critical: don't pad. "55.12" with --location-privacy low
        // must NOT become "55.120".
        let s = coord("55.12").for_api(LocationPrivacy::Low);
        assert_eq!(s, "55.12");
        assert!(!s.contains("55.120"));
    }

    #[test]
    fn test_parse_latitude_captures_precision() {
        let c = parse_latitude("60.4").unwrap();
        assert_eq!(c.value, 60.4);
        assert_eq!(c.user_decimal_places, 1);

        let c = parse_latitude("60.382517").unwrap();
        assert_eq!(c.user_decimal_places, 6);

        let c = parse_latitude("-89").unwrap();
        assert_eq!(c.user_decimal_places, 0);

        // Out of range still rejected.
        assert!(parse_latitude("91").is_err());
        assert!(parse_latitude("-91").is_err());
        assert!(parse_latitude("not a number").is_err());
    }
}
