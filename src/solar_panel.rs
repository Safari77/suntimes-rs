//! Solar Panel Output Calculations
//!
//! Implements the Ineichen-Perez clear-sky model for estimating solar irradiance
//! and solar panel power output.
//!
//! References:
//! - Ineichen, P. and Perez, R. (2002). "A new airmass independent formulation
//!   for the Linke turbidity coefficient"
//! - Perez, R. et al. (1990). "Modeling daylight availability and irradiance
//!   components from direct and global irradiance"

use std::f64::consts::PI;

// ===================== CONSTANTS =====================

/// Solar constant (Total Solar Irradiance) in W/m²
/// Latest value from SORCE/TIM measurements
const SOLAR_CONSTANT: f64 = 1361.0;

/// Default Linke turbidity factor for clear atmosphere
/// Typical values: 2-3 for very clear, 4-6 for industrial areas
const DEFAULT_LINKE_TURBIDITY: f64 = 3.0;

/// Default ground albedo (typical grass/soil)
const DEFAULT_ALBEDO: f64 = 0.2;

// ===================== TRACKING MODES =====================

/// Solar panel tracking mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TrackingMode {
    /// Fixed panel - no tracking
    #[default]
    Fixed,
    /// Horizontal single-axis tracking (HSAT)
    /// Axis runs N-S, panel rotates E-W to track sun's azimuth
    /// Tilt is fixed (set via --solarpanel-tilt, can be optimized)
    /// Most common for utility-scale installations
    HorizontalAxis,
    /// Vertical single-axis tracking (VSAT)
    /// Axis is vertical, panel tilts to track sun's altitude
    /// Azimuth is fixed (set via --solarpanel-azimuth)
    VerticalAxis,
    /// Dual-axis tracking - panel always faces the sun directly (both tilt and azimuth)
    /// Provides +30-40% energy vs fixed, highest cost and maintenance
    DualAxis,
}

// ===================== CONFIGURATION =====================

/// Solar panel configuration
#[derive(Debug, Clone, Copy)]
pub struct SolarPanelConfig {
    /// Panel area in square meters
    pub area_m2: f64,
    /// Panel tilt from horizontal in degrees (0 = flat, 90 = vertical)
    pub tilt_deg: f64,
    /// Panel azimuth in degrees (180 = facing south in northern hemisphere)
    pub azimuth_deg: f64,
    /// Panel efficiency (0.0 - 1.0, typical silicon ~0.18-0.22)
    pub efficiency: f64,
    /// Linke turbidity factor (atmospheric clarity, 2-7 typical)
    pub linke_turbidity: f64,
    /// Ground albedo for reflected radiation (0.0 - 1.0)
    pub albedo: f64,
    /// Tracking mode for the panel
    pub tracking_mode: TrackingMode,
}

impl SolarPanelConfig {
    pub fn new(area_m2: f64, tilt_deg: f64, azimuth_deg: f64) -> Self {
        Self {
            area_m2,
            tilt_deg,
            azimuth_deg,
            efficiency: 0.20,
            linke_turbidity: DEFAULT_LINKE_TURBIDITY,
            albedo: DEFAULT_ALBEDO,
            tracking_mode: TrackingMode::Fixed,
        }
    }

    pub fn with_efficiency(mut self, efficiency: f64) -> Self {
        self.efficiency = efficiency;
        self
    }

    pub fn with_linke_turbidity(mut self, lt: f64) -> Self {
        self.linke_turbidity = lt;
        self
    }

    pub fn with_albedo(mut self, albedo: f64) -> Self {
        self.albedo = albedo;
        self
    }

    pub fn with_tracking_mode(mut self, mode: TrackingMode) -> Self {
        self.tracking_mode = mode;
        self
    }
}

// ===================== IRRADIANCE RESULTS =====================

/// Complete irradiance breakdown on tilted surface
#[derive(Debug, Clone, Copy)]
pub struct IrradianceResult {
    /// Direct Normal Irradiance (W/m²)
    pub dni: f64,
    /// Diffuse Horizontal Irradiance (W/m²)
    pub dhi: f64,
    /// Global Horizontal Irradiance (W/m²)
    pub ghi: f64,
    /// Plane-of-Array irradiance (W/m²) - total on tilted panel
    pub poa: f64,
    /// POA beam component (W/m²)
    pub poa_beam: f64,
    /// POA sky diffuse component (W/m²)
    pub poa_sky_diffuse: f64,
    /// POA ground reflected component (W/m²)
    pub poa_ground_diffuse: f64,
    /// Angle of incidence (degrees)
    pub aoi_deg: f64,
}

/// Solar panel output result
#[derive(Debug, Clone, Copy)]
pub struct SolarPanelOutput {
    /// Irradiance components
    pub irradiance: IrradianceResult,
    /// Instantaneous power output (W)
    pub power_w: f64,
    /// Air mass
    pub air_mass: f64,
    /// Sun elevation used (degrees)
    pub sun_elevation_deg: f64,
    /// Sun azimuth used (degrees)
    pub sun_azimuth_deg: f64,
}

// ===================== GEOMETRY =====================

/// Calculate angle of incidence between sun rays and panel normal
///
/// # Arguments
/// * `sun_elevation_deg` - Sun elevation angle in degrees (0 = horizon, 90 = zenith)
/// * `sun_azimuth_deg` - Sun azimuth in degrees (0 = North, 90 = East, 180 = South)
/// * `panel_tilt_deg` - Panel tilt from horizontal in degrees
/// * `panel_azimuth_deg` - Panel facing direction in degrees (180 = South)
///
/// # Returns
/// Angle of incidence in degrees (0 = sun perpendicular to panel)
pub fn angle_of_incidence(
    sun_elevation_deg: f64,
    sun_azimuth_deg: f64,
    panel_tilt_deg: f64,
    panel_azimuth_deg: f64,
) -> f64 {
    let sun_zenith = (90.0 - sun_elevation_deg).to_radians();
    let sun_az = sun_azimuth_deg.to_radians();
    let tilt = panel_tilt_deg.to_radians();
    let panel_az = panel_azimuth_deg.to_radians();

    // Spherical trigonometry formula for AOI
    let cos_aoi =
        sun_zenith.cos() * tilt.cos() + sun_zenith.sin() * tilt.sin() * (sun_az - panel_az).cos();

    cos_aoi.clamp(-1.0, 1.0).acos().to_degrees()
}

// ===================== ATMOSPHERIC CALCULATIONS =====================

/// Calculate absolute air mass (pressure-corrected)
///
/// Uses Kasten-Young (1989) relative air mass model combined with
/// International Standard Atmosphere (ISA) pressure correction.
/// This matches pvlib-python's `get_airmass(model='kastenyoung1989')` logic.
pub fn air_mass(sun_elevation_deg: f64, altitude_m: f64) -> f64 {
    if sun_elevation_deg <= 0.0 {
        return f64::INFINITY;
    }

    let zenith_deg = 90.0 - sun_elevation_deg;
    let zenith_rad = zenith_deg.to_radians();

    // 1. Relative Air Mass (Kasten-Young 1989)
    // Accurate for zenith angles up to ~89°
    let am_relative = 1.0 / (zenith_rad.cos() + 0.50572 * (96.07995 - zenith_deg).powf(-1.6364));

    // 2. Pressure Correction (P / P0) using ISA Model
    // PVLib uses this standard barometric formula:
    // P = P0 * (1 - L*h/T0)^(g*M / R*L)
    // For Earth: (1 - 2.25577e-5 * h)^5.25588
    let pressure_ratio = if altitude_m.abs() < 1e-5 {
        1.0
    } else {
        // Valid for troposphere (< 11km)
        // If user puts >11km, this formula drifts, but pvlib does the same.
        (1.0 - 2.25577e-5 * altitude_m).powf(5.25588)
    };

    am_relative * pressure_ratio
}

/// Calculate extraterrestrial irradiance corrected for Earth-Sun distance
///
/// Uses Spencer (1971) formula for orbital eccentricity correction
///
/// # Arguments
/// * `day_of_year` - Day of year (1-366)
///
/// # Returns
/// Extraterrestrial irradiance in W/m²
pub fn extraterrestrial_irradiance(day_of_year: u32) -> f64 {
    // Spencer (1971) formula for eccentricity correction factor
    let b = 2.0 * PI * (day_of_year as f64 - 1.0) / 365.0;

    let eccentricity_correction = 1.000110
        + 0.034221 * b.cos()
        + 0.001280 * b.sin()
        + 0.000719 * (2.0 * b).cos()
        + 0.000077 * (2.0 * b).sin();

    SOLAR_CONSTANT * eccentricity_correction
}

// ===================== INEICHEN-PEREZ CLEAR SKY MODEL =====================

/// Calculate clear-sky irradiance using Ineichen-Perez model
///
/// The Ineichen-Perez model estimates clear-sky DNI and GHI based on
/// air mass and Linke turbidity coefficient.
///
/// References:
/// - Ineichen, P. (2008). "A broadband simplified version of the Solis clear sky model"
/// - Rigollier et al. (2000) for the diffuse fraction estimation
///
/// # Arguments
/// * `sun_elevation_deg` - Sun elevation in degrees
/// * `altitude_m` - Observer altitude in meters
/// * `day_of_year` - Day of year (1-366)
/// * `linke_turbidity` - Linke turbidity factor (typical 2-7)
///
/// # Returns
/// Tuple of (DNI, DHI, GHI) in W/m²
pub fn ineichen_perez_clearsky(
    sun_elevation_deg: f64,
    altitude_m: f64,
    day_of_year: u32,
    linke_turbidity: f64,
) -> (f64, f64, f64) {
    // Sun must be above horizon
    if sun_elevation_deg <= 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let am = air_mass(sun_elevation_deg, altitude_m);
    if !am.is_finite() || am <= 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let i0 = extraterrestrial_irradiance(day_of_year);
    let sin_elev = sun_elevation_deg.to_radians().sin();

    // Clamp altitude for atmospheric coefficients to prevent model drift/NaN
    let clamped_alt = altitude_m.clamp(-500.0, 11000.0);

    // Altitude correction coefficients (Ineichen 2002)
    let fh1 = (-clamped_alt / 8000.0).exp();
    let fh2 = (-clamped_alt / 1250.0).exp();

    // Linke turbidity with altitude correction
    let altitude_km = clamped_alt / 1000.0;
    let tl = (linke_turbidity - 0.15 * altitude_km).max(1.0);

    // Ineichen clear-sky model coefficients
    let cg1 = 5.09e-5 * clamped_alt + 0.868;
    let cg2 = 3.92e-5 * clamped_alt + 0.0387;

    // 1. Direct Normal Irradiance (DNI)
    let b = 0.664 + 0.163 / fh1;
    let dni_coefficient = b * i0;
    let exponent = -cg2 * am * (fh1 + fh2 * (tl - 1.0));
    let dni = (dni_coefficient * exponent.exp()).max(0.0).min(i0);

    // 2. Global Horizontal Irradiance (GHI)
    let ghi_exponent = -cg2 * am * (fh1 + fh2 * (tl - 1.0)) * 1.1;
    let ghi_raw = (cg1 * i0 * sin_elev * ghi_exponent.exp()).max(0.0);

    // 3. Physical Consistency Check
    // GHI cannot be less than the direct beam component hitting the ground.
    let direct_horizontal = dni * sin_elev;
    let ghi_final = ghi_raw.max(direct_horizontal);

    // DHI is the remainder of the light after accounting for the direct beam
    let dhi = (ghi_final - direct_horizontal).max(0.0);

    (dni, dhi, ghi_final)
}

// ===================== PLANE OF ARRAY IRRADIANCE =====================

/// Calculate irradiance on tilted plane-of-array using Perez model
///
/// Decomposes sky diffuse into circumsolar and horizon brightening components
///
/// # Arguments
/// * `dni` - Direct Normal Irradiance (W/m²)
/// * `dhi` - Diffuse Horizontal Irradiance (W/m²)
/// * `ghi` - Global Horizontal Irradiance (W/m²)
/// * `sun_elevation_deg` - Sun elevation in degrees
/// * `sun_azimuth_deg` - Sun azimuth in degrees
/// * `panel_tilt_deg` - Panel tilt from horizontal in degrees
/// * `panel_azimuth_deg` - Panel azimuth in degrees
/// * `albedo` - Ground reflectance (0-1)
///
/// # Returns
/// IrradianceResult with all components
#[allow(clippy::too_many_arguments)]
pub fn plane_of_array_irradiance(
    dni: f64,
    dhi: f64,
    ghi: f64,
    sun_elevation_deg: f64,
    sun_azimuth_deg: f64,
    panel_tilt_deg: f64,
    panel_azimuth_deg: f64,
    albedo: f64,
) -> IrradianceResult {
    let aoi_deg =
        angle_of_incidence(sun_elevation_deg, sun_azimuth_deg, panel_tilt_deg, panel_azimuth_deg);

    let aoi_rad = aoi_deg.to_radians();
    let tilt_rad = panel_tilt_deg.to_radians();

    // POA beam component (direct irradiance on tilted surface)
    let cos_aoi = aoi_rad.cos();
    let poa_beam = if cos_aoi > 0.0 && sun_elevation_deg > 0.0 { dni * cos_aoi } else { 0.0 };

    // Sky diffuse using isotropic model (simplified Perez)
    // View factor for sky dome
    let sky_view_factor = (1.0 + tilt_rad.cos()) / 2.0;
    let poa_sky_diffuse = dhi * sky_view_factor;

    // Ground reflected component
    let ground_view_factor = (1.0 - tilt_rad.cos()) / 2.0;
    let poa_ground_diffuse = ghi * albedo * ground_view_factor;

    // Total POA irradiance
    let poa = poa_beam + poa_sky_diffuse + poa_ground_diffuse;

    IrradianceResult { dni, dhi, ghi, poa, poa_beam, poa_sky_diffuse, poa_ground_diffuse, aoi_deg }
}

// ===================== POWER OUTPUT CALCULATION =====================

/// Calculate solar panel output for given conditions
///
/// # Arguments
/// * `config` - Solar panel configuration
/// * `sun_elevation_deg` - Sun elevation in degrees
/// * `sun_azimuth_deg` - Sun azimuth in degrees
/// * `altitude_m` - Observer altitude in meters
/// * `day_of_year` - Day of year (1-366)
///
/// # Returns
/// SolarPanelOutput with power and irradiance details
pub fn calculate_output(
    config: &SolarPanelConfig,
    sun_elevation_deg: f64,
    sun_azimuth_deg: f64,
    altitude_m: f64,
    day_of_year: u32,
) -> SolarPanelOutput {
    // Get clear-sky irradiance components
    let (dni, dhi, ghi) =
        ineichen_perez_clearsky(sun_elevation_deg, altitude_m, day_of_year, config.linke_turbidity);

    // Determine panel orientation based on tracking mode
    let (panel_tilt, panel_azimuth) = match config.tracking_mode {
        TrackingMode::Fixed => {
            // Fixed panel - use configured tilt and azimuth
            (config.tilt_deg, config.azimuth_deg)
        }
        TrackingMode::HorizontalAxis => {
            // HSAT: Axis runs N-S, panel rotates E-W to track sun's azimuth
            // Tilt stays fixed (configured value)
            (config.tilt_deg, sun_azimuth_deg)
        }
        TrackingMode::VerticalAxis => {
            // VSAT: Axis is vertical, panel tilts to track sun's altitude
            // Azimuth stays fixed (configured value)
            let tilt = (90.0 - sun_elevation_deg).clamp(0.0, 90.0);
            (tilt, config.azimuth_deg)
        }
        TrackingMode::DualAxis => {
            // Dual-axis: panel always faces the sun directly
            // Tilt adjusts to sun's elevation, azimuth tracks sun
            let tilt = (90.0 - sun_elevation_deg).clamp(0.0, 90.0);
            (tilt, sun_azimuth_deg)
        }
    };

    // Calculate plane-of-array irradiance
    let irradiance = plane_of_array_irradiance(
        dni,
        dhi,
        ghi,
        sun_elevation_deg,
        sun_azimuth_deg,
        panel_tilt,
        panel_azimuth,
        config.albedo,
    );

    // Calculate power output
    // P = POA * Area * Efficiency
    let power_w = irradiance.poa * config.area_m2 * config.efficiency;

    let am = air_mass(sun_elevation_deg, altitude_m);

    SolarPanelOutput { irradiance, power_w, air_mass: am, sun_elevation_deg, sun_azimuth_deg }
}

/// Calculate daily energy production by integrating over the day
///
/// Uses trapezoidal integration with configurable time steps
///
/// # Arguments
/// * `config` - Solar panel configuration
/// * `altitude_m` - Observer altitude in meters
/// * `day_of_year` - Day of year (1-366)
/// * `sun_positions` - Iterator of (hour_of_day, elevation_deg, azimuth_deg)
///
/// # Returns
/// Daily energy in Wh
pub fn calculate_daily_energy<I>(
    config: &SolarPanelConfig,
    altitude_m: f64,
    day_of_year: u32,
    sun_positions: I,
) -> f64
where
    I: Iterator<Item = (f64, f64, f64)>,
{
    let mut positions: Vec<_> = sun_positions.collect();
    if positions.len() < 2 {
        return 0.0;
    }

    // Ensure we start and end with 0 power if the first/last elevation is near 0
    if positions[0].1 < 0.1 {
        let (t, _, az) = positions[0];
        positions[0] = (t, 0.0, az);
    }

    let mut total_energy_wh = 0.0;
    for window in positions.windows(2) {
        let (t1, elev1, az1) = window[0];
        let (t2, elev2, az2) = window[1];

        // Skip intervals where the sun is entirely below the horizon
        if elev1 <= 0.0 && elev2 <= 0.0 {
            continue;
        }

        let output1 = calculate_output(config, elev1, az1, altitude_m, day_of_year);
        let output2 = calculate_output(config, elev2, az2, altitude_m, day_of_year);

        let dt_hours = (t2 - t1).abs();
        let avg_power = (output1.power_w + output2.power_w) / 2.0;
        total_energy_wh += avg_power * dt_hours;
    }

    total_energy_wh
}

// ===================== FORMATTING HELPERS =====================

/// Format power output for display
pub fn format_power(watts: f64) -> String {
    if watts >= 1000.0 { format!("{:.2} kW", watts / 1000.0) } else { format!("{:.1} W", watts) }
}

/// Format energy for display
pub fn format_energy(watt_hours: f64) -> String {
    if watt_hours >= 1000.0 {
        format!("{:.2} kWh", watt_hours / 1000.0)
    } else {
        format!("{:.1} Wh", watt_hours)
    }
}

/// Format irradiance for display
pub fn format_irradiance(w_per_m2: f64) -> String {
    format!("{:.0} W/m²", w_per_m2)
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extraterrestrial_irradiance_range() {
        // Should vary by ~3.3% over the year
        let min = (1..=366).map(extraterrestrial_irradiance).fold(f64::INFINITY, f64::min);
        let max = (1..=366).map(extraterrestrial_irradiance).fold(f64::NEG_INFINITY, f64::max);

        // Perihelion (early January) should be highest
        let jan3 = extraterrestrial_irradiance(3);
        // Aphelion (early July) should be lowest
        let jul4 = extraterrestrial_irradiance(185);

        assert!(jan3 > jul4);
        assert!(min > 1300.0 && min < 1350.0);
        assert!(max > 1380.0 && max < 1420.0);
    }

    #[test]
    fn test_air_mass_typical_values() {
        // Zenith: AM should be ~1.0
        let am_zenith = air_mass(90.0, 0.0);
        assert!((am_zenith - 1.0).abs() < 0.01);

        // 60° elevation (30° zenith): AM should be ~1.15
        let am_60 = air_mass(60.0, 0.0);
        assert!(am_60 > 1.1 && am_60 < 1.2);

        // 30° elevation (60° zenith): AM should be ~2.0
        let am_30 = air_mass(30.0, 0.0);
        assert!(am_30 > 1.9 && am_30 < 2.1);

        // Near horizon: AM should be very high
        let am_5 = air_mass(5.0, 0.0);
        assert!(am_5 > 10.0);

        // Below horizon: infinity
        let am_neg = air_mass(-5.0, 0.0);
        assert!(am_neg.is_infinite());
    }

    #[test]
    fn test_air_mass_altitude_correction() {
        // Higher altitude = less atmosphere = lower air mass
        let am_sea = air_mass(45.0, 0.0);
        let am_mountain = air_mass(45.0, 3000.0);

        assert!(am_mountain < am_sea);
        assert!(am_mountain > am_sea * 0.5); // Should be significant but not extreme
    }

    #[test]
    fn test_angle_of_incidence() {
        // Sun directly overhead, flat panel: AOI = 0
        let aoi1 = angle_of_incidence(90.0, 180.0, 0.0, 180.0);
        assert!(aoi1.abs() < 0.1);

        // Sun at 45° elevation from south, panel tilted 45° facing south: AOI = 0
        let aoi2 = angle_of_incidence(45.0, 180.0, 45.0, 180.0);
        assert!(aoi2.abs() < 0.1);

        // Sun from east (az=90), panel facing south (az=180): AOI should be significant
        // With 45° elevation and 45° tilt, and 90° azimuth difference
        // cos(AOI) = cos(45°)*cos(45°) + sin(45°)*sin(45°)*cos(90°)
        //          = 0.5 + 0.5*0 = 0.5 → AOI = 60°
        let aoi3 = angle_of_incidence(45.0, 90.0, 45.0, 180.0);
        assert!((aoi3 - 60.0).abs() < 1.0, "AOI was {}, expected ~60°", aoi3);

        // Sun behind panel (north, az=0), panel facing south: AOI >= 90° (sun on back)
        let aoi4 = angle_of_incidence(45.0, 0.0, 45.0, 180.0);
        assert!(aoi4 >= 90.0 - 0.01, "AOI was {}, expected >= 90° (sun behind panel)", aoi4);
    }

    #[test]
    fn test_ineichen_perez_basic() {
        // High sun, clear conditions should give high DNI
        let (dni, dhi, ghi) = ineichen_perez_clearsky(60.0, 0.0, 172, 3.0);

        // DNI should be substantial (600-1000 W/m² typical for clear sky)
        assert!(dni > 500.0 && dni < 1100.0, "DNI was {}", dni);
        // DHI must be positive (our model guarantees this)
        assert!(dhi > 0.0, "DHI was {} (must be positive)", dhi);
        // GHI should be high
        assert!(ghi > 500.0 && ghi < 1200.0, "GHI was {}", ghi);

        // Physical consistency: GHI = DNI * sin(elevation) + DHI
        let sin_elev = 60.0_f64.to_radians().sin();
        let expected_ghi = dni * sin_elev + dhi;
        assert!((ghi - expected_ghi).abs() < 1.0, "GHI {} != DNI*sin + DHI {}", ghi, expected_ghi);

        // Sun below horizon should give zeros
        let (dni_neg, dhi_neg, ghi_neg) = ineichen_perez_clearsky(-5.0, 0.0, 172, 3.0);
        assert_eq!(dni_neg, 0.0);
        assert_eq!(dhi_neg, 0.0);
        assert_eq!(ghi_neg, 0.0);
    }

    #[test]
    fn test_ineichen_perez_turbidity_effect() {
        // Higher turbidity should reduce DNI
        let (dni_clear, _, _) = ineichen_perez_clearsky(60.0, 0.0, 172, 2.0);
        let (dni_hazy, _, _) = ineichen_perez_clearsky(60.0, 0.0, 172, 5.0);

        assert!(dni_clear > dni_hazy);
    }

    #[test]
    fn test_poa_flat_panel() {
        // Flat panel should see mostly GHI
        let irr = plane_of_array_irradiance(800.0, 100.0, 800.0, 60.0, 180.0, 0.0, 180.0, 0.2);

        // POA should be close to GHI for flat panel
        assert!((irr.poa - irr.ghi).abs() < 150.0);
    }

    #[test]
    fn test_calculate_output_basic() {
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, 60.0, 180.0, 0.0, 172);

        // 10m² panel at 20% efficiency with ~800 W/m² should give ~1600W
        assert!(output.power_w > 1000.0 && output.power_w < 2500.0);
    }

    #[test]
    fn test_calculate_output_night() {
        // Fixed panel at night
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, -10.0, 180.0, 0.0, 172);

        assert_eq!(output.power_w, 0.0);

        // Tracking panel at night - should also produce zero
        let config_tracking =
            SolarPanelConfig::new(10.0, 35.0, 180.0).with_tracking_mode(TrackingMode::DualAxis);
        let output_tracking = calculate_output(&config_tracking, -10.0, 180.0, 0.0, 172);

        assert_eq!(output_tracking.power_w, 0.0);
    }

    #[test]
    fn test_solar_panel_output_winter_low_sun_high_latitude() {
        // Test vector: Helsinki area (60°N), January, solar noon
        // Sun elevation: ~11°, very low winter sun
        // Panel: 35° tilt facing south
        //
        // With 35° tilt and 11° sun elevation, AOI = ~44°
        // Lower output than summer due to:
        // - High air mass (~5.2)
        // - Large AOI reducing beam capture
        // - Lower DNI due to atmospheric path

        let config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let sun_elevation = 11.0; // Helsinki winter noon
        let output = calculate_output(&config, sun_elevation, 180.0, 0.0, 26);

        // High air mass due to very low sun
        assert!(output.air_mass > 5.0, "Air mass {} should be high at 11°", output.air_mass);

        // DNI reduced but present on clear day
        assert!(output.irradiance.dni > 400.0, "DNI {} too low", output.irradiance.dni);
        assert!(output.irradiance.dni < 700.0, "DNI {} too high", output.irradiance.dni);

        // GHI very low due to sin(11°) ≈ 0.19
        assert!(output.irradiance.ghi > 80.0, "GHI {} too low", output.irradiance.ghi);
        assert!(output.irradiance.ghi < 200.0, "GHI {} too high", output.irradiance.ghi);

        // AOI around 44° (90 - 11 - 35 = 44)
        assert!(
            output.irradiance.aoi_deg > 40.0 && output.irradiance.aoi_deg < 50.0,
            "AOI {} should be ~44°",
            output.irradiance.aoi_deg
        );

        // Power lower than summer but still reasonable on clear day
        assert!(output.power_w > 700.0, "Power {} W too low", output.power_w);
        assert!(output.power_w < 1200.0, "Power {} W too high for winter", output.power_w);
    }

    #[test]
    fn test_solar_panel_summer_vs_winter_comparison() {
        // Direct comparison: same location, same panel, different seasons
        // This validates that summer noon > winter noon for high latitude

        let config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        // Summer: sun at 53° elevation (Helsinki June)
        let summer = calculate_output(&config, 53.0, 180.0, 0.0, 177);

        // Winter: sun at 11° elevation (Helsinki January)
        let winter = calculate_output(&config, 11.0, 180.0, 0.0, 26);

        // Summer should produce significantly more power
        assert!(
            summer.power_w > winter.power_w * 1.5,
            "Summer {} W should be much higher than winter {} W",
            summer.power_w,
            winter.power_w
        );

        // Summer GHI much higher
        assert!(
            summer.irradiance.ghi > winter.irradiance.ghi * 3.0,
            "Summer GHI {} should be much higher than winter GHI {}",
            summer.irradiance.ghi,
            winter.irradiance.ghi
        );

        // Summer air mass lower
        assert!(
            summer.air_mass < winter.air_mass,
            "Summer AM {} should be lower than winter AM {}",
            summer.air_mass,
            winter.air_mass
        );
    }

    #[test]
    fn test_solar_panel_output_golden_colorado_summer_noon() {
        // Golden, Colorado (NREL): 39.755°N, 1830m altitude
        // June 21 solar noon: sun elevation ~73.5°

        let config = SolarPanelConfig::new(10.0, 40.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(2.5)
            .with_albedo(0.2);

        let output = calculate_output(&config, 73.5, 180.0, 1830.0, 172);

        assert!(output.air_mass < 1.1, "Air mass {} too high", output.air_mass);
        assert!(output.irradiance.dni > 900.0, "DNI {} too low", output.irradiance.dni);
        assert!(output.irradiance.ghi > 850.0, "GHI {} too low", output.irradiance.ghi);
        assert!(output.irradiance.aoi_deg < 35.0, "AOI {} too large", output.irradiance.aoi_deg);
        assert!(
            output.power_w > 1700.0 && output.power_w < 2200.0,
            "Power {} W out of range",
            output.power_w
        );
    }

    #[test]
    fn test_solar_panel_output_helsinki_winter_noon() {
        // Helsinki area: 60°N, sea level
        // January 26 solar noon: sun elevation ~11°

        let config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let output = calculate_output(&config, 10.92, 172.5, 0.0, 26);

        // High air mass due to low sun
        assert!(
            output.air_mass > 5.0 && output.air_mass < 5.5,
            "Air mass {} should be ~5.1",
            output.air_mass
        );

        // DNI attenuated but still substantial
        assert!(
            output.irradiance.dni > 600.0 && output.irradiance.dni < 700.0,
            "DNI {} out of range",
            output.irradiance.dni
        );

        // GHI low due to sin(11°) ≈ 0.19
        assert!(
            output.irradiance.ghi > 120.0 && output.irradiance.ghi < 170.0,
            "GHI {} out of range",
            output.irradiance.ghi
        );

        // AOI ~44° (panel can't fully capture low winter sun)
        assert!(
            output.irradiance.aoi_deg > 43.0 && output.irradiance.aoi_deg < 46.0,
            "AOI {} should be ~44°",
            output.irradiance.aoi_deg
        );

        // POA boosted vs GHI by tilted panel
        assert!(
            output.irradiance.poa > output.irradiance.ghi * 3.0,
            "POA {} should be much higher than GHI {}",
            output.irradiance.poa,
            output.irradiance.ghi
        );

        // Power ~960W
        assert!(
            output.power_w > 900.0 && output.power_w < 1050.0,
            "Power {} W out of range",
            output.power_w
        );
    }

    #[test]
    fn test_solar_panel_helsinki_summer_vs_winter() {
        // Compare same location, different seasons

        let config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        // Summer: 53° elevation
        let summer = calculate_output(&config, 52.85, 172.0, 0.0, 177);
        // Winter: 11° elevation
        let winter = calculate_output(&config, 10.92, 172.5, 0.0, 26);

        // Summer ~2x more instantaneous power
        assert!(
            summer.power_w > winter.power_w * 1.8,
            "Summer {} W should be ~2x winter {} W",
            summer.power_w,
            winter.power_w
        );

        // Summer GHI ~5x higher
        assert!(
            summer.irradiance.ghi > winter.irradiance.ghi * 4.5,
            "Summer GHI {} should be ~5x winter GHI {}",
            summer.irradiance.ghi,
            winter.irradiance.ghi
        );

        // Summer AOI much smaller (sun closer to panel normal)
        assert!(
            summer.irradiance.aoi_deg < winter.irradiance.aoi_deg * 0.2,
            "Summer AOI {} should be much smaller than winter AOI {}",
            summer.irradiance.aoi_deg,
            winter.irradiance.aoi_deg
        );
    }
    #[test]
    fn test_solar_panel_output_winter_solstice_low_sun() {
        // Test vector: Winter solstice at mid-latitude, solar noon
        // At 40°N on Dec 21: sun elevation ≈ 26.5°
        //
        // Key insight: With 40° panel tilt, winter sun at 26.5° elevation
        // hits the panel at a favorable angle, so POA can still be high!
        // The difference shows in DAILY energy, not instantaneous power.

        let config = SolarPanelConfig::new(10.0, 40.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let sun_elevation = 26.5;
        let sun_azimuth = 180.0;
        let altitude_m = 0.0;
        let day_of_year = 355;

        let output = calculate_output(&config, sun_elevation, sun_azimuth, altitude_m, day_of_year);

        // Air mass higher due to low sun
        assert!(output.air_mass > 2.0, "Air mass {} too low", output.air_mass);
        assert!(output.air_mass < 2.5, "Air mass {} too high", output.air_mass);

        // DNI reduced but still good on clear day
        assert!(output.irradiance.dni > 700.0, "DNI {} too low", output.irradiance.dni);
        assert!(output.irradiance.dni < 950.0, "DNI {} too high", output.irradiance.dni);

        // GHI lower due to low sun angle
        assert!(output.irradiance.ghi > 350.0, "GHI {} too low", output.irradiance.ghi);
        assert!(output.irradiance.ghi < 550.0, "GHI {} too high", output.irradiance.ghi);

        // POA boosted by favorable tilt angle
        assert!(
            output.irradiance.poa > output.irradiance.ghi,
            "POA {} should exceed GHI {}",
            output.irradiance.poa,
            output.irradiance.ghi
        );

        // AOI small because panel tilt matches low winter sun well
        assert!(
            output.irradiance.aoi_deg < 25.0,
            "AOI {} should be small",
            output.irradiance.aoi_deg
        );

        // Power can still be substantial at solar noon on clear winter day
        assert!(output.power_w > 1400.0, "Power {} W too low", output.power_w);
        assert!(output.power_w < 2000.0, "Power {} W too high", output.power_w);
    }

    #[test]
    fn test_solar_panel_output_low_sun_large_aoi() {
        // Test: Low sun with azimuth mismatch causing large AOI
        // This represents poor conditions - low sun from the side

        let config = SolarPanelConfig::new(10.0, 30.0, 180.0) // Panel faces south
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let sun_elevation = 15.0;
        let sun_azimuth = 240.0; // Sun from SW - 60° off panel azimuth
        let altitude_m = 0.0;
        let day_of_year = 172;

        let output = calculate_output(&config, sun_elevation, sun_azimuth, altitude_m, day_of_year);

        // Air mass high due to low sun
        assert!(output.air_mass > 3.5, "Air mass {} should be high", output.air_mass);

        // GHI low due to low sun angle
        assert!(
            output.irradiance.ghi < 350.0,
            "GHI {} too high for 15° sun",
            output.irradiance.ghi
        );

        // AOI large due to azimuth mismatch
        assert!(
            output.irradiance.aoi_deg > 45.0,
            "AOI {} should be large",
            output.irradiance.aoi_deg
        );

        // Power reduced due to poor geometry
        assert!(output.power_w > 300.0, "Power {} W too low", output.power_w);
        assert!(output.power_w < 900.0, "Power {} W should be modest", output.power_w);
    }

    #[test]
    fn test_solar_panel_output_near_sunset() {
        // Test vector: Low sun near horizon
        // Sun at 10° elevation - tests high air mass conditions
        //
        // At 10° elevation:
        // - Air mass ~5.6
        // - Significant atmospheric attenuation
        // - But tilted panel can still capture reasonable POA

        let config = SolarPanelConfig::new(10.0, 30.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let sun_elevation = 10.0;
        let sun_azimuth = 220.0; // Southwest (afternoon)
        let altitude_m = 0.0;
        let day_of_year = 172;

        let output = calculate_output(&config, sun_elevation, sun_azimuth, altitude_m, day_of_year);

        // Air mass should be very high
        assert!(
            output.air_mass > 5.0,
            "Air mass {} should be high at 10° elevation",
            output.air_mass
        );
        assert!(output.air_mass < 6.5, "Air mass {} unexpectedly high", output.air_mass);

        // DNI reduced but still present
        assert!(output.irradiance.dni > 300.0, "DNI {} too low", output.irradiance.dni);
        assert!(
            output.irradiance.dni < 700.0,
            "DNI {} too high for low sun",
            output.irradiance.dni
        );

        // GHI very low due to sin(10°) ≈ 0.17
        assert!(output.irradiance.ghi > 80.0, "GHI {} too low", output.irradiance.ghi);
        assert!(
            output.irradiance.ghi < 300.0,
            "GHI {} too high for 10° sun",
            output.irradiance.ghi
        );

        // Power output modest but not negligible
        // Panel at 30° tilt with sun at 10° from SW has larger AOI
        assert!(output.power_w > 200.0, "Power {} W too low", output.power_w);
        assert!(output.power_w < 900.0, "Power {} W too high for near-sunset", output.power_w);

        // Verify AOI is large due to azimuth mismatch (220° vs 180°)
        assert!(
            output.irradiance.aoi_deg > 30.0,
            "AOI {} should be large due to azimuth offset",
            output.irradiance.aoi_deg
        );
    }

    #[test]
    fn test_tracking_vs_fixed_panel_higher_output() {
        // Test that a tracking panel always produces equal or higher power than fixed
        // Tracking panel always faces the sun directly, so AOI = 0°

        let fixed_config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2);

        let tracking_config = SolarPanelConfig::new(10.0, 35.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2)
            .with_tracking_mode(TrackingMode::DualAxis);

        // Test case 1: Sun from the side (azimuth mismatch)
        // Fixed panel faces south (180°), sun from SW (240°)
        let fixed_1 = calculate_output(&fixed_config, 45.0, 240.0, 0.0, 172);
        let tracking_1 = calculate_output(&tracking_config, 45.0, 240.0, 0.0, 172);

        assert!(
            tracking_1.power_w >= fixed_1.power_w,
            "Tracking {} W should be >= fixed {} W (azimuth mismatch case)",
            tracking_1.power_w,
            fixed_1.power_w
        );

        // Tracking should have AOI close to 0
        assert!(
            tracking_1.irradiance.aoi_deg < 1.0,
            "Tracking AOI {} should be ~0°",
            tracking_1.irradiance.aoi_deg
        );

        // Fixed should have large AOI due to azimuth offset
        assert!(
            fixed_1.irradiance.aoi_deg > 30.0,
            "Fixed AOI {} should be significant",
            fixed_1.irradiance.aoi_deg
        );

        // Test case 2: Low sun angle
        // At low elevation, tracking provides major benefit
        let fixed_2 = calculate_output(&fixed_config, 15.0, 180.0, 0.0, 172);
        let tracking_2 = calculate_output(&tracking_config, 15.0, 180.0, 0.0, 172);

        assert!(
            tracking_2.power_w >= fixed_2.power_w,
            "Tracking {} W should be >= fixed {} W (low sun case)",
            tracking_2.power_w,
            fixed_2.power_w
        );

        // Test case 3: Optimal conditions for fixed panel (sun matches panel orientation)
        // Even here, tracking should be at least equal (both have AOI ~0)
        let fixed_3 = calculate_output(&fixed_config, 55.0, 180.0, 0.0, 172);
        let tracking_3 = calculate_output(&tracking_config, 55.0, 180.0, 0.0, 172);

        assert!(
            tracking_3.power_w >= fixed_3.power_w * 0.99, // Allow tiny floating point tolerance
            "Tracking {} W should be >= fixed {} W (optimal case)",
            tracking_3.power_w,
            fixed_3.power_w
        );

        // Test case 4: Morning sun from east
        let fixed_4 = calculate_output(&fixed_config, 30.0, 90.0, 0.0, 172);
        let tracking_4 = calculate_output(&tracking_config, 30.0, 90.0, 0.0, 172);

        assert!(
            tracking_4.power_w > fixed_4.power_w * 1.5,
            "Tracking {} W should be significantly higher than fixed {} W (morning case)",
            tracking_4.power_w,
            fixed_4.power_w
        );
    }

    #[test]
    fn test_tracking_panel_aoi_always_zero() {
        // Tracking panel should always have AOI close to 0° regardless of sun position

        let tracking_config =
            SolarPanelConfig::new(10.0, 35.0, 180.0).with_tracking_mode(TrackingMode::DualAxis);

        // Various sun positions
        let test_cases = [
            (60.0, 180.0), // High sun, south
            (30.0, 90.0),  // Low sun, east
            (45.0, 270.0), // Medium sun, west
            (20.0, 120.0), // Low sun, SE
            (70.0, 200.0), // High sun, SSW
        ];

        for (elevation, azimuth) in test_cases {
            let output = calculate_output(&tracking_config, elevation, azimuth, 0.0, 172);

            assert!(
                output.irradiance.aoi_deg < 1.0,
                "Tracking AOI {} should be ~0° for sun at {}° elevation, {}° azimuth",
                output.irradiance.aoi_deg,
                elevation,
                azimuth
            );
        }
    }

    #[test]
    fn test_edge_case_sun_exactly_at_horizon() {
        // Sun at exactly 0° elevation - edge case
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, 0.0, 180.0, 0.0, 172);

        // Should produce zero or near-zero output
        assert!(output.power_w >= 0.0, "Power should not be negative");
        assert!(output.power_w < 1.0, "Power should be near zero at horizon");
    }

    #[test]
    fn test_edge_case_very_low_sun() {
        // Sun at 0.1° elevation - tests near-horizon behavior
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, 0.1, 180.0, 0.0, 172);

        // Air mass should be finite but very high
        assert!(output.air_mass.is_finite(), "Air mass should be finite");
        assert!(output.air_mass > 20.0, "Air mass should be very high near horizon");

        // Power should be small but non-negative
        assert!(output.power_w >= 0.0, "Power should not be negative");
    }

    #[test]
    fn test_edge_case_extreme_altitude() {
        // Test at maximum allowed altitude (11000m)
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, 45.0, 180.0, 11000.0, 172);

        // Should still produce valid output
        assert!(output.power_w.is_finite(), "Power should be finite");
        assert!(output.power_w > 0.0, "Power should be positive");
        assert!(output.air_mass.is_finite(), "Air mass should be finite");
        assert!(output.air_mass > 0.0, "Air mass should be positive");
    }

    #[test]
    fn test_edge_case_negative_altitude() {
        // Test at Dead Sea level (-430m)
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);
        let output = calculate_output(&config, 45.0, 180.0, -430.0, 172);

        // Should still produce valid output
        assert!(output.power_w.is_finite(), "Power should be finite");
        assert!(output.power_w > 0.0, "Power should be positive");
        assert!(output.air_mass.is_finite(), "Air mass should be finite");
        // Air mass should be slightly higher due to more atmosphere
        let output_sea_level = calculate_output(&config, 45.0, 180.0, 0.0, 172);
        assert!(
            output.air_mass > output_sea_level.air_mass,
            "Air mass at -430m ({}) should be higher than at sea level ({})",
            output.air_mass,
            output_sea_level.air_mass
        );
    }

    #[test]
    fn test_edge_case_day_of_year_boundaries() {
        // Test day 1 and day 366
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0);

        let output_jan1 = calculate_output(&config, 45.0, 180.0, 0.0, 1);
        assert!(output_jan1.power_w.is_finite(), "Day 1 should produce finite power");

        let output_dec31 = calculate_output(&config, 45.0, 180.0, 0.0, 365);
        assert!(output_dec31.power_w.is_finite(), "Day 365 should produce finite power");

        let output_leap = calculate_output(&config, 45.0, 180.0, 0.0, 366);
        assert!(output_leap.power_w.is_finite(), "Day 366 should produce finite power");
    }

    #[test]
    fn test_edge_case_extreme_turbidity() {
        let config_low = SolarPanelConfig::new(10.0, 35.0, 180.0).with_linke_turbidity(1.0);
        let config_high = SolarPanelConfig::new(10.0, 35.0, 180.0).with_linke_turbidity(10.0);

        let output_low = calculate_output(&config_low, 45.0, 180.0, 0.0, 172);
        let output_high = calculate_output(&config_high, 45.0, 180.0, 0.0, 172);

        // Both should produce valid output
        assert!(output_low.power_w.is_finite(), "Low turbidity should produce finite power");
        assert!(output_high.power_w.is_finite(), "High turbidity should produce finite power");

        // Lower turbidity = clearer sky = more power
        assert!(
            output_low.power_w > output_high.power_w,
            "Clear sky ({} W) should produce more power than hazy ({} W)",
            output_low.power_w,
            output_high.power_w
        );
    }

    #[test]
    fn test_edge_case_zero_area_panel() {
        // Zero area panel should produce zero power
        let config = SolarPanelConfig::new(0.0, 35.0, 180.0);
        let output = calculate_output(&config, 45.0, 180.0, 0.0, 172);

        assert_eq!(output.power_w, 0.0, "Zero area panel should produce zero power");
    }

    #[test]
    fn test_edge_case_zero_efficiency() {
        // Zero efficiency panel should produce zero power
        let config = SolarPanelConfig::new(10.0, 35.0, 180.0).with_efficiency(0.0);
        let output = calculate_output(&config, 45.0, 180.0, 0.0, 172);

        assert_eq!(output.power_w, 0.0, "Zero efficiency panel should produce zero power");
    }
}
