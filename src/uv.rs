//! UV Index Calculation Module
//!
//! Based on empirical models for Ozone column estimation and UV transmission.

/// Global monthly ozone climatology (Dobson Units) dataset based on modern OMI/TOMS satellite data
/// Rows: Latitude zones from -90 to +90 (increments of 20°)
/// Cols: Months Jan to Dec
const OZONE_GRID: [[f64; 12]; 9] = [
    // -90 to -70 (South Pole)
    [280.0, 270.0, 260.0, 250.0, 250.0, 250.0, 250.0, 220.0, 150.0, 140.0, 180.0, 250.0],
    // -70 to -50
    [300.0, 290.0, 290.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 330.0, 320.0, 310.0],
    // -50 to -30
    [280.0, 275.0, 275.0, 280.0, 290.0, 310.0, 330.0, 340.0, 350.0, 340.0, 320.0, 295.0],
    // -30 to -10 (Subtropical SH)
    [265.0, 260.0, 260.0, 260.0, 265.0, 270.0, 280.0, 285.0, 290.0, 285.0, 275.0, 270.0],
    // -10 to +10 (Equatorial)
    [250.0, 255.0, 260.0, 260.0, 260.0, 255.0, 255.0, 260.0, 260.0, 260.0, 255.0, 250.0],
    // +10 to +30 (Subtropical NH)
    [265.0, 275.0, 285.0, 295.0, 300.0, 295.0, 285.0, 280.0, 275.0, 265.0, 260.0, 260.0],
    // +30 to +50 (Mid-latitude NH)
    [320.0, 340.0, 360.0, 365.0, 340.0, 315.0, 300.0, 290.0, 285.0, 285.0, 295.0, 310.0],
    // +50 to +70
    [360.0, 380.0, 400.0, 395.0, 365.0, 335.0, 315.0, 300.0, 290.0, 290.0, 310.0, 330.0],
    // +70 to +90 (North Pole)
    [380.0, 400.0, 420.0, 410.0, 380.0, 345.0, 320.0, 300.0, 290.0, 290.0, 310.0, 340.0],
];

/// Estimate total ozone column in Dobson Units (DU)
///
/// Uses a monthly gridded climatology for better regional accuracy,
/// with a tropospheric altitude correction on top.
///
/// Note: This function is called from inner loops in yearly optimization, so
/// it deliberately performs no I/O. Callers that want to surface the
/// ozone-hole caveat should check [`is_ozone_hole_season`] once at the top
/// level and print their own message.
pub fn estimate_ozone_column(latitude: f64, month: u32, altitude_m: f64) -> f64 {
    // 1. Map latitude to grid row (-90..90 -> 0..8)
    let row_idx = ((latitude + 90.0) / 20.0).floor() as usize;
    let row_idx = row_idx.min(8);

    // 2. Map month to grid column (1..12 -> 0..11)
    let col_idx = (month.clamp(1, 12) - 1) as usize;

    // 3. Lookup base value
    let mut ozone_du = OZONE_GRID[row_idx][col_idx];

    // 4. Tropospheric Altitude Correction
    // Cities at high altitudes have less of the ozone column above them.
    // Clamped at 11km (Troposphere limit).
    let altitude_km = (altitude_m / 1000.0).clamp(0.0, 11.0);
    let elevation_correction = 1.0 - (0.035 * altitude_km);
    ozone_du *= elevation_correction;

    ozone_du.clamp(150.0, 500.0)
}

/// True if the climatology-driven UV index at this latitude & month is likely
/// to underestimate reality due to seasonal Antarctic ozone depletion.
///
/// Southern Hemisphere, poleward of 60°S, September through November.
pub fn is_ozone_hole_season(latitude: f64, month: u32) -> bool {
    latitude < -60.0 && (9..=11).contains(&month)
}

/// Calculate Angstrom Aerosol Optical Depth at specific wavelength
fn angstrom_aod(aod0: f64, lambda0: f64, lambda1: f64, alpha: f64) -> f64 {
    aod0 * (lambda1 / lambda0).powf(-alpha)
}

/// Calculates clear-sky UV Index using an empirical Solar Zenith Angle parameterization adjusted
/// for ozone, aerosols, altitude, and albedo.
pub fn calculate_uv_index(
    solar_elevation: f64,
    airmass: f64,
    latitude: f64,
    month: u32,
    altitude_m: f64,
) -> f64 {
    if solar_elevation <= 0.0 {
        return 0.0;
    }

    let ozone_du = estimate_ozone_column(latitude, month, altitude_m);

    // Solar zenith angle conversion
    let mu = solar_elevation.to_radians().sin();

    // Standard clear-sky empirical UVI model
    let uv_ozone = 12.5 * mu.powf(2.42) * (300.0 / ozone_du).powf(1.2);
    let aod_uv = angstrom_aod(0.04, 550.0, 310.0, 1.3);
    let aerosol_trans = (-0.5 * aod_uv * airmass).exp();
    let albedo_factor = 1.02;

    // Altitude correction
    let altitude_km = altitude_m / 1000.0;
    let altitude_correction = 1.0 + (0.10 * altitude_km);

    let uv_index = uv_ozone * aerosol_trans * albedo_factor * altitude_correction;
    uv_index.clamp(0.0, 16.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ozone_estimation_seasonal_bounds() {
        // NH Mid-latitudes: April (Month 4) should have higher ozone than October (Month 10)
        let spring = estimate_ozone_column(45.0, 4, 0.0);
        let fall = estimate_ozone_column(45.0, 10, 0.0);
        assert!(spring > fall, "NH Spring ozone ({}) should be > Fall ({})", spring, fall);

        // Tropics should be consistently lower than mid-latitudes
        let tropics = estimate_ozone_column(0.0, 1, 0.0);
        assert!(tropics < 300.0);
    }

    #[test]
    fn test_uvi_benchmarks() {
        // Benchmark 1: Tropical Noon (Zenith sun, high GHI, low airmass, low ozone)
        // Result should be "Extreme" (11+)
        let uvi_tropics = calculate_uv_index(90.0, 1.0, 0.0, 3, 0.0);
        assert!(uvi_tropics > 11.0, "Tropical noon UVI should be extreme, got {}", uvi_tropics);

        // Benchmark 2: Mid-latitude Summer Noon (45N, ~68 elev, June/Month 6)
        // Result should be "Very High" (8-10)
        let uvi_summer = calculate_uv_index(68.0, 1.08, 45.0, 6, 0.0);
        assert!(uvi_summer > 7.0 && uvi_summer < 11.0, "Summer UVI {} out of range", uvi_summer);

        // Benchmark 3: Mid-latitude Winter Noon (45N, ~22 elev, Dec/Month 12)
        // Result should be "Low" (< 2)
        let uvi_winter = calculate_uv_index(22.0, 2.6, 45.0, 12, 0.0);
        assert!(uvi_winter < 2.0, "Winter UVI {} must be low (< 2.0)", uvi_winter);
    }

    #[test]
    fn test_stratospheric_uv_resilience() {
        // Sea level: Month 6, altitude 0m
        let uvi_sea = calculate_uv_index(60.0, 1.15, 45.0, 6, 0.0);
        // Stratosphere (11km): Less ozone above city
        let uvi_11km = calculate_uv_index(60.0, 0.25, 45.0, 6, 11000.0);

        assert!(
            uvi_11km > uvi_sea * 1.5,
            "Stratospheric UVI ({}) should be much higher than sea level ({})",
            uvi_11km,
            uvi_sea
        );
    }

    #[test]
    fn test_polar_day_lapland() {
        // Lapland (~67.4°N) June (Month 6)
        let latitude = 67.4;
        let month = 6;
        let sun_elevation = 46.0;
        let airmass = 1.39;

        let uvi = calculate_uv_index(sun_elevation, airmass, latitude, month, 200.0);

        // Expected: Moderate (3-6)
        assert!(uvi >= 3.0 && uvi <= 6.0, "Lapland summer UVI {} should be moderate", uvi);
    }

    #[test]
    fn test_extreme_uv_high_altitude_quito() {
        // Quito (~0.2°S) at ~2850m elevation in March (Month 3)
        let latitude = -0.2;
        let month = 3;
        let sun_elevation = 90.0;
        let airmass = 0.72;

        let uvi = calculate_uv_index(sun_elevation, airmass, latitude, month, 2850.0);

        // Quito has extreme UV (>14) due to low ozone and altitude correction
        assert!(uvi > 14.0, "Quito equinox UVI {} should be extreme (>14)", uvi);
    }

    #[test]
    fn test_australia_summer_perth() {
        // Perth (~31.9°S) in January (Month 1)
        let latitude = -31.9;
        let month = 1;
        let sun_elevation = 81.5;
        let airmass = 1.01;

        let uvi = calculate_uv_index(sun_elevation, airmass, latitude, month, 30.0);

        // Australian summer UVI is notoriously extreme (11+)
        assert!(uvi >= 11.0, "Australian summer UVI {} should be extreme (>=11)", uvi);
    }

    #[test]
    fn test_zero_uvi_at_night() {
        assert_eq!(calculate_uv_index(-5.0, 0.0, 10.0, 6, 0.0), 0.0);
        assert_eq!(calculate_uv_index(0.0, 10.0, 40.0, 6, 0.0), 0.0);
    }
}
