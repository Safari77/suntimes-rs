//! UV Index Calculation Module
//!
//! Based on empirical models for Ozone column estimation and UV transmission.

/// Global monthly ozone climatology (Dobson Units)
/// Rows: Latitude zones from -90 to +90 (increments of 20°)
/// Cols: Months Jan to Dec
const OZONE_GRID: [[f64; 12]; 9] = [
    // -90 to -70 (South Pole)
    [300.0, 280.0, 260.0, 250.0, 260.0, 280.0, 300.0, 320.0, 280.0, 220.0, 250.0, 280.0],
    // -70 to -50
    [320.0, 310.0, 300.0, 290.0, 300.0, 320.0, 340.0, 360.0, 350.0, 330.0, 330.0, 320.0],
    // -50 to -30
    [300.0, 290.0, 285.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 320.0, 310.0, 305.0],
    // -30 to -10 (Subtropical SH)
    [270.0, 265.0, 260.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 285.0, 280.0, 275.0],
    // -10 to +10 (Equatorial)
    [260.0, 260.0, 265.0, 270.0, 270.0, 265.0, 260.0, 260.0, 260.0, 265.0, 265.0, 260.0],
    // +10 to +30 (Subtropical NH)
    [270.0, 280.0, 290.0, 300.0, 300.0, 290.0, 280.0, 275.0, 270.0, 265.0, 265.0, 270.0],
    // +30 to +50 (Mid-latitude NH)
    [320.0, 340.0, 360.0, 370.0, 360.0, 340.0, 320.0, 310.0, 300.0, 290.0, 300.0, 310.0],
    // +50 to +70
    [360.0, 390.0, 420.0, 430.0, 410.0, 380.0, 350.0, 330.0, 320.0, 310.0, 320.0, 340.0],
    // +70 to +90 (North Pole)
    [380.0, 410.0, 440.0, 450.0, 420.0, 380.0, 350.0, 330.0, 320.0, 330.0, 350.0, 370.0],
];

/// Estimate total ozone column in Dobson Units (DU)
///
/// Uses a monthly gridded climatology for better regional accuracy.
/// Now includes altitude correction and data-mismatch warnings.
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

    // 5. Warning Logic: Specific Polar Shifts
    // Warn users in the Southern Hemisphere during Ozone Hole season (Sept-Nov)
    if latitude < -60.0 && (9..=11).contains(&month) {
        eprintln!(
            "Warning: Significant ozone depletion possible at lat {:.1} during month {}. Climatology may overestimate UVI.",
            latitude, month
        );
    }

    ozone_du.clamp(150.0, 500.0)
}

/// Calculate Angstrom Aerosol Optical Depth at specific wavelength
fn angstrom_aod(aod0: f64, lambda0: f64, lambda1: f64, alpha: f64) -> f64 {
    aod0 * (lambda1 / lambda0).powf(-alpha)
}

/// Updated UV Index calculation now accepts altitude to handle ozone correction internally.
pub fn calculate_uv_index(
    solar_elevation: f64,
    clearsky_ghi: f64,
    airmass: f64,
    latitude: f64,
    month: u32,
    altitude_m: f64,
) -> f64 {
    if solar_elevation <= 0.0 {
        return 0.0;
    }

    // Ozone logic is now fully contained within uv.rs
    let ozone_du = estimate_ozone_column(latitude, month, altitude_m);
    let uv_transmission = calculate_uv_transmission(90.0 - solar_elevation, ozone_du, airmass);

    // Base fraction and non-linear elevation adjustment
    let base_fraction = 0.00038;
    let elevation_factor = solar_elevation.to_radians().sin().powf(0.8).max(0.15);

    let uv_base = clearsky_ghi * base_fraction * elevation_factor;
    let uv_ozone = uv_base * uv_transmission;

    // Aerosol correction (Ångström scaling)
    let aod_uv = angstrom_aod(0.1, 550.0, 310.0, 1.3);
    let aerosol_trans = (-0.5 * aod_uv * airmass).exp();

    let uv_index = (uv_ozone * aerosol_trans * 1.02) / 0.025;
    uv_index.clamp(0.0, 16.0)
}

/// Calculate UV transmission through ozone layer using empirically validated coefficients.
fn calculate_uv_transmission(solar_zenith_angle: f64, ozone_du: f64, airmass: f64) -> f64 {
    if solar_zenith_angle >= 90.0 {
        return 0.0;
    }

    // Slant path correction
    let ozone_path = ozone_du * airmass;

    // UV-B (280-315nm) - 0.08 per 100 DU
    let k_uvb = 0.08;
    let transmission_uvb = (-k_uvb * ozone_path / 100.0).exp();

    // UV-A (315-400nm) - 0.01 per 100 DU
    let k_uva = 0.01;
    let transmission_uva = (-k_uva * ozone_path / 100.0).exp();

    // Erythemal action spectrum weighting (UV-B contributes ~75%)
    0.75 * transmission_uvb + 0.25 * transmission_uva
}

#[cfg(test)]
mod tests {
    use super::*;
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
            let uvi_tropics = calculate_uv_index(90.0, 1050.0, 1.0, 0.0, 3, 0.0);
            assert!(uvi_tropics > 11.0, "Tropical noon UVI should be extreme, got {}", uvi_tropics);

            // Benchmark 2: Mid-latitude Summer Noon (45N, ~68 elev, June/Month 6)
            // Result should be "Very High" (8-10)
            let uvi_summer = calculate_uv_index(68.0, 950.0, 1.08, 45.0, 6, 0.0);
            assert!(
                uvi_summer > 7.0 && uvi_summer < 11.0,
                "Summer UVI {} out of range",
                uvi_summer
            );

            // Benchmark 3: Mid-latitude Winter Noon (45N, ~22 elev, Dec/Month 12)
            // Result should be "Low" (< 2)
            let uvi_winter = calculate_uv_index(22.0, 400.0, 2.6, 45.0, 12, 0.0);
            assert!(uvi_winter < 2.0, "Winter UVI {} must be low (< 2.0)", uvi_winter);
        }

        #[test]
        fn test_stratospheric_uv_resilience() {
            // Sea level: Month 6, altitude 0m
            let uvi_sea = calculate_uv_index(60.0, 800.0, 1.15, 45.0, 6, 0.0);
            // Stratosphere (11km): Less ozone above city
            let uvi_11km = calculate_uv_index(60.0, 1000.0, 0.25, 45.0, 6, 11000.0);

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
            let ghi = 750.0;
            let airmass = 1.39;

            let uvi = calculate_uv_index(sun_elevation, ghi, airmass, latitude, month, 200.0);

            // Expected: Moderate (3-6)
            assert!(uvi >= 3.0 && uvi <= 6.0, "Lapland summer UVI {} should be moderate", uvi);
        }

        #[test]
        fn test_extreme_uv_high_altitude_quito() {
            // Quito (~0.2°S) at ~2850m elevation in March (Month 3)
            let latitude = -0.2;
            let month = 3;
            let sun_elevation = 90.0;
            let ghi = 1150.0;
            let airmass = 0.72;

            let uvi = calculate_uv_index(sun_elevation, ghi, airmass, latitude, month, 2850.0);

            // Quito has extreme UV (>14) due to low ozone and altitude correction
            assert!(uvi > 14.0, "Quito equinox UVI {} should be extreme (>14)", uvi);
        }

        #[test]
        fn test_australia_summer_perth() {
            // Perth (~31.9°S) in January (Month 1)
            let latitude = -31.9;
            let month = 1;
            let sun_elevation = 81.5;
            let ghi = 1080.0;
            let airmass = 1.01;

            let uvi = calculate_uv_index(sun_elevation, ghi, airmass, latitude, month, 30.0);

            // Australian summer UVI is notoriously extreme (11+)
            assert!(uvi >= 11.0, "Australian summer UVI {} should be extreme (>=11)", uvi);
        }

        #[test]
        fn test_zero_uvi_at_night() {
            assert_eq!(calculate_uv_index(-5.0, 0.0, 0.0, 10.0, 6, 0.0), 0.0);
            assert_eq!(calculate_uv_index(0.0, 10.0, 10.0, 40.0, 6, 0.0), 0.0);
        }
    }
}
