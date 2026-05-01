//! UV Index Calculation Module
//!
//! Based on empirical models for Ozone column estimation and UV transmission.

/// Global monthly ozone climatology (Dobson Units) dataset based on modern OMI/TOMS satellite data
/// Rows: 9 latitude band centers from -80° to +80° (every 20°);
///       values represent the band-center monthly mean and are bilinearly interpolated.
/// Cols: Months Jan to Dec (interpreted as monthly mean / mid-month)
const OZONE_GRID: [[f64; 12]; 9] = [
    // Center -80° (covers -90 to -70, South Pole)
    [280.0, 270.0, 260.0, 250.0, 250.0, 250.0, 250.0, 220.0, 150.0, 140.0, 180.0, 250.0],
    // Center -60°
    [300.0, 290.0, 290.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 330.0, 320.0, 310.0],
    // Center -40°
    [280.0, 275.0, 275.0, 280.0, 290.0, 310.0, 330.0, 340.0, 350.0, 340.0, 320.0, 295.0],
    // Center -20° (Subtropical SH)
    [265.0, 260.0, 260.0, 260.0, 265.0, 270.0, 280.0, 285.0, 290.0, 285.0, 275.0, 270.0],
    // Center 0° (Equatorial)
    [250.0, 255.0, 260.0, 260.0, 260.0, 255.0, 255.0, 260.0, 260.0, 260.0, 255.0, 250.0],
    // Center +20° (Subtropical NH)
    [265.0, 275.0, 285.0, 295.0, 300.0, 295.0, 285.0, 280.0, 275.0, 265.0, 260.0, 260.0],
    // Center +40° (Mid-latitude NH)
    [320.0, 340.0, 360.0, 365.0, 340.0, 315.0, 300.0, 290.0, 285.0, 285.0, 295.0, 310.0],
    // Center +60°
    [360.0, 380.0, 400.0, 395.0, 365.0, 335.0, 315.0, 300.0, 290.0, 290.0, 310.0, 330.0],
    // Center +80° (covers +70 to +90, North Pole)
    [380.0, 400.0, 420.0, 410.0, 380.0, 345.0, 320.0, 300.0, 290.0, 290.0, 310.0, 340.0],
];

/// Bilinear interpolation across the ozone climatology grid.
///
/// `frac_month` is a continuous column-index where integer values land exactly on
/// table columns: 0.0 = January's monthly mean, 5.0 = June's monthly mean. Halfway
/// values (e.g. 5.5) interpolate between adjacent monthly means (chronologically:
/// the boundary between those months). Wraps cyclically (11.5 = Dec → Jan).
fn lookup_ozone_bilinear(latitude: f64, frac_month: f64) -> f64 {
    // Latitude → fractional row. Band centers are at -80, -60, ..., +80 (rows 0..8);
    // outside ±80 we clamp to the polar band.
    let frac_row = ((latitude + 80.0) / 20.0).clamp(0.0, 8.0);
    let row0 = frac_row.floor() as usize;
    let row1 = (row0 + 1).min(8);
    let row_t = frac_row - row0 as f64;

    // Month → fractional column with cyclic wrap (Dec ↔ Jan).
    let fc = frac_month.rem_euclid(12.0);
    let col0 = fc.floor() as usize;
    let col1 = (col0 + 1) % 12;
    let col_t = fc - col0 as f64;

    let v00 = OZONE_GRID[row0][col0];
    let v01 = OZONE_GRID[row0][col1];
    let v10 = OZONE_GRID[row1][col0];
    let v11 = OZONE_GRID[row1][col1];

    // Interpolate along month first, then latitude.
    let v0 = v00 * (1.0 - col_t) + v01 * col_t;
    let v1 = v10 * (1.0 - col_t) + v11 * col_t;
    v0 * (1.0 - row_t) + v1 * row_t
}

/// Estimate total ozone column in Dobson Units (DU)
///
/// Uses bilinear interpolation over a monthly gridded climatology, with a
/// tropospheric altitude correction on top.
///
/// Note: This function is called from inner loops in yearly optimization, so
/// it deliberately performs no I/O. Callers that want to surface the
/// ozone-hole caveat should check [`is_ozone_hole_season`] once at the top
/// level and print their own message.
pub fn estimate_ozone_column(latitude: f64, month: u32, altitude_m: f64) -> f64 {
    // 1. Map integer month to its column (mid-month value); latitude is interpolated
    //    bilinearly against neighboring band centers. Callers wanting smooth time
    //    interpolation can drive lookup_ozone_bilinear with a fractional month.
    let frac_month = (month.clamp(1, 12) - 1) as f64;
    let mut ozone_du = lookup_ozone_bilinear(latitude, frac_month);

    // 2. Tropospheric Altitude Correction.
    // ~90% of the ozone column is stratospheric, so a city at altitude only loses
    // the small tropospheric fraction above it (~10%). A rate of ~1%/km matches
    // that. Clamped at 11 km (tropopause). The previous 3.5%/km double-counted
    // with the +%/km altitude UV boost in calculate_uv_index.
    let altitude_km = (altitude_m / 1000.0).clamp(0.0, 11.0);
    let elevation_correction = 1.0 - (0.010 * altitude_km);
    ozone_du *= elevation_correction;

    // Floor at 100 DU lets the climatology's deepest Antarctic spring values
    // (e.g. 140 DU October at 80°S) through; the previous 150 DU floor was
    // suppressing them.
    ozone_du.clamp(100.0, 500.0)
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

/// Aerosol-specific relative airmass (Kasten 1965).
///
/// Aerosols are concentrated in the boundary layer, so their slant-path airmass
/// differs from the standard atmospheric airmass — most noticeably at low solar
/// elevations. Caller must ensure `solar_elevation_deg > 0`.
///
/// Reference: Kasten F., Archiv für Meteorologie, Geophysik und Bioklimatologie B,
/// 14, 206-223 (1965).
fn aerosol_airmass(solar_elevation_deg: f64) -> f64 {
    let h = solar_elevation_deg;
    1.0 / (h.to_radians().sin() + 0.0548 * (h + 2.65).powf(-1.452))
}

/// Calculates clear-sky UV Index using an empirical Solar Zenith Angle parameterization adjusted
/// for ozone, aerosols, altitude, and albedo.
pub fn calculate_uv_index(solar_elevation: f64, latitude: f64, month: u32, altitude_m: f64) -> f64 {
    if solar_elevation <= 0.0 {
        return 0.0;
    }

    let ozone_du = estimate_ozone_column(latitude, month, altitude_m);
    let altitude_km = altitude_m / 1000.0;

    // Solar zenith angle conversion
    let mu = solar_elevation.to_radians().sin();

    // Standard clear-sky empirical UVI model (mu^2.42 bundles geometric
    // dependence including approximate ozone slant path scaling).
    let uv_ozone = 12.5 * mu.powf(2.42) * (300.0 / ozone_du).powf(1.2);

    // Aerosol attenuation. Two physical effects with altitude:
    //   1. Aerosol-specific airmass (Kasten 1965), more accurate than the
    //      generic atmospheric airmass especially at low sun.
    //   2. AOD itself decreases roughly exponentially with altitude (boundary
    //      -layer aerosol scale height ~1.5 km), since most aerosols sit
    //      below the observer.
    // The 0.5 prefactor is the conventional UV "effective AOD" reduction
    // accounting for the diffuse forward-scattered contribution that still
    // reaches the surface.
    const AEROSOL_SCALE_HEIGHT_KM: f64 = 1.5;
    let aod_uv_sea_level = angstrom_aod(0.04, 550.0, 310.0, 1.3);
    let aod_uv = aod_uv_sea_level * (-altitude_km / AEROSOL_SCALE_HEIGHT_KM).exp();
    let m_aerosol = aerosol_airmass(solar_elevation);
    let aerosol_trans = (-0.5 * aod_uv * m_aerosol).exp();

    let albedo_factor = 1.02;

    // Altitude correction for residual physics not captured by the explicit
    // ozone and aerosol terms above (mainly reduced Rayleigh scattering above
    // the observer). Softened to ~5%/km from the previous +10%/km, which had
    // been double-counting the ozone and aerosol altitude effects now modeled
    // separately above.
    let altitude_correction = 1.0 + (0.05 * altitude_km);

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
        let uvi_tropics = calculate_uv_index(90.0, 0.0, 3, 0.0);
        assert!(uvi_tropics > 11.0, "Tropical noon UVI should be extreme, got {}", uvi_tropics);

        // Benchmark 2: Mid-latitude Summer Noon (45N, ~68 elev, June/Month 6)
        // Result should be "Very High" (8-10)
        let uvi_summer = calculate_uv_index(68.0, 45.0, 6, 0.0);
        assert!(uvi_summer > 7.0 && uvi_summer < 11.0, "Summer UVI {} out of range", uvi_summer);

        // Benchmark 3: Mid-latitude Winter Noon (45N, ~22 elev, Dec/Month 12)
        // Result should be "Low" (< 2)
        let uvi_winter = calculate_uv_index(22.0, 45.0, 12, 0.0);
        assert!(uvi_winter < 2.0, "Winter UVI {} must be low (< 2.0)", uvi_winter);
    }

    #[test]
    fn test_stratospheric_uv_resilience() {
        // Sea level: Month 6, altitude 0m
        let uvi_sea = calculate_uv_index(60.0, 45.0, 6, 0.0);
        // Stratosphere (11km): Less ozone above city
        let uvi_11km = calculate_uv_index(60.0, 45.0, 6, 11000.0);

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

        let uvi = calculate_uv_index(sun_elevation, latitude, month, 200.0);

        // Expected: Moderate (3-6)
        assert!(uvi >= 3.0 && uvi <= 6.0, "Lapland summer UVI {} should be moderate", uvi);
    }

    #[test]
    fn test_extreme_uv_high_altitude_quito() {
        // Quito (~0.2°S) at ~2850m elevation in March (Month 3)
        let latitude = -0.2;
        let month = 3;
        let sun_elevation = 90.0;

        let uvi = calculate_uv_index(sun_elevation, latitude, month, 2850.0);

        // Quito has extreme UV (>14) due to low ozone and altitude correction
        assert!(uvi > 14.0, "Quito equinox UVI {} should be extreme (>14)", uvi);
    }

    #[test]
    fn test_australia_summer_perth() {
        // Perth (~31.9°S) in January (Month 1)
        let latitude = -31.9;
        let month = 1;
        let sun_elevation = 81.5;

        let uvi = calculate_uv_index(sun_elevation, latitude, month, 30.0);

        // Australian summer UVI is notoriously extreme (11+)
        assert!(uvi >= 11.0, "Australian summer UVI {} should be extreme (>=11)", uvi);
    }

    #[test]
    fn test_zero_uvi_at_night() {
        assert_eq!(calculate_uv_index(-5.0, 10.0, 6, 0.0), 0.0);
        assert_eq!(calculate_uv_index(0.0, 40.0, 6, 0.0), 0.0);
    }

    #[test]
    fn test_latitude_continuity() {
        // Bilinear interpolation should produce a smooth gradient across the old
        // 20°-band boundaries. With the previous floor() lookup, lat 49.9 vs 50.1
        // jumped from row 6 to row 7 → ~25 DU step in summer. Now should be smooth.
        let a = estimate_ozone_column(49.9, 6, 0.0);
        let b = estimate_ozone_column(50.1, 6, 0.0);
        assert!((a - b).abs() < 1.0, "Latitude transition not smooth: {} vs {}", a, b);
    }

    #[test]
    fn test_antarctic_ozone_hole_value_not_clamped() {
        // 80°S, October — the climatology has 140 DU here. The previous clamp at
        // 150 DU was incorrectly raising it; the new 100 DU floor lets it through.
        let ozone = estimate_ozone_column(-80.0, 10, 0.0);
        assert!(
            ozone < 145.0,
            "Antarctic October ozone {} should reflect 140 DU table value",
            ozone
        );
    }
}
