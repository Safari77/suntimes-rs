//! Geographic and Geometry Module
//!
//! Provides WGS84 Earth radius calculations and horizon dip angle computations
//! for accurate sunrise/sunset calculations at different altitudes.

// ===================== CONSTANTS =====================

/// WGS84 semi-major axis (equatorial radius) in meters
pub const A_EQUATOR: f64 = 6_378_137.0;

/// WGS84 semi-minor axis (polar radius) in meters
pub const B_POLAR: f64 = 6_356_752.314245;

/// Solar apparent radius in degrees (angular semi-diameter)
pub const SOLAR_RADIUS_DEG: f64 = 0.266;

// ===================== GEOMETRY FUNCTIONS =====================

/// Calculate Earth radius at a given latitude using WGS84 ellipsoid model.
///
/// # Arguments
/// * `lat_deg` - Latitude in degrees (-90 to 90)
///
/// # Returns
/// Earth radius in meters at the specified latitude
pub fn earth_radius_wgs84(lat_deg: f64) -> f64 {
    let phi = lat_deg.to_radians();
    let cos = phi.cos();
    let sin = phi.sin();
    let a2 = A_EQUATOR * A_EQUATOR;
    let b2 = B_POLAR * B_POLAR;
    let numerator = a2 * a2 * cos * cos + b2 * b2 * sin * sin;
    let denominator = (A_EQUATOR * cos).powi(2) + (B_POLAR * sin).powi(2);
    (numerator / denominator).sqrt()
}

/// Calculate horizon dip angle due to observer altitude.
///
/// The dip angle represents how far below the geometric horizon
/// the apparent horizon appears due to the observer's elevation.
///
/// Positive altitude (above sea level) gives positive dip (horizon below geometric level).
/// Negative altitude (below sea level, e.g., Dead Sea) gives negative dip.
///
/// # Arguments
/// * `lat_deg` - Observer latitude in degrees
/// * `h` - Observer altitude in meters (can be negative)
///
/// # Returns
/// Horizon dip angle in degrees
pub fn horizon_dip_deg(lat_deg: f64, h: f64) -> f64 {
    if h.abs() < 1e-5 {
        return 0.0;
    }
    let r = earth_radius_wgs84(lat_deg);
    // Use magnitude for geometry, apply sign afterward
    let ratio = r / (r + h.abs());
    let clamped = ratio.clamp(-1.0, 1.0);
    let dip = clamped.acos().to_degrees();
    if h > 0.0 { dip } else { -dip }
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_earth_radius_wgs84_reasonable() {
        let r_equator = earth_radius_wgs84(0.0);
        let r_pole = earth_radius_wgs84(90.0);

        // Expected magnitudes
        assert!(r_equator > 6_370_000.0 && r_equator < 6_380_000.0);
        assert!(r_pole > 6_350_000.0 && r_pole < 6_360_000.0);

        // Equatorial radius must be larger
        assert!(r_equator > r_pole);
    }

    #[test]
    fn test_horizon_dip_sign_and_magnitude() {
        let lat = 31.0;

        let dip_0 = horizon_dip_deg(lat, 0.0);
        let dip_pos = horizon_dip_deg(lat, 1000.0);
        let dip_neg = horizon_dip_deg(lat, -450.0);

        // Zero altitude → zero dip
        assert!(dip_0.abs() < 1e-6);

        // Positive altitude → positive dip
        assert!(dip_pos > 0.5 && dip_pos < 1.2);

        // Negative altitude → negative dip
        assert!(dip_neg < -0.4 && dip_neg > -0.9);
    }
}
