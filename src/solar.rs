//! Solar Position Calculation Module
//!
//! Provides the core solar calculation context and position solving algorithms.
//! Uses the NREL SPA (Solar Position Algorithm) for high-precision calculations.

use chrono::{DateTime, Duration, TimeZone};
use chrono_tz::Tz;
use solar_positioning::{
    Horizon, spa,
    types::{RefractionCorrection, SunriseResult},
};

// ===================== TYPES =====================

/// Time and azimuth of a sun event (sunrise/sunset)
pub type SunEvent = (DateTime<Tz>, f64);

/// Sunrise and sunset pair
pub type SunEvents = (Option<SunEvent>, Option<SunEvent>);

// ===================== SOLAR CALCULATION CONTEXT =====================

/// Context for solar position calculations.
///
/// Encapsulates all parameters needed for consistent solar calculations
/// including observer location, altitude, atmospheric correction, and
/// target elevation angle.
#[derive(Clone, Copy)]
pub struct SolarCalc {
    /// Observer latitude in degrees
    pub lat: f64,
    /// Observer longitude in degrees
    pub lon: f64,
    /// Observer altitude in meters
    pub alt: f64,
    /// Delta-T correction for TT-UT1 difference
    pub delta_t: f64,
    /// Atmospheric refraction correction
    pub refr: Option<RefractionCorrection>,
    /// Target sun elevation angle in degrees
    pub target: f64,
}

impl SolarCalc {
    /// Calculate elevation error at a given time.
    ///
    /// Returns the difference between current sun elevation and target.
    pub fn elevation_error(&self, t: DateTime<Tz>) -> f64 {
        spa::solar_position(t, self.lat, self.lon, self.alt, self.delta_t, self.refr)
            .unwrap()
            .elevation_angle()
            - self.target
    }

    /// Get the solar position at a given time.
    pub fn position(&self, t: DateTime<Tz>) -> solar_positioning::SolarPosition {
        spa::solar_position(t, self.lat, self.lon, self.alt, self.delta_t, self.refr).unwrap()
    }

    /// Solve for the time when sun crosses target elevation using bisection.
    ///
    /// # Arguments
    /// * `a` - Start of search interval
    /// * `b` - End of search interval
    ///
    /// # Returns
    /// The crossing time and azimuth, or None if no crossing exists
    pub fn solve_root(&self, mut a: DateTime<Tz>, mut b: DateTime<Tz>) -> Option<SunEvent> {
        let mut fa = self.elevation_error(a);
        let fb = self.elevation_error(b);

        // Guard against NaN from invalid inputs
        if !fa.is_finite() || !fb.is_finite() {
            return None;
        }

        if fa.signum() == fb.signum() {
            return None;
        }

        for _ in 0..60 {
            let m = a + (b - a) / 2;
            let fm = self.elevation_error(m);

            if !fm.is_finite() {
                return None;
            }

            if fm.abs() < 1e-7 {
                let az = self.position(m).azimuth();
                return Some((m, az));
            }

            if fm.signum() == fa.signum() {
                a = m;
                fa = fm;
            } else {
                b = m;
            }
        }

        let az = self.position(a).azimuth();
        Some((a, az))
    }

    /// Solve for sunrise and sunset from solar noon.
    ///
    /// # Arguments
    /// * `noon` - Solar transit time
    ///
    /// # Returns
    /// Tuple of (sunrise, sunset) events
    pub fn solve_from_noon(&self, noon: DateTime<Tz>) -> SunEvents {
        let span = Duration::hours(12);
        (self.solve_root(noon - span, noon), self.solve_root(noon, noon + span))
    }

    /// Get the solar transit result for a given date.
    ///
    /// # Arguments
    /// * `date` - Date to calculate transit for
    ///
    /// # Returns
    /// SunriseResult containing transit information
    pub fn get_transit(&self, date: DateTime<Tz>) -> Option<SunriseResult<DateTime<Tz>>> {
        spa::sunrise_sunset_for_horizon(
            date,
            self.lat,
            self.lon,
            self.delta_t,
            Horizon::SunriseSunset,
        )
        .ok()
    }

    /// Extract the transit time from a SunriseResult.
    pub fn extract_transit_time(&self, res: &SunriseResult<DateTime<Tz>>) -> DateTime<Tz> {
        match res {
            SunriseResult::RegularDay { transit, .. } => *transit,
            SunriseResult::AllDay { transit } => *transit,
            SunriseResult::AllNight { transit } => *transit,
        }
    }

    /// Find the next sunrise or sunset event after a given time.
    ///
    /// Searches up to 370 days ahead for polar regions.
    ///
    /// # Arguments
    /// * `start` - Time to start searching from
    ///
    /// # Returns
    /// Tuple of (event_name, event_time), or None if not found
    pub fn find_next_event(&self, start: DateTime<Tz>) -> Option<(String, DateTime<Tz>)> {
        let tz = start.timezone();
        let mut current_naive = start.date_naive();

        for _ in 0..370 {
            // Safe Anchor: Find a valid time on this calendar day.
            // We prefer 12:00:00 as it is the most stable reference for solar events.
            let d = match tz.from_local_datetime(&current_naive.and_hms_opt(12, 0, 0)?) {
                chrono::LocalResult::Single(t) => t,
                chrono::LocalResult::Ambiguous(t, _) => t,
                chrono::LocalResult::None => {
                    // If 12:00 doesn't exist (very rare), try start of day logic
                    tz.from_local_datetime(&current_naive.and_hms_opt(0, 0, 0)?)
                        .earliest()
                        .unwrap_or(start)
                }
            };

            if let Some(transit_res) = self.get_transit(d) {
                let transit = self.extract_transit_time(&transit_res);
                let (sr, ss) = self.solve_from_noon(transit);

                let mut events = Vec::new();
                if let Some((t, _)) = sr
                    && t > start
                {
                    events.push(("Sunrise", t));
                }
                if let Some((t, _)) = ss
                    && t > start
                {
                    events.push(("Sunset", t));
                }

                if let Some((kind, t)) = events.into_iter().min_by_key(|(_, t)| *t) {
                    return Some((kind.into(), t));
                }
            }

            // Move to the next calendar day safely
            current_naive = current_naive.succ_opt()?;
        }
        None
    }
}

// ===================== HELPER FUNCTIONS =====================

/// Calculate day length from sunrise and sunset events.
///
/// # Arguments
/// * `sr` - Sunrise event (time and azimuth)
/// * `ss` - Sunset event (time and azimuth)
///
/// # Returns
/// Day length in seconds, or None if either event is missing
pub fn day_length(sr: &Option<SunEvent>, ss: &Option<SunEvent>) -> Option<i64> {
    Some((ss.as_ref()?.0 - sr.as_ref()?.0).num_seconds())
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SOLAR_RADIUS_DEG;
    use crate::geo::horizon_dip_deg;
    use chrono::TimeZone;
    use chrono_tz::Europe::Helsinki;
    use chrono_tz::UTC;
    use solar_positioning::time::DeltaT;

    #[test]
    fn test_dead_sea_sunrise_shift() {
        let tz = Helsinki;
        let date = tz.with_ymd_and_hms(2025, 7, 1, 0, 0, 0).unwrap();
        let lat = 31.0;
        let lon = 35.4;

        let delta_t: f64 = DeltaT::estimate_from_date(2025, 7).unwrap();
        let refraction = Some(RefractionCorrection::standard());

        let base_alt = -SOLAR_RADIUS_DEG;

        // --- Sea level ---
        let calc0 = SolarCalc {
            lat,
            lon,
            alt: 0.0,
            delta_t,
            refr: refraction,
            target: base_alt - horizon_dip_deg(lat, 0.0),
        };

        let noon0_res = calc0.get_transit(date).unwrap();
        let noon0 = calc0.extract_transit_time(&noon0_res);
        let (sr0, _) = calc0.solve_from_noon(noon0);
        let sr0 = sr0.unwrap().0;

        // --- Dead Sea ---
        let calc_ds = SolarCalc {
            lat,
            lon,
            alt: -450.0,
            delta_t,
            refr: refraction,
            target: base_alt - horizon_dip_deg(lat, -450.0),
        };

        let (sr_ds, _) = calc_ds.solve_from_noon(noon0);
        let sr_ds = sr_ds.unwrap().0;

        let shift = (sr_ds - sr0).num_seconds();

        // NOAA + JPL both show ~3–5 minutes later sunrise
        assert!(
            shift >= 180 && shift <= 300,
            "Dead Sea sunrise shift out of expected range: {} s",
            shift
        );
    }

    #[test]
    fn test_civil_twilight_ignores_horizon_dip() {
        let tz = Helsinki;
        let lat = 31.0;
        let lon = 35.4;
        let delta_t = DeltaT::estimate_from_date(2025, 7).unwrap();
        let refraction = Some(RefractionCorrection::standard());

        let date = tz.with_ymd_and_hms(2025, 7, 1, 0, 0, 0).unwrap();

        // Civil twilight altitude
        let target_alt = -6.0;

        // Get solar noon (same for both)
        let res = spa::sunrise_sunset_for_horizon(date, lat, lon, delta_t, Horizon::CivilTwilight)
            .unwrap();

        let noon = match res {
            SunriseResult::RegularDay { transit, .. } => transit,
            _ => panic!("Expected civil twilight on this date"),
        };

        // Low altitude observer
        let calc_low =
            SolarCalc { lat, lon, alt: 0.0, delta_t, refr: refraction, target: target_alt };
        let (tw_low, _) = calc_low.solve_from_noon(noon);
        let tw_low = tw_low.unwrap().0;

        // High altitude observer
        let calc_high =
            SolarCalc { lat, lon, alt: 2000.0, delta_t, refr: refraction, target: target_alt };
        let (tw_high, _) = calc_high.solve_from_noon(noon);
        let tw_high = tw_high.unwrap().0;

        let diff = (tw_high - tw_low).num_seconds().abs();

        // Twilight must not shift with altitude
        assert!(diff <= 1, "Civil twilight shifted by {} seconds due to altitude", diff);
    }

    #[test]
    fn test_midnight_sun_exists_in_may_tromso() {
        use chrono_tz::Europe::Oslo;

        let tz = Oslo;
        let lat = 69.6492;
        let lon = 18.9553;
        let delta_t = DeltaT::estimate_from_date(2025, 5).unwrap();

        let mut found_midnight_sun = false;

        // Scan a safe window
        for day in 15..25 {
            let date = tz.with_ymd_and_hms(2025, 5, day, 0, 0, 0).unwrap();

            let res =
                spa::sunrise_sunset_for_horizon(date, lat, lon, delta_t, Horizon::SunriseSunset)
                    .unwrap();

            if matches!(res, SunriseResult::AllDay { .. }) {
                found_midnight_sun = true;
                break;
            }
        }

        assert!(found_midnight_sun, "Expected at least one midnight-sun day in mid-May");
    }

    #[test]
    fn test_equator_never_has_polar_night() {
        let tz = UTC;
        let lat = 0.0;
        let lon = 0.0;
        let year = 2025;

        let refraction = Some(RefractionCorrection::standard());
        let target_alt = -SOLAR_RADIUS_DEG;

        for month in 1..=12 {
            for day in 1..=31 {
                // Skip invalid dates safely
                let date = match tz.with_ymd_and_hms(year, month, day, 0, 0, 0) {
                    chrono::LocalResult::Single(d) => d,
                    _ => continue, // Invalid or ambiguous date (DST etc.)
                };

                let delta_t = DeltaT::estimate_from_date(year, month).unwrap();

                let calc =
                    SolarCalc { lat, lon, alt: 0.0, delta_t, refr: refraction, target: target_alt };

                let noon_res = calc.get_transit(date).expect("Failed to get transit at equator");
                let noon = calc.extract_transit_time(&noon_res);

                let (sr, ss) = calc.solve_from_noon(noon);

                assert!(
                    sr.is_some() && ss.is_some(),
                    "Missing sunrise/sunset at equator on {}",
                    date.date_naive()
                );
            }
        }
    }

    #[test]
    fn test_day_length_continuity_across_centuries() {
        let lat = 45.0;
        let lon = 0.0;
        let refraction = Some(RefractionCorrection::standard());
        let target_alt = -SOLAR_RADIUS_DEG;

        // ΔT model is valid from year -500 onward
        let years = [-500, -200, 0, 500, 1000, 1800, 1900, 2000, 2001, 2100, 3000, 5000];

        for &year in &years {
            let date = UTC.with_ymd_and_hms(year, 6, 21, 0, 0, 0).unwrap();

            let delta_t = match DeltaT::estimate_from_date(year, 6) {
                Ok(dt) => dt,
                Err(_) => continue, // Explicitly respect ΔT limits
            };

            let calc =
                SolarCalc { lat, lon, alt: 0.0, delta_t, refr: refraction, target: target_alt };

            let noon_res = match calc.get_transit(date) {
                Some(n) => n,
                None => continue,
            };
            let noon = calc.extract_transit_time(&noon_res);

            let (sr, ss) = calc.solve_from_noon(noon);

            let sr = sr.unwrap().0;
            let ss = ss.unwrap().0;

            let day_len = (ss - sr).num_seconds();

            assert!(day_len > 0);
            assert!(day_len < 24 * 3600);
        }
    }
}
