//! Solar Position Calculation Module
//!
//! Provides the core solar calculation context and position solving algorithms.
//! Uses the NREL SPA (Solar Position Algorithm) for high-precision calculations.

use chrono::{DateTime, Datelike, Duration, NaiveDate, TimeZone, Timelike};
use chrono_tz::Tz;
use solar_positioning::{
    Horizon,
    error::Error as SpaError,
    spa,
    types::{RefractionCorrection, SunriseResult},
};

// ===================== TYPES =====================

/// Time and azimuth of a sun event (sunrise/sunset)
pub type SunEvent = (DateTime<Tz>, f64);

/// Sunrise and sunset pair
pub type SunEvents = (Option<SunEvent>, Option<SunEvent>);

/// One yearly extreme: the event itself plus the clock reading used to rank it.
#[derive(Clone, Copy, Debug)]
pub struct ExtremeEvent {
    /// Time of the event in the observer's time zone
    pub time: DateTime<Tz>,
    /// Azimuth of the sun at the event, in degrees
    pub azimuth: f64,
    /// Where the event sits on the local clock face, in seconds counted from
    /// 00:00 of the calendar day it was calculated for.
    ///
    /// Normally within `0..86_400`, and then exactly equal to the wall-clock
    /// reading — including across DST changes, since it is derived from the
    /// local time fields and not from an elapsed-time subtraction.
    ///
    /// Close to the polar circles an event can spill over midnight: a sunrise
    /// at 23:50 the previous evening, or a sunset at 00:30 the next morning.
    /// Those get negative values, or values above 86_400, so that the ordering
    /// stays continuous instead of wrapping around and reporting a sunset just
    /// after midnight as the "earliest" of the year.
    pub clock_seconds: i64,
}

/// The four extreme sunrise/sunset times of a calendar year, plus statistics
/// about the scan that produced them.
#[derive(Clone, Copy, Debug, Default)]
pub struct YearlyExtremes {
    /// Earliest sunrise by local clock time
    pub earliest_sunrise: Option<ExtremeEvent>,
    /// Latest sunrise by local clock time
    pub latest_sunrise: Option<ExtremeEvent>,
    /// Earliest sunset by local clock time
    pub earliest_sunset: Option<ExtremeEvent>,
    /// Latest sunset by local clock time
    pub latest_sunset: Option<ExtremeEvent>,
    /// Calendar days actually scanned: 365 (366 in a leap year), minus any
    /// day that does not exist in this time zone
    pub days_scanned: u32,
    /// Days on which the sun crossed the target elevation going up
    pub days_with_sunrise: u32,
    /// Days on which the sun crossed the target elevation going down
    pub days_with_sunset: u32,
}

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
    /// Compute the solar position and the elevation residual (`elevation - target`)
    /// in a single SPA evaluation.
    ///
    /// This is the fundamental cached primitive used by the bisection search —
    /// we call into SPA exactly once per time point and reuse the result for both
    /// the residual test and (at convergence) the azimuth readout.
    fn position_and_error(
        &self,
        t: DateTime<Tz>,
    ) -> Result<(solar_positioning::SolarPosition, f64), SpaError> {
        let pos = spa::solar_position(t, self.lat, self.lon, self.alt, self.delta_t, self.refr)?;
        let err = pos.elevation_angle() - self.target;
        Ok((pos, err))
    }

    /// Get the solar position at a given time.
    pub fn position(&self, t: DateTime<Tz>) -> Result<solar_positioning::SolarPosition, SpaError> {
        spa::solar_position(t, self.lat, self.lon, self.alt, self.delta_t, self.refr)
    }

    /// Solve for the time when sun crosses target elevation using bisection.
    ///
    /// # Arguments
    /// * `a` - Start of search interval
    /// * `b` - End of search interval
    ///
    /// # Returns
    /// The crossing time and azimuth, or `Ok(None)` if no crossing exists in
    /// the interval. Returns `Err` only on underlying SPA failure (invalid
    /// inputs or out-of-range date).
    pub fn solve_root(
        &self,
        mut a: DateTime<Tz>,
        mut b: DateTime<Tz>,
    ) -> Result<Option<SunEvent>, SpaError> {
        // Cache both endpoints: we need `pa`'s azimuth at the final fallback,
        // and `pb`'s azimuth if the endpoint itself is already the root.
        let (mut pa, mut fa) = self.position_and_error(a)?;
        let (pb, fb) = self.position_and_error(b)?;

        // Guard against NaN from invalid inputs
        if !fa.is_finite() || !fb.is_finite() {
            return Ok(None);
        }

        // If endpoints are already legitimate crossing points, return them.
        if fa.abs() < 1e-7 {
            return Ok(Some((a, pa.azimuth())));
        }
        if fb.abs() < 1e-7 {
            return Ok(Some((b, pb.azimuth())));
        }

        if fa.signum() == fb.signum() {
            return Ok(None);
        }

        for _ in 0..60 {
            // chrono Durations are integer nanoseconds: once the bracket is a
            // single nanosecond wide the midpoint stops moving, and the
            // remaining iterations repeat the same SPA call for nothing.
            if b - a <= Duration::nanoseconds(1) {
                break;
            }

            let m = a + (b - a) / 2;
            // One SPA call per midpoint — no extra call at convergence.
            let (pm, fm) = self.position_and_error(m)?;

            if !fm.is_finite() {
                return Ok(None);
            }

            if fm.abs() < 1e-7 {
                return Ok(Some((m, pm.azimuth())));
            }

            if fm.signum() == fa.signum() {
                a = m;
                pa = pm; // carry the cached position forward
                fa = fm;
            } else {
                b = m;
            }
        }

        // Bisection did not converge within 60 iterations — return current best
        // using the cached position at `a`. No extra SPA call.
        Ok(Some((a, pa.azimuth())))
    }

    /// Solve for sunrise and sunset from solar noon.
    ///
    /// # Arguments
    /// * `noon` - Solar transit time
    ///
    /// # Returns
    /// Tuple of (sunrise, sunset) events
    pub fn solve_from_noon(&self, noon: DateTime<Tz>) -> Result<SunEvents, SpaError> {
        let span = Duration::hours(12);
        let sr = self.solve_root(noon - span, noon)?;
        let ss = self.solve_root(noon, noon + span)?;
        Ok((sr, ss))
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
    /// `Ok(Some((event_name, event_time)))` if found, `Ok(None)` if no event
    /// within 370 days, `Err` on SPA failure.
    pub fn find_next_event(
        &self,
        start: DateTime<Tz>,
    ) -> Result<Option<(String, DateTime<Tz>)>, SpaError> {
        let tz = start.timezone();
        let mut current_naive = start.date_naive();

        for _ in 0..370 {
            // Safe Anchor: Find a valid time on this calendar day.
            // We prefer 12:00:00 as it is the most stable reference for solar events.
            let d = match tz.from_local_datetime(&current_naive.and_hms_opt(12, 0, 0).unwrap()) {
                chrono::LocalResult::Single(t) => t,
                chrono::LocalResult::Ambiguous(t, _) => t,
                chrono::LocalResult::None => {
                    // If 12:00 doesn't exist (very rare), fallback to the first
                    // valid instant of the day using robust local-time probing.
                    match crate::time::start_of_day_opt(tz, current_naive) {
                        Some(t) => t,
                        None => start, // ultimate fallback if the entire day is somehow skipped
                    }
                }
            };

            if let Some(transit_res) = self.get_transit(d) {
                let transit = self.extract_transit_time(&transit_res);
                let (sr, ss) = self.solve_from_noon(transit)?;

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
                    return Ok(Some((kind.into(), t)));
                }
            }

            // Move to the next calendar day safely
            let Some(next) = current_naive.succ_opt() else {
                return Ok(None);
            };
            current_naive = next;
        }
        Ok(None)
    }

    /// Scan a whole calendar year for the extreme sunrise and sunset times.
    ///
    /// Every day is solved with exactly the same machinery the single-day
    /// output uses (`get_transit` + `solve_from_noon`), so the extremes always
    /// agree with what the program prints for those dates, and they follow the
    /// active target elevation — with `--civil` this yields the extremes of
    /// civil dawn and dusk rather than of sunrise and sunset.
    ///
    /// Events are ranked by their position on the local clock face, not by UTC
    /// instant: "earliest sunrise" is what an observer reading a clock means by
    /// it, and a DST change genuinely does move sunrise on that clock. See
    /// [`ExtremeEvent::clock_seconds`] for events that spill across midnight.
    ///
    /// A consequence worth knowing about: in a zone whose DST transition falls
    /// far enough from the solstice, an extreme can land on the transition
    /// boundary rather than near the solstice. In Sydney *both* sunrise
    /// extremes do — the last morning before the October change (05:28) and
    /// the last one before the April change (07:10) — because either side of
    /// those cliffs the clock reading moves by a whole hour, more than the
    /// seasonal swing can recover. The sunsets there behave seasonally. Both
    /// are the correct answer to "earliest and latest by the clock". In zones
    /// whose transitions sit near the equinoxes — the EU, for one — the
    /// question does not arise.
    ///
    /// The search covers 1 January to 31 December of `year` inclusive. In the
    /// southern hemisphere the summer extremes straddle New Year, so the latest
    /// sunset of a calendar year normally lands in early January and the
    /// earliest sunrise in early December.
    ///
    /// `delta_t` is held fixed at the value in `self` for the whole year. Its
    /// drift over twelve months is well under a second and applies almost
    /// equally to every day, so it cannot move which day comes out extreme.
    ///
    /// # Arguments
    /// * `tz` - Time zone whose clock the events are read on
    /// * `year` - Calendar year to scan
    ///
    /// # Returns
    /// The four extremes plus scan statistics. A field is `None` when that
    /// event never occurs during the year (polar day or polar night).
    ///
    /// # Errors
    /// Propagates an SPA failure from the underlying root solve. Days whose
    /// transit cannot be computed, and days that do not exist locally at all,
    /// are skipped rather than treated as errors.
    ///
    /// # Cost
    /// Roughly 45 000 SPA evaluations for a full year — fast in a release
    /// build, noticeably slower in a debug build. This is why it sits behind
    /// its own flag rather than running on every invocation.
    pub fn yearly_extremes(&self, tz: Tz, year: i32) -> Result<YearlyExtremes, SpaError> {
        let mut ex = YearlyExtremes::default();

        let Some(jan_first) = NaiveDate::from_ymd_opt(year, 1, 1) else {
            return Ok(ex); // year outside chrono's representable range
        };

        for day in jan_first.iter_days().take_while(|d| d.year() == year) {
            // Anchor on local noon, the most stable reference for the transit.
            let Some(anchor) = crate::time::noon_of_day_opt(tz, day) else {
                continue; // the whole calendar day is skipped in this zone
            };
            ex.days_scanned += 1;

            let Some(transit_res) = self.get_transit(anchor) else {
                continue; // SPA could not place the transit for this day
            };
            let transit = self.extract_transit_time(&transit_res);
            let (sr, ss) = self.solve_from_noon(transit)?;

            if let Some((t, az)) = sr {
                ex.days_with_sunrise += 1;
                let cand =
                    ExtremeEvent { time: t, azimuth: az, clock_seconds: clock_seconds(t, day) };
                keep_earlier(&mut ex.earliest_sunrise, cand);
                keep_later(&mut ex.latest_sunrise, cand);
            }

            if let Some((t, az)) = ss {
                ex.days_with_sunset += 1;
                let cand =
                    ExtremeEvent { time: t, azimuth: az, clock_seconds: clock_seconds(t, day) };
                keep_earlier(&mut ex.earliest_sunset, cand);
                keep_later(&mut ex.latest_sunset, cand);
            }
        }

        Ok(ex)
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

/// Position of `t` on the local clock face, in seconds from 00:00 of `day`.
///
/// Built from the local date and time fields rather than from an elapsed-time
/// subtraction, so on a DST transition day the result is still the reading a
/// clock on the wall would show. Events landing on a neighbouring calendar day
/// are carried by whole days, giving values outside `0..86_400` — see
/// [`ExtremeEvent::clock_seconds`].
fn clock_seconds(t: DateTime<Tz>, day: NaiveDate) -> i64 {
    let day_shift = (t.date_naive() - day).num_days();
    day_shift * 86_400 + i64::from(t.num_seconds_from_midnight())
}

/// Keep `cand` if it is earlier on the clock than what `slot` already holds.
///
/// Ties keep the incumbent, so the first date in the year wins. Around a
/// solstice consecutive days differ by only a second or two, but never by
/// exactly zero at the resolution the bisection works to.
fn keep_earlier(slot: &mut Option<ExtremeEvent>, cand: ExtremeEvent) {
    let better = match slot {
        Some(current) => cand.clock_seconds < current.clock_seconds,
        None => true,
    };
    if better {
        *slot = Some(cand);
    }
}

/// Keep `cand` if it is later on the clock than what `slot` already holds.
///
/// Ties keep the incumbent, matching [`keep_earlier`].
fn keep_later(slot: &mut Option<ExtremeEvent>, cand: ExtremeEvent) {
    let better = match slot {
        Some(current) => cand.clock_seconds > current.clock_seconds,
        None => true,
    };
    if better {
        *slot = Some(cand);
    }
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
        let (sr0, _) = calc0.solve_from_noon(noon0).unwrap();
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

        let (sr_ds, _) = calc_ds.solve_from_noon(noon0).unwrap();
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
        let (tw_low, _) = calc_low.solve_from_noon(noon).unwrap();
        let tw_low = tw_low.unwrap().0;

        // High altitude observer
        let calc_high =
            SolarCalc { lat, lon, alt: 2000.0, delta_t, refr: refraction, target: target_alt };
        let (tw_high, _) = calc_high.solve_from_noon(noon).unwrap();
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

                let (sr, ss) = calc.solve_from_noon(noon).expect("SPA failed at equator");

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

            let (sr, ss) = match calc.solve_from_noon(noon) {
                Ok(pair) => pair,
                Err(_) => continue, // SPA out-of-range for extreme years
            };

            let sr = sr.unwrap().0;
            let ss = ss.unwrap().0;

            let day_len = (ss - sr).num_seconds();

            assert!(day_len > 0);
            assert!(day_len < 24 * 3600);
        }
    }

    /// Build a sea-level NOAA-style calculator for the yearly extreme tests.
    fn extremes_calc(lat: f64, lon: f64, year: i32) -> SolarCalc {
        let delta_t = DeltaT::estimate_from_date(year, 6).unwrap();
        SolarCalc {
            lat,
            lon,
            alt: 0.0,
            delta_t,
            refr: Some(RefractionCorrection::standard()),
            target: -SOLAR_RADIUS_DEG - horizon_dip_deg(lat, 0.0),
        }
    }

    #[test]
    fn test_yearly_extremes_straddle_the_solstices() {
        let lat = 60.4;
        let lon = 25.1;
        let ex = extremes_calc(lat, lon, 2025).yearly_extremes(Helsinki, 2025).unwrap();

        // Southern Finland: every day of the year has both events.
        assert_eq!(ex.days_scanned, 365);
        assert_eq!(ex.days_with_sunrise, 365);
        assert_eq!(ex.days_with_sunset, 365);

        let er = ex.earliest_sunrise.unwrap().time;
        let ls = ex.latest_sunset.unwrap().time;
        let lr = ex.latest_sunrise.unwrap().time;
        let es = ex.earliest_sunset.unwrap().time;

        // The equation of time makes the four extremes straddle the solstices
        // instead of landing on them: the earliest sunrise comes days before
        // the June solstice and the latest sunset days after it, while the
        // earliest sunset precedes the December solstice and the latest
        // sunrise follows it.
        assert_eq!(er.month(), 6, "earliest sunrise should be in June, got {}", er);
        assert!(er.day() < 21, "earliest sunrise should precede the solstice, got {}", er);

        assert_eq!(ls.month(), 6, "latest sunset should be in June, got {}", ls);
        assert!(ls.day() > 21, "latest sunset should follow the solstice, got {}", ls);

        assert_eq!(es.month(), 12, "earliest sunset should be in December, got {}", es);
        assert!(es.day() < 21, "earliest sunset should precede the solstice, got {}", es);

        assert_eq!(lr.month(), 12, "latest sunrise should be in December, got {}", lr);
        assert!(lr.day() > 21, "latest sunrise should follow the solstice, got {}", lr);

        // Sanity: the extremes are ordered on the clock, and at this latitude
        // no event crosses midnight, so the ranking key is exactly the
        // wall-clock reading.
        let (er_key, lr_key) =
            (ex.earliest_sunrise.unwrap().clock_seconds, ex.latest_sunrise.unwrap().clock_seconds);
        let (es_key, ls_key) =
            (ex.earliest_sunset.unwrap().clock_seconds, ex.latest_sunset.unwrap().clock_seconds);

        assert!(er_key < lr_key);
        assert!(es_key < ls_key);
        assert_eq!(er_key, i64::from(er.num_seconds_from_midnight()));
        assert_eq!(ls_key, i64::from(ls.num_seconds_from_midnight()));
    }

    #[test]
    fn test_yearly_extremes_southern_hemisphere() {
        use chrono_tz::Australia::Sydney;

        let ex = extremes_calc(-33.87, 151.21, 2025).yearly_extremes(Sydney, 2025).unwrap();

        let er = ex.earliest_sunrise.unwrap().time;
        let ls = ex.latest_sunset.unwrap().time;
        let lr = ex.latest_sunrise.unwrap().time;
        let es = ex.earliest_sunset.unwrap().time;

        // Sunsets follow the seasons, flipped: the summer maximum spills past
        // New Year — the scan is a calendar year, so it lands in early January
        // rather than December — and the winter minimum sits in mid-June. Both
        // are broad plateaus (Jan 6-9 all read 20:09), so only the month is
        // worth asserting; which day inside the plateau wins is decided by
        // seconds and is model-sensitive.
        assert!(matches!(ls.month(), 12 | 1), "latest sunset should be in Dec/Jan, got {}", ls);
        assert_eq!(es.month(), 6, "earliest sunset should be in June, got {}", es);

        // Sunrises do not follow the seasons here. Sydney's DST runs October
        // to April and both sunrise extremes land on its edges, not near a
        // solstice:
        //
        //   Oct 4, 05:28 AEST - last morning before the clocks go forward.
        //                       From Oct 5 the same sunrise reads an hour
        //                       later, and December only recovers to 05:37.
        //   Apr 5, 07:10 AEDT - last morning before they go back. The June
        //                       maximum only reaches 07:01 AEST.
        //
        // Each wins by about nine minutes over the seasonal candidate, and by
        // a further minute over the day beside it, so the dates are stable.
        // Ranking by clock time instead of by UTC instant is the whole point
        // of the feature: for a DST zone these are the right answers, not
        // artefacts to be assert-ed away.
        assert_eq!(
            (er.month(), er.day()),
            (10, 4),
            "earliest sunrise should be the last morning before DST starts, got {}",
            er
        );
        assert_eq!(
            (lr.month(), lr.day()),
            (4, 5),
            "latest sunrise should be the last morning before DST ends, got {}",
            lr
        );
    }

    #[test]
    fn test_yearly_extremes_polar_gaps_are_skipped() {
        use chrono_tz::Europe::Oslo;

        // Tromsø: polar night around midwinter, midnight sun around midsummer.
        let ex = extremes_calc(69.6492, 18.9553, 2025).yearly_extremes(Oslo, 2025).unwrap();

        assert_eq!(ex.days_scanned, 365);
        assert!(
            ex.days_with_sunrise < 365 && ex.days_with_sunset < 365,
            "expected days without sunrise/sunset at 69.6°N, got {} / {}",
            ex.days_with_sunrise,
            ex.days_with_sunset
        );
        // Well over half the year still has ordinary sunrises and sunsets.
        assert!(ex.days_with_sunrise > 200 && ex.days_with_sunset > 200);

        // All four extremes still exist — the polar gaps are skipped, not
        // treated as missing data or as errors.
        assert!(ex.earliest_sunrise.is_some());
        assert!(ex.latest_sunrise.is_some());
        assert!(ex.earliest_sunset.is_some());
        assert!(ex.latest_sunset.is_some());
    }

    #[test]
    fn test_clock_seconds_carries_across_midnight() {
        use chrono::NaiveDate;

        let day = NaiveDate::from_ymd_opt(2025, 6, 21).unwrap();

        // Same day: the key is the plain wall-clock reading.
        let noon = Helsinki.with_ymd_and_hms(2025, 6, 21, 12, 0, 0).unwrap();
        assert_eq!(clock_seconds(noon, day), 12 * 3600);

        // Sunset spilling into the next morning stays "late", not "early".
        let after_midnight = Helsinki.with_ymd_and_hms(2025, 6, 22, 0, 30, 0).unwrap();
        assert_eq!(clock_seconds(after_midnight, day), 24 * 3600 + 30 * 60);

        // Sunrise on the previous evening stays "early", not "late".
        let before_midnight = Helsinki.with_ymd_and_hms(2025, 6, 20, 23, 50, 0).unwrap();
        assert_eq!(clock_seconds(before_midnight, day), -10 * 60);
    }
}
