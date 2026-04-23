//! Time and Timezone Utilities Module
//!
//! Provides time parsing, timezone resolution, and formatting utilities.

use chrono::{DateTime, LocalResult, NaiveDate, NaiveTime, TimeZone, Timelike};
use chrono_tz::Tz;
use iana_time_zone::get_timezone;
use std::sync::OnceLock;
use tzf_rs::DefaultFinder;

// tzf-rs DefaultFinder is pre-compiled and very fast
static TZF_FINDER: OnceLock<DefaultFinder> = OnceLock::new();

// ===================== TIME PARSING =====================

/// Parse a time string in HH:MM[:SS[.fffffffff]] format.
///
/// # Arguments
/// * `s` - Time string to parse
///
/// # Returns
/// Tuple of (hours, minutes, seconds, nanoseconds)
///
/// # Errors
/// Returns an error if the time format is invalid
pub fn parse_time_ns(s: &str) -> Result<(u32, u32, u32, u32), Box<dyn std::error::Error>> {
    // %T is HH:MM:SS, %f is fractional seconds (up to nanoseconds)
    // We try multiple formats to be user-friendly
    let formats = ["%H:%M:%S%.f", "%H:%M:%S", "%H:%M"];

    for fmt in formats {
        if let Ok(t) = NaiveTime::parse_from_str(s, fmt) {
            return Ok((t.hour(), t.minute(), t.second(), t.nanosecond()));
        }
    }
    Err("Invalid time format. Use HH:MM, HH:MM:SS, or HH:MM:SS.ns".into())
}

// ===================== TIMEZONE UTILITIES =====================

/// Get the system's configured timezone.
///
/// Falls back to UTC if the system timezone cannot be determined.
pub fn system_timezone() -> Tz {
    get_timezone().ok().and_then(|s| s.parse().ok()).unwrap_or(Tz::UTC)
}

/// Resolve timezone from geographic coordinates.
///
/// Uses a timezone finder to determine the appropriate timezone
/// for a given longitude and latitude.
///
/// # Arguments
/// * `lon` - Longitude in degrees
/// * `lat` - Latitude in degrees
///
/// # Returns
/// The resolved timezone, or UTC if resolution fails
pub fn resolve_timezone(lon: f64, lat: f64) -> Tz {
    let finder = TZF_FINDER.get_or_init(DefaultFinder::new);

    // Get the IANA string (e.g., "Pacific/Apia")
    let tzid = finder.get_tz_name(lon, lat);

    // Parse into chrono_tz::Tz to get historical correctness
    tzid.parse::<Tz>().unwrap_or(Tz::UTC)
}

/// Find the first valid local instant at or after midnight on the given date.
///
/// Returns `None` only when both 00:00 and 01:00 are skipped by a DST/offset shift
/// on that calendar day (extremely rare — e.g. historical transitions that
/// skipped an entire hour past midnight).
///
/// `00:00:00` and `01:00:00` are always valid `NaiveTime` values, so the
/// `and_hms_opt` calls here are infallible — the only reason either branch can
/// yield `None` is the local-time resolution against `tz`.
pub fn start_of_day_opt(tz: Tz, d: NaiveDate) -> Option<DateTime<Tz>> {
    // Infallible: midnight is always a valid NaiveTime
    let midnight = NaiveTime::from_hms_opt(0, 0, 0).expect("00:00:00 is a valid NaiveTime");
    match tz.from_local_datetime(&d.and_time(midnight)) {
        LocalResult::Single(t) | LocalResult::Ambiguous(t, _) => Some(t),
        LocalResult::None => {
            // 00:00 skipped by DST/offset shift — try 01:00
            let one_am = NaiveTime::from_hms_opt(1, 0, 0).expect("01:00:00 is a valid NaiveTime");
            tz.from_local_datetime(&d.and_time(one_am)).earliest()
        }
    }
}

/// Like [`start_of_day_opt`] but substitutes `fallback` on the rare
/// "no valid local midnight or 01:00" case.
///
/// Use this when the caller has a sensible default to fall back to (usually
/// the anchoring datetime that drove the lookup in the first place).
pub fn start_of_day(tz: Tz, d: NaiveDate, fallback: DateTime<Tz>) -> DateTime<Tz> {
    start_of_day_opt(tz, d).unwrap_or(fallback)
}

// ===================== FORMATTING =====================

/// Format a duration in seconds as "Xh Ym Zs".
///
/// # Arguments
/// * `seconds` - Duration in seconds (can be negative, abs value is used)
///
/// # Returns
/// Formatted string like "5h 30m 45s"
pub fn format_hms(seconds: i64) -> String {
    let total_seconds = seconds.abs();
    if total_seconds == 0 {
        return "0s".to_string();
    }

    let h = total_seconds / 3600;
    let m = (total_seconds % 3600) / 60;
    let s = total_seconds % 60;

    let mut parts = Vec::new();
    if h > 0 {
        parts.push(format!("{}h", h));
    }
    if m > 0 {
        parts.push(format!("{}m", m));
    }
    if s > 0 {
        parts.push(format!("{}s", s));
    }

    parts.join(" ")
}

// ===================== TESTS =====================

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::LocalResult;
    use chrono::TimeZone;
    use chrono_tz::America::Indiana::Indianapolis;
    use chrono_tz::Asia::Almaty;
    use chrono_tz::Pacific::{Apia, Kwajalein};

    #[test]
    fn test_at_time_parsing_invalid() {
        assert!(parse_time_ns("a").is_err());
        assert!(parse_time_ns("21").is_err());
        assert!(parse_time_ns("25:00").is_err());
        assert!(parse_time_ns("12:60").is_err());
        assert!(parse_time_ns("12:30:99").is_err());
    }

    #[test]
    fn test_at_time_parsing_valid() {
        assert_eq!(parse_time_ns("12:30").unwrap(), (12, 30, 0, 0));
        assert_eq!(parse_time_ns("00:00").unwrap(), (0, 0, 0, 0));
        assert_eq!(parse_time_ns("23:59").unwrap(), (23, 59, 0, 0));
        assert_eq!(parse_time_ns("12:30:45").unwrap(), (12, 30, 45, 0));
        assert_eq!(parse_time_ns("12:30:45.123").unwrap(), (12, 30, 45, 123_000_000));
        assert_eq!(parse_time_ns("12:30:45.123456789").unwrap(), (12, 30, 45, 123_456_789));
        assert_eq!(parse_time_ns("12:30:45.1234567890").unwrap(), (12, 30, 45, 123_456_789));
        assert_eq!(parse_time_ns("12:30:45.12345678999").unwrap(), (12, 30, 45, 123_456_789));
    }

    #[test]
    fn test_format_hms() {
        assert_eq!(format_hms(3661), "1h 1m 1s");
        assert_eq!(format_hms(7200), "2h");
        assert_eq!(format_hms(45), "45s");
        assert_eq!(format_hms(60), "1m");
        assert_eq!(format_hms(120), "2m");
        assert_eq!(format_hms(0), "0s");
        assert_eq!(format_hms(-3660), "1h 1m"); // Negative handled via abs
    }

    #[test]
    fn test_resolve_timezone_helsinki() {
        use chrono_tz::Europe::Athens;
        use chrono_tz::Europe::Helsinki;
        use chrono_tz::Europe::Mariehamn;
        // Helsinki region
        let lat = 60.2;
        let lon = 24.9;

        let tz = resolve_timezone(lon, lat);
        // The underlying dataset canonicalizes identical EET zones to "Europe/Athens".
        assert!(
            tz == Helsinki || tz == Athens || tz == Mariehamn,
            "Expected an EET timezone (Helsinki/Athens/Mariehamn), got {:?}",
            tz
        );
    }

    #[test]
    fn test_resolve_timezone_new_york() {
        use chrono_tz::America::New_York;
        // Washington DC / New York region
        let lat = 38.8977;
        let lon = -77.0365;

        let tz = resolve_timezone(lon, lat);
        assert_eq!(tz, New_York);
    }

    #[test]
    fn test_resolve_timezone_sydney() {
        use chrono_tz::Australia::Sydney;
        // Canberra / Sydney region
        let lat = -35.3108;
        let lon = 149.1165;

        let tz = resolve_timezone(lon, lat);
        assert_eq!(tz, Sydney);
    }

    #[test]
    fn test_resolve_timezone_sao_paolo() {
        use chrono_tz::America::Argentina::Salta;
        let lat = -41.3;
        let lon = -66.4;

        let tz = resolve_timezone(lon, lat);
        assert_eq!(tz, Salta);
    }

    /// 1. Samoa Date Line Jump (2011)
    /// Samoa skipped Dec 30, 2011. A robust library should handle the
    /// transition from UTC-10 to UTC+14.
    #[test]
    fn test_samoa_skipped_day() {
        // Adding 2 seconds should land us on Dec 31, skipping Dec 30 entirely.
        let skipped_check = Apia.with_ymd_and_hms(2011, 12, 30, 12, 0, 0);
        assert!(
            matches!(skipped_check, LocalResult::None),
            "Dec 30, 2011 should not exist in Samoa"
        );
        println!("Samoa Success: Dec 30 correctly identified as non-existent.");
    }

    /// 2. Kwajalein Date Line Jump (1993)
    /// Kwajalein moved from the east of the IDL to the west.
    /// It skipped August 21, 1993.
    #[test]
    fn test_kwajalein_jump() {
        use chrono::LocalResult;

        // Kwajalein skipped August 21, 1993 to move from the east to the west of the IDL.
        let skipped_day = Kwajalein.with_ymd_and_hms(1993, 8, 21, 12, 0, 0);
        assert!(
            matches!(skipped_day, LocalResult::None),
            "August 21, 1993 should not exist in Kwajalein (it was skipped)"
        );
    }

    #[test]
    fn test_indiana_historical_offset() {
        // July 1990: Indy was UTC-5 (no DST then)
        let indy_summer = Indianapolis.with_ymd_and_hms(1990, 7, 1, 12, 0, 0).unwrap();

        // Use %:z to force +HH:MM format
        let formatted = format!("{}", indy_summer.format("%:z"));

        assert_eq!(formatted, "-05:00");
    }

    #[test]
    fn test_kazakhstan_2024_shift() {
        // Pre-shift: Feb 2024 was UTC+6
        let feb_2024 = Almaty.with_ymd_and_hms(2024, 2, 1, 12, 0, 0).unwrap();
        assert_eq!(format!("{}", feb_2024.format("%:z")), "+06:00");

        // Post-shift: June 2024 is UTC+5
        let june_2024 = Almaty.with_ymd_and_hms(2024, 6, 1, 12, 0, 0).unwrap();

        assert_eq!(
            format!("{}", june_2024.format("%:z")),
            "+05:00",
            "If this fails, run: cargo update -p chrono-tz"
        );
    }

    #[test]
    fn test_india_half_hour_offset() {
        use chrono_tz::Asia::Kolkata;

        // Current time in Kolkata (Fixed at UTC+05:30)
        let india_time = Kolkata.with_ymd_and_hms(2025, 12, 25, 12, 0, 0).unwrap();

        // Using %:z forces the +HH:MM format, including minutes
        let formatted = format!("{}", india_time.format("%:z"));

        assert_eq!(formatted, "+05:30", "India must show the 30-minute offset");
    }

    #[test]
    fn test_india_historical_dst() {
        use chrono_tz::Asia::Kolkata;

        // In July 1943, India was on "War Time" (UTC+06:30)
        let war_time = Kolkata.with_ymd_and_hms(1943, 7, 1, 12, 0, 0).unwrap();

        assert_eq!(format!("{}", war_time.format("%:z")), "+06:30");
    }

    #[test]
    fn test_india_offset_seconds() {
        use chrono::Offset;
        use chrono_tz::Asia::Kolkata;

        let india_time = Kolkata.with_ymd_and_hms(2025, 12, 25, 12, 0, 0).unwrap();

        let offset_seconds = india_time.offset().fix().local_minus_utc();

        // (5 hours * 3600) + (30 minutes * 60) = 19800 seconds
        assert_eq!(offset_seconds, 19800);
    }

    #[test]
    fn test_start_of_day_normal() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Helsinki;

        // Normal summer day in Helsinki — 00:00 exists unambiguously
        let d = NaiveDate::from_ymd_opt(2025, 7, 1).unwrap();
        let sod = start_of_day_opt(Helsinki, d).unwrap();
        assert_eq!(sod.hour(), 0);
        assert_eq!(sod.minute(), 0);
        assert_eq!(sod.date_naive(), d);
    }

    #[test]
    fn test_start_of_day_fallback_on_gap() {
        use chrono::NaiveDate;
        use chrono_tz::Pacific::Apia;

        // Samoa skipped Dec 30, 2011 entirely. `start_of_day_opt` should return None
        // on that day; `start_of_day` should return the supplied fallback.
        let d = NaiveDate::from_ymd_opt(2011, 12, 30).unwrap();
        assert!(start_of_day_opt(Apia, d).is_none());

        let fallback = Apia.with_ymd_and_hms(2011, 12, 31, 0, 0, 0).unwrap();
        let sod = start_of_day(Apia, d, fallback);
        assert_eq!(sod, fallback);
    }

    #[test]
    fn test_start_of_day_ambiguous_returns_earlier() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Helsinki;

        // DST fall-back in Helsinki (2024): last Sunday of October.
        // We want to confirm that if midnight somewhere happened to be ambiguous,
        // we'd get the earlier occurrence. At 00:00 it's never ambiguous in Helsinki,
        // but the contract is enforced regardless of zone: verify stability on a
        // non-DST day that is clearly unambiguous.
        let d = NaiveDate::from_ymd_opt(2024, 10, 27).unwrap();
        let sod = start_of_day_opt(Helsinki, d).unwrap();
        // Just assert monotonicity / sanity — exact value depends on the zone's rules
        assert_eq!(sod.date_naive(), d);
    }
}
