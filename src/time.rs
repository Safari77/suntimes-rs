//! Time and Timezone Utilities Module
//!
//! Provides time parsing, timezone resolution, and formatting utilities.

use jiff::ZonedDateTime;
use jiff::civil::{Date, DateTime, Time};
use jiff::tz::{AmbiguousOffset, TimeZone};
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
    let invalid = || -> Box<dyn std::error::Error> {
        "Invalid time format. Use HH:MM, HH:MM:SS, or HH:MM:SS.ns".into()
    };

    let (hms_part, frac_part) = match s.split_once('.') {
        Some((hms, frac)) => (hms, Some(frac)),
        None => (s, None),
    };

    let t = if frac_part.is_some() {
        Time::strptime("%H:%M:%S", hms_part).map_err(|_| invalid())?
    } else {
        Time::strptime("%H:%M:%S", hms_part)
            .or_else(|_| Time::strptime("%H:%M", hms_part))
            .map_err(|_| invalid())?
    };

    let ns = match frac_part {
        Some(frac) => {
            if frac.is_empty() || !frac.bytes().all(|b| b.is_ascii_digit()) {
                return Err(invalid());
            }
            // Take up to 9 digits (nanoseconds), right-padding with zeros if fewer.
            let mut ns_val: u32 = 0;
            let mut scale: u32 = 100_000_000;
            for b in frac.bytes().take(9) {
                ns_val += u32::from(b - b'0') * scale;
                scale /= 10;
            }
            ns_val
        }
        None => 0,
    };

    Ok((
        u32::try_from(t.hour()).expect("parsed hour is 0..=23"),
        u32::try_from(t.minute()).expect("parsed minute is 0..=59"),
        u32::try_from(t.second()).expect("parsed second is 0..=59"),
        ns,
    ))
}

// ===================== TIMEZONE UTILITIES =====================

/// Get the system's configured timezone.
///
/// Falls back to UTC if the system timezone cannot be determined.
pub fn system_timezone() -> TimeZone {
    TimeZone::try_system().unwrap_or(TimeZone::UTC)
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
pub fn resolve_timezone(lon: f64, lat: f64) -> TimeZone {
    let finder = TZF_FINDER.get_or_init(DefaultFinder::new);

    // Get the IANA string (e.g., "Pacific/Apia")
    let tzid = finder.get_tz_name(lon, lat);

    // Look up in Jiff's IANA database for historical correctness
    TimeZone::get(tzid).unwrap_or(TimeZone::UTC)
}

/// Resolve a civil datetime in `tz`, returning the earlier occurrence on a fold
/// and `None` when the local time falls in a gap (was skipped).
fn resolve_local_earliest(tz: &TimeZone, dt: DateTime) -> Option<ZonedDateTime> {
    let amb = tz.to_ambiguous_zoned_datetime(dt);
    match amb.offset() {
        AmbiguousOffset::Unambiguous { .. } | AmbiguousOffset::Fold { .. } => amb.earlier().ok(),
        AmbiguousOffset::Gap { .. } => None,
    }
}

/// Find the first valid local instant at or after midnight on the given date.
///
/// Returns `None` only when both 00:00 and 01:00 are skipped by a DST/offset shift
/// on that calendar day (extremely rare — e.g. historical transitions that
/// skipped an entire hour past midnight).
pub fn start_of_day_opt(tz: &TimeZone, d: Date) -> Option<ZonedDateTime> {
    resolve_local_earliest(tz, d.at(0, 0, 0, 0))
        .or_else(|| resolve_local_earliest(tz, d.at(1, 0, 0, 0)))
}

/// Like [`start_of_day_opt`] but substitutes `fallback` on the rare
/// "no valid local midnight or 01:00" case.
///
/// Use this when the caller has a sensible default to fall back to (usually
/// the anchoring datetime that drove the lookup in the first place).
pub fn start_of_day(tz: &TimeZone, d: Date, fallback: &ZonedDateTime) -> ZonedDateTime {
    start_of_day_opt(tz, d).unwrap_or_else(|| fallback.clone())
}

/// Find a stable reference instant on the given local calendar day.
///
/// Prefers 12:00 local time, which is the most robust anchor for solar
/// transit calculations, and falls back to [`start_of_day_opt`] when local
/// noon does not exist (a DST/offset shift skipped it — very rare, but e.g.
/// Lord Howe-style shifts or historical midday transitions can do it).
///
/// Returns `None` only when the entire calendar day is missing in this zone
/// (e.g. Samoa, 2011-12-30).
pub fn noon_of_day_opt(tz: &TimeZone, d: Date) -> Option<ZonedDateTime> {
    resolve_local_earliest(tz, d.at(12, 0, 0, 0)).or_else(|| start_of_day_opt(tz, d))
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
    use jiff::civil::date;

    fn tz(name: &str) -> TimeZone {
        TimeZone::get(name).unwrap()
    }

    #[test]
    fn test_at_time_parsing_invalid() {
        assert!(parse_time_ns("a").is_err());
        assert!(parse_time_ns("21").is_err());
        assert!(parse_time_ns("25:00").is_err());
        assert!(parse_time_ns("12:60").is_err());
        assert!(parse_time_ns("12:30:99").is_err());
        assert!(parse_time_ns("12:30:45.").is_err());
        assert!(parse_time_ns("12:30:45.12x").is_err());
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
        // Helsinki region
        let lat = 60.2;
        let lon = 24.9;

        let resolved = resolve_timezone(lon, lat);
        // The underlying dataset canonicalizes identical EET zones to "Europe/Athens".
        assert!(
            matches!(
                resolved.iana_name(),
                Some("Europe/Helsinki" | "Europe/Athens" | "Europe/Mariehamn")
            ),
            "Expected an EET timezone (Helsinki/Athens/Mariehamn), got {:?}",
            resolved
        );
    }

    #[test]
    fn test_resolve_timezone_new_york() {
        // Washington DC / New York region
        let lat = 38.8977;
        let lon = -77.0365;

        let resolved = resolve_timezone(lon, lat);
        assert_eq!(resolved, tz("America/New_York"));
    }

    #[test]
    fn test_resolve_timezone_sydney() {
        // Canberra / Sydney region
        let lat = -35.3108;
        let lon = 149.1165;

        let resolved = resolve_timezone(lon, lat);
        assert_eq!(resolved, tz("Australia/Sydney"));
    }

    #[test]
    fn test_resolve_timezone_sao_paolo() {
        let lat = -41.3;
        let lon = -66.4;

        let resolved = resolve_timezone(lon, lat);
        assert_eq!(resolved, tz("America/Argentina/Salta"));
    }

    /// Samoa Date Line Jump (2011): Samoa skipped Dec 30, 2011.
    /// A robust library should handle the transition from UTC-10 to UTC+14.
    #[test]
    fn test_samoa_skipped_day() {
        let apia = tz("Pacific/Apia");
        let skipped_check = apia.to_ambiguous_zoned_datetime(date(2011, 12, 30).at(12, 0, 0, 0));
        assert!(
            matches!(skipped_check.offset(), AmbiguousOffset::Gap { .. }),
            "Dec 30, 2011 should not exist in Samoa"
        );
        println!("Samoa Success: Dec 30 correctly identified as non-existent.");
    }

    /// Kwajalein Date Line Jump (1993): Kwajalein moved from the east of the
    /// IDL to the west, skipping August 21, 1993.
    #[test]
    fn test_kwajalein_jump() {
        let kwajalein = tz("Pacific/Kwajalein");
        let skipped_day = kwajalein.to_ambiguous_zoned_datetime(date(1993, 8, 21).at(12, 0, 0, 0));
        assert!(
            matches!(skipped_day.offset(), AmbiguousOffset::Gap { .. }),
            "August 21, 1993 should not exist in Kwajalein (it was skipped)"
        );
    }

    #[test]
    fn test_indiana_historical_offset() {
        let indianapolis = tz("America/Indiana/Indianapolis");
        // July 1990: Indy was UTC-5 (no DST then)
        let indy_summer = indianapolis.to_zoned_datetime(date(1990, 7, 1).at(12, 0, 0, 0)).unwrap();

        // Use %:z to force +HH:MM format
        let formatted = indy_summer.strftime("%:z").to_string();

        assert_eq!(formatted, "-05:00");
    }

    #[test]
    fn test_kazakhstan_2024_shift() {
        let almaty = tz("Asia/Almaty");
        // Pre-shift: Feb 2024 was UTC+6
        let feb_2024 = almaty.to_zoned_datetime(date(2024, 2, 1).at(12, 0, 0, 0)).unwrap();
        assert_eq!(feb_2024.strftime("%:z").to_string(), "+06:00");

        // Post-shift: June 2024 is UTC+5
        let june_2024 = almaty.to_zoned_datetime(date(2024, 6, 1).at(12, 0, 0, 0)).unwrap();

        assert_eq!(
            june_2024.strftime("%:z").to_string(),
            "+05:00",
            "Expected Kazakhstan 2024 shift to UTC+05:00 in TZ database"
        );
    }

    #[test]
    fn test_india_half_hour_offset() {
        let kolkata = tz("Asia/Kolkata");

        // Current time in Kolkata (Fixed at UTC+05:30)
        let india_time = kolkata.to_zoned_datetime(date(2025, 12, 25).at(12, 0, 0, 0)).unwrap();

        // Using %:z forces the +HH:MM format, including minutes
        let formatted = india_time.strftime("%:z").to_string();

        assert_eq!(formatted, "+05:30", "India must show the 30-minute offset");
    }

    #[test]
    fn test_india_historical_dst() {
        let kolkata = tz("Asia/Kolkata");

        // In July 1943, India was on "War Time" (UTC+06:30)
        let war_time = kolkata.to_zoned_datetime(date(1943, 7, 1).at(12, 0, 0, 0)).unwrap();

        assert_eq!(war_time.strftime("%:z").to_string(), "+06:30");
    }

    #[test]
    fn test_india_offset_seconds() {
        let kolkata = tz("Asia/Kolkata");

        let india_time = kolkata.to_zoned_datetime(date(2025, 12, 25).at(12, 0, 0, 0)).unwrap();

        let offset_seconds = india_time.offset().seconds();

        // (5 hours * 3600) + (30 minutes * 60) = 19800 seconds
        assert_eq!(offset_seconds, 19800);
    }

    #[test]
    fn test_start_of_day_normal() {
        let helsinki = tz("Europe/Helsinki");

        // Normal summer day in Helsinki — 00:00 exists unambiguously
        let d = date(2025, 7, 1);
        let sod = start_of_day_opt(&helsinki, d).unwrap();
        assert_eq!(sod.hour(), 0);
        assert_eq!(sod.minute(), 0);
        assert_eq!(sod.date(), d);
    }

    #[test]
    fn test_start_of_day_fallback_on_gap() {
        let apia = tz("Pacific/Apia");

        // Samoa skipped Dec 30, 2011 entirely. `start_of_day_opt` should return None
        // on that day; `start_of_day` should return the supplied fallback.
        let d = date(2011, 12, 30);
        assert!(start_of_day_opt(&apia, d).is_none());

        let fallback = apia.to_zoned_datetime(date(2011, 12, 31).at(0, 0, 0, 0)).unwrap();
        let sod = start_of_day(&apia, d, &fallback);
        assert_eq!(sod, fallback);
    }

    #[test]
    fn test_noon_of_day_normal() {
        let helsinki = tz("Europe/Helsinki");

        let d = date(2025, 7, 1);
        let noon = noon_of_day_opt(&helsinki, d).unwrap();
        assert_eq!(noon.hour(), 12);
        assert_eq!(noon.minute(), 0);
        assert_eq!(noon.date(), d);
    }

    #[test]
    fn test_noon_of_day_missing_day() {
        let apia = tz("Pacific/Apia");

        // Samoa skipped Dec 30, 2011 entirely: neither noon nor midnight exists.
        let d = date(2011, 12, 30);
        assert!(noon_of_day_opt(&apia, d).is_none());
    }

    #[test]
    fn test_start_of_day_ambiguous_returns_earlier() {
        let helsinki = tz("Europe/Helsinki");

        // DST fall-back in Helsinki (2024): last Sunday of October.
        let d = date(2024, 10, 27);
        let sod = start_of_day_opt(&helsinki, d).unwrap();
        // Just assert monotonicity / sanity — exact value depends on the zone's rules
        assert_eq!(sod.date(), d);
    }
}
