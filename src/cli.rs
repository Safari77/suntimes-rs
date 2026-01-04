//! Command-Line Interface Module
//!
//! Handles argument parsing and validation for the suntimes-rs application.

use clap::Parser;
use serde::Deserialize;

// ===================== CLI =====================

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Args {
    #[arg(required_unless_present = "show_build_info")]

    /// Observer latitude in decimal degrees (-90 to 90)
    #[arg(long, allow_hyphen_values = true, value_parser = parse_latitude, env = "ARGOS_SUNTIMES_LATITUDE")]
    pub latitude: f64,
    /// Observer longitude in decimal degrees (-180 to 180)
    #[arg(long, allow_hyphen_values = true, value_parser = parse_longitude, env = "ARGOS_SUNTIMES_LONGITUDE")]
    pub longitude: f64,
    /// Time zone to use ("system", "location", or IANA time zone name)
    #[arg(long, default_value = "system", env = "ARGOS_SUNTIMES_TIMEZONE")]
    pub timezone: String,
    /// Observer altitude above mean sea level (meters, may be negative)
    /// Valid range: -500m (Dead Sea) to 11000m (Troposphere limit for ISA formula)
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true, value_parser = parse_altitude, env = "ARGOS_SUNTIMES_ALTITUDE")]
    pub altitude: f64,

    /// Sun position model to use
    #[arg(long, default_value = "noaa", value_parser = ["noaa", "horizons", "physical"], env = "ARGOS_SUNTIMES_MODEL")]
    pub model: String,
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

    /// Show build info from Cargo.lock at time of building
    #[arg(long)]
    pub show_build_info: bool,

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

fn parse_latitude(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-90.0..=90.0).contains(&v) {
        return Err(format!("Latitude must be between -90 and 90, got {}", v));
    }
    Ok(v)
}

fn parse_longitude(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-180.0..=180.0).contains(&v) {
        return Err(format!("Longitude must be between -180 and 180, got {}", v));
    }
    Ok(v)
}

fn parse_altitude(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("Invalid number: {}", s))?;
    if !(-500.0..=11000.0).contains(&v) {
        return Err(format!("Altitude must be between -500 and 11000 meters, got {}", v));
    }
    Ok(v)
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
