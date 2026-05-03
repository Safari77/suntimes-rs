//! Output Formatting Module
//!
//! Provides formatting functions for terminal and Argos (GNOME Shell) output.

use chrono::DateTime;
use chrono_tz::Tz;
use solar_positioning::types::SunriseResult;

use crate::air_quality::{self, AirQualityResponse};
use crate::solar::SunEvent;
use crate::solar_panel::{self, SolarPanelOutput, TrackingMode};
use crate::time::format_hms;

// ===================== ARGOS OUTPUT =====================

/// Print output in Argos (GNOME Shell extension) format.
///
/// This format is designed for the Argos shell extension which displays
/// information in the top panel of GNOME Shell.
///
/// # Arguments
/// * `now_pos` - Current sun position
/// * `sr` - Sunrise event (time and azimuth)
/// * `ss` - Sunset event (time and azimuth)
/// * `transit_res` - Transit result (for polar day/night detection)
/// * `transit` - Solar noon time
/// * `transit_pos` - Sun position at solar noon
/// * `len_today` - Day length in seconds
/// * `len_tomorrow` - Tomorrow's day length in seconds
/// * `solar_panel_output` - Optional solar panel output
/// * `daily_energy` - Optional daily energy estimate
/// * `uv_data` - Optional UV index data (current, max)
/// * `air_quality` - Optional air-quality / pollen response
/// * `show_aqi` - Whether to render the Air Quality section
/// * `show_pollen` - Whether to render the Pollen Warning section
#[allow(clippy::too_many_arguments)]
pub fn print_argos(
    now_pos: &solar_positioning::SolarPosition,
    sr: &Option<SunEvent>,
    ss: &Option<SunEvent>,
    transit_res: &SunriseResult<DateTime<Tz>>,
    transit: DateTime<Tz>,
    transit_pos: &solar_positioning::SolarPosition,
    len_today: Option<i64>,
    len_tomorrow: Option<i64>,
    solar_panel_output: Option<&SolarPanelOutput>,
    daily_energy: Option<f64>,
    uv_data: Option<(f64, f64)>, // (Current UV, Max UV)
    air_quality: Option<&AirQualityResponse>,
    show_aqi: bool,
    show_pollen: bool,
) {
    let elevation = now_pos.elevation_angle();
    let sunchar = if elevation > 3.0 {
        "🟡"
    } else if elevation >= -3.0 {
        "🟠"
    } else if elevation >= -18.0 {
        "🌚"
    } else {
        "🌠"
    };
    // Section 1: Current sun position (menu bar item)
    if let Some(output) = solar_panel_output {
        println!(
            "{} Alt: {:.2}° | ⚡ {}",
            sunchar,
            elevation,
            solar_panel::format_power(output.power_w)
        );
    } else {
        println!("{} Alt: {:.2}° Az: {:.2}°", sunchar, elevation, now_pos.azimuth());
    }
    println!("---");

    // Section 2: Today's sun events
    if let Some((t, az)) = sr {
        println!("Sunrise: {} (Az: {:.2}°)", t.format("%H:%M:%S"), az);
    }
    println!(
        "Max Altitude: {:.2}° at {}",
        transit_pos.elevation_angle(),
        transit.format("%H:%M:%S")
    );
    if let Some((t, az)) = ss {
        println!("Sunset: {} (Az: {:.2}°)", t.format("%H:%M:%S"), az);
    }
    println!("---");

    // Section 3: Day length comparison
    if let Some(today) = len_today {
        println!("Daylight: {}", format_hms(today));
        if let Some(tomorrow) = len_tomorrow {
            let diff = tomorrow - today;
            if diff.abs() == 0 {
                println!("Tomorrow: same length");
            } else {
                println!(
                    "Tomorrow: {} {}",
                    format_hms(diff),
                    if diff > 0 { "longer" } else { "shorter" }
                );
            }
        }
    } else {
        match transit_res {
            SunriseResult::AllDay { .. } => println!("Polar Day (Midnight Sun)."),
            SunriseResult::AllNight { .. } => println!("Polar Night."),
            _ => println!("Sun does not cross target altitude today."),
        }
    }

    // Section 4: UV Data
    if let Some((curr, max)) = uv_data {
        println!("---");
        println!("🟣 UV Index");
        println!("Current: {:.1}", curr);
        println!("Max today: {:.1}", max);
    }

    // Section 5: Air Quality / Pollen (Argos format)
    if let Some(aq) = air_quality {
        print_argos_air_quality(aq, show_aqi, show_pollen);
    }

    // Section 6: Solar panel output (if configured)
    if let Some(output) = solar_panel_output {
        println!("---");
        println!("⚡ Solar Panel");
        println!("Power: {}", solar_panel::format_power(output.power_w));
        println!("POA: {}", solar_panel::format_irradiance(output.irradiance.poa));
        println!("AOI: {:.1}°", output.irradiance.aoi_deg);
        if let Some(energy) = daily_energy {
            println!("Est. Daily: {}", solar_panel::format_energy(energy));
        }
    }
}

/// Render the Air Quality and Pollen Warning sections in Argos format.
///
/// Argos uses `---` separators between groups; each section we emit is
/// preceded by one. If neither flag has any printable content (e.g. the
/// API was queried but returned only trace pollen and `--aqi` is off),
/// nothing is printed and no orphan separator appears.
fn print_argos_air_quality(aq: &AirQualityResponse, show_aqi: bool, show_pollen: bool) {
    if show_aqi {
        println!("---");
        println!("Air Quality Details:");
        if let Some(eaqi) = aq.current.european_aqi {
            println!("European AQI: {:.0} ({})", eaqi, air_quality::european_aqi_category(eaqi));
        }
        if let Some(v) = aq.current.pm2_5 {
            println!("PM2.5: {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.pm10 {
            println!("PM10: {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.ozone {
            println!("O3: {:.0} μg/m³", v);
        }
        if let Some(v) = aq.current.nitrogen_dioxide {
            println!("NO2: {:.0} μg/m³", v);
        }
        if let Some(v) = aq.current.sulphur_dioxide {
            println!("SO2: {:.1} μg/m³", v);
        }
        // Dust is hidden on normal days — only show it when something
        // unusual is happening (Saharan dust transport, sandstorms, etc.)
        if let Some(v) = aq.current.dust
            && v > air_quality::DUST_DISPLAY_THRESHOLD
        {
            println!("Dust: {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.aerosol_optical_depth {
            println!("AOD: {:.2} ({})", v, air_quality::aerosol_optical_depth_category(v));
        }
    }

    if show_pollen {
        let readings = aq.current.pollen_readings();
        if !readings.is_empty() {
            println!("---");
            println!("Pollen Warning:");
            for r in &readings {
                println!("{}: {} ({:.0} grains/m³)", r.display_name, r.severity.label(), r.value);
            }
        }
    }
}

// ===================== TERMINAL OUTPUT =====================

/// Print detailed terminal output for sun events and conditions.
///
/// # Arguments
/// * `sr` - Sunrise event (time and azimuth)
/// * `ss` - Sunset event (time and azimuth)
/// * `transit_res` - Transit result for polar day/night detection
/// * `transit` - Solar noon time
/// * `transit_pos` - Sun position at solar noon
/// * `len_today` - Day length in seconds
/// * `len_tomorrow` - Tomorrow's day length in seconds
/// * `calc` - Solar calculation context
/// * `date` - Current date
#[allow(clippy::too_many_arguments)]
pub fn print_sun_events(
    sr: &Option<SunEvent>,
    ss: &Option<SunEvent>,
    transit_res: &SunriseResult<DateTime<Tz>>,
    transit: DateTime<Tz>,
    transit_pos: &solar_positioning::SolarPosition,
    len_today: Option<i64>,
    len_tomorrow: Option<i64>,
    calc: &crate::solar::SolarCalc,
    date: DateTime<Tz>,
) {
    match (sr, ss) {
        (Some((t1, a1)), Some((t2, a2))) => {
            println!("Sunrise     : {}", t1.format("%H:%M:%S %Z"));
            println!("  Azimuth   : {:8.3}°", a1);
            println!(
                "Max Altitude: {:8.3}° at {}",
                transit_pos.elevation_angle(),
                transit.format("%H:%M:%S")
            );
            println!("Sunset      : {}", t2.format("%H:%M:%S %Z"));
            println!("  Azimuth   : {:8.3}°", a2);

            if let Some(len_today) = len_today {
                println!("Daylight    : {}", format_hms(len_today));

                if let Some(len_tomorrow) = len_tomorrow {
                    let diff = len_tomorrow - len_today;
                    if diff.abs() == 0 {
                        println!("Tomorrow day is same length");
                    } else {
                        println!(
                            "Tomorrow day is {} {}",
                            format_hms(diff),
                            if diff > 0 { "longer" } else { "shorter" }
                        );
                    }
                }
            }
        }
        _ => {
            match transit_res {
                SunriseResult::AllDay { .. } => println!("Polar Day (Midnight Sun)."),
                SunriseResult::AllNight { .. } => println!("Polar Night."),
                _ => println!("Sun does not cross target altitude today."),
            }

            println!(
                "Max Altitude: {:8.3}° at {}",
                transit_pos.elevation_angle(),
                transit.format("%H:%M:%S")
            );

            // `find_next_event` now returns Result. The "next event" line is
            // informational; a real SPA error would already have surfaced in the
            // main computation above, so we silently ignore it here.
            if let Ok(Some((kind, t))) = calc.find_next_event(date) {
                println!("Next {} on {} at {}", kind, t.date_naive(), t.format("%H:%M:%S %Z"));
            }
        }
    }
}

/// Print UV index information.
///
/// # Arguments
/// * `panel_output_label` - Label for conditions (e.g., "now", "at 12:00:00")
/// * `uv_current` - Current UV index
/// * `uv_max` - Maximum UV index for the day
pub fn print_uv_info(panel_output_label: &str, uv_current: f64, uv_max: f64) {
    println!();
    println!("=== UV Index ===");
    println!("Conditions ({})", panel_output_label);
    println!("  Current       : {:.1}", uv_current);
    println!("  Max Today     : {:.1} (at solar noon)", uv_max);
}

/// Print Air Quality and Pollen sections in terminal format.
///
/// Mirrors the Argos sections, but uses the `=== Title ===` heading
/// style consistent with the rest of the terminal output.
///
/// # Arguments
/// * `aq` - Fetched air-quality response (cached or fresh)
/// * `show_aqi` - Whether to render the Air Quality block
/// * `show_pollen` - Whether to render the Pollen block (skipped when
///   no species are above their trace threshold)
pub fn print_air_quality_terminal(aq: &AirQualityResponse, show_aqi: bool, show_pollen: bool) {
    if show_aqi {
        println!();
        println!("=== Air Quality ===");
        if let Some(eaqi) = aq.current.european_aqi {
            println!("European AQI : {:.0} ({})", eaqi, air_quality::european_aqi_category(eaqi));
        }
        if let Some(v) = aq.current.pm2_5 {
            println!("PM2.5        : {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.pm10 {
            println!("PM10         : {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.ozone {
            println!("O3           : {:.0} μg/m³", v);
        }
        if let Some(v) = aq.current.nitrogen_dioxide {
            println!("NO2          : {:.0} μg/m³", v);
        }
        if let Some(v) = aq.current.sulphur_dioxide {
            println!("SO2          : {:.1} μg/m³", v);
        }
        // Hidden when below background — only show during dust events.
        if let Some(v) = aq.current.dust
            && v > air_quality::DUST_DISPLAY_THRESHOLD
        {
            println!("Dust         : {:.1} μg/m³", v);
        }
        if let Some(v) = aq.current.aerosol_optical_depth {
            println!(
                "AOD          : {:.2} ({})",
                v,
                air_quality::aerosol_optical_depth_category(v)
            );
        }
    }

    if show_pollen {
        let readings = aq.current.pollen_readings();
        if !readings.is_empty() {
            println!();
            println!("=== Pollen ===");
            for r in &readings {
                println!(
                    "{:<8}: {:<10} ({:.0} grains/m³)",
                    r.display_name,
                    r.severity.label(),
                    r.value
                );
            }
        }
    }
}

/// Print solar panel configuration and output.
///
/// # Arguments
/// * `config` - Solar panel configuration
/// * `output` - Current solar panel output
/// * `panel_output_label` - Label for conditions
/// * `daily_energy` - Optional daily energy estimate
pub fn print_solar_panel_info(
    config: &solar_panel::SolarPanelConfig,
    output: &SolarPanelOutput,
    panel_output_label: &str,
    daily_energy: Option<f64>,
) {
    println!();
    println!("=== Solar Panel Output (Ineichen-Perez Clear-Sky Model) ===");

    // Display panel configuration based on tracking mode
    match config.tracking_mode {
        TrackingMode::Fixed => {
            println!(
                "Panel     : {:.2} m² @ {:.0}° tilt, {:.0}° azimuth (fixed)",
                config.area_m2, config.tilt_deg, config.azimuth_deg
            );
        }
        TrackingMode::HorizontalAxis => {
            println!(
                "Panel     : {:.2} m² @ {:.0}° tilt (horizontal single-axis, tracks azimuth)",
                config.area_m2, config.tilt_deg
            );
        }
        TrackingMode::VerticalAxis => {
            println!(
                "Panel     : {:.2} m² @ {:.0}° azimuth (vertical single-axis, tracks altitude)",
                config.area_m2, config.azimuth_deg
            );
        }
        TrackingMode::DualAxis => {
            println!("Panel     : {:.2} m² (dual-axis tracking)", config.area_m2);
        }
    }

    println!(
        "Efficiency: {:.1}% | Turbidity: {:.1} | Albedo: {:.2}",
        config.efficiency * 100.0,
        config.linke_turbidity,
        config.albedo
    );
    println!();

    println!("Conditions ({}):", panel_output_label);
    println!("  Sun elevation : {:8.2}°", output.sun_elevation_deg);
    println!("  Sun azimuth   : {:8.2}°", output.sun_azimuth_deg);
    println!("  Air mass      : {:8.2}", output.air_mass);
    println!("  Angle of inc. : {:8.2}°", output.irradiance.aoi_deg);
    println!();
    println!("Irradiance (clear-sky):");
    println!("  DNI (direct)  : {}", solar_panel::format_irradiance(output.irradiance.dni));
    println!("  DHI (diffuse) : {}", solar_panel::format_irradiance(output.irradiance.dhi));
    println!("  GHI (global)  : {}", solar_panel::format_irradiance(output.irradiance.ghi));
    println!();
    println!("Plane-of-Array (tilted surface):");
    println!("  POA total     : {}", solar_panel::format_irradiance(output.irradiance.poa));
    println!("  POA beam      : {}", solar_panel::format_irradiance(output.irradiance.poa_beam));
    println!(
        "  POA sky diff. : {}",
        solar_panel::format_irradiance(output.irradiance.poa_sky_diffuse)
    );
    println!(
        "  POA ground    : {}",
        solar_panel::format_irradiance(output.irradiance.poa_ground_diffuse)
    );
    println!();
    println!("Power output    : {}", solar_panel::format_power(output.power_w));

    if let Some(energy) = daily_energy {
        println!("Daily estimate  : {} (clear-sky)", solar_panel::format_energy(energy));
    }
}

/// Print sun position at a specific time.
///
/// # Arguments
/// * `dt` - DateTime for the position
/// * `pos` - Sun position
/// * `uv_current` - Current UV index
/// * `uv_max` - Maximum UV index for the day
pub fn print_sun_position_at_time(
    dt: DateTime<Tz>,
    pos: &solar_positioning::SolarPosition,
    uv_current: f64,
    uv_max: f64,
) {
    println!("Sun position at {}:", dt.format("%H:%M:%S%.f %Z"));
    println!("  Azimuth       : {:8.3}°", pos.azimuth());
    println!("  Altitude      : {:8.3}°", pos.elevation_angle());
    println!("  Zenith angle  : {:8.3}°", pos.zenith_angle());
    println!("  UV Index      : {:6.1} (Max today: {:.1})", uv_current, uv_max);
}

/// Print solar panel output at a specific time.
///
/// # Arguments
/// * `dt` - DateTime for the output
/// * `output` - Solar panel output
pub fn print_solar_panel_at_time(dt: DateTime<Tz>, output: &SolarPanelOutput) {
    println!();
    println!("Solar panel at {}:", dt.format("%H:%M:%S"));
    println!("  Angle of inc. : {:8.2}°", output.irradiance.aoi_deg);
    println!("  POA irradiance: {}", solar_panel::format_irradiance(output.irradiance.poa));
    println!("  Power output  : {}", solar_panel::format_power(output.power_w));
}

/// Print yearly optimization results with adjustment schedule.
///
/// # Arguments
/// * `result` - The yearly optimization result
/// * `year` - The year for date formatting
/// * `constraints` - The constraints that were applied
pub fn print_yearly_optimization(
    result: &crate::optimize::YearlyOptimizationResult,
    year: i32,
    constraints: &crate::optimize::OptimizationConstraints,
    is_hsat: bool,
) {
    use chrono::NaiveDate;

    println!();
    if is_hsat {
        println!("=== Yearly HSAT Optimization (Seasonal Tilt Adjustments) ===");
        println!("Year: {}", year);
        println!("Number of tilt adjustment periods: {}", result.periods.len());
        println!();
        println!("HSAT mode: Azimuth tracks sun throughout day, tilt adjusted seasonally");
    } else {
        println!("=== Yearly Solar Panel Optimization (Seasonal Tilt Adjustments) ===");
        println!("Year: {}", year);
        println!("Number of adjustment periods: {}", result.periods.len());
        println!();

        // Explain that azimuth is fixed (real-world behavior)
        if !result.periods.is_empty() {
            let fixed_az = result.periods[0].azimuth_deg;
            let direction = match fixed_az {
                a if (135.0..=225.0).contains(&a) => "south",
                a if (45.0..=135.0).contains(&a) => "east",
                a if (225.0..=315.0).contains(&a) => "west",
                _ => "north",
            };
            println!("Panel orientation: {:.0}° azimuth (facing {})", fixed_az, direction);
            println!("Note: Only tilt is adjusted seasonally (standard practice)");
        }
    }
    println!();

    // Print constraints if non-default
    let default = crate::optimize::OptimizationConstraints::default();
    if (constraints.tilt_min - default.tilt_min).abs() > 0.1
        || (constraints.tilt_max - default.tilt_max).abs() > 0.1
    {
        println!(
            "Tilt range constraint: {:.1}° - {:.1}°",
            constraints.tilt_min, constraints.tilt_max
        );
        println!();
    }

    println!("Adjustment Schedule:");
    println!("{:-<68}", "");
    println!(
        "{:<4} {:>12} {:>12} {:>12} {:>14}",
        "#", "Start Date", "End Date", "Tilt (°)", "Period (kWh)"
    );
    println!("{:-<68}", "");

    for (i, period) in result.periods.iter().enumerate() {
        let start_date = NaiveDate::from_yo_opt(year, period.start_day)
            .map(|d| d.format("%b %d").to_string())
            .unwrap_or_else(|| format!("Day {}", period.start_day));
        let end_date = NaiveDate::from_yo_opt(year, period.end_day)
            .map(|d| d.format("%b %d").to_string())
            .unwrap_or_else(|| format!("Day {}", period.end_day));

        println!(
            "{:<4} {:>12} {:>12} {:>12.1} {:>14}",
            i + 1,
            start_date,
            end_date,
            period.tilt_deg,
            solar_panel::format_energy(period.period_energy_wh),
        );
    }
    println!("{:-<68}", "");
    println!();

    // Summary statistics
    let total_kwh = result.total_energy_wh / 1000.0;
    let fixed_kwh = result.fixed_optimal_energy_wh / 1000.0;
    let improvement = if fixed_kwh > 0.0 { ((total_kwh / fixed_kwh) - 1.0) * 100.0 } else { 0.0 };

    println!("Summary:");
    println!(
        "  Total yearly energy (with {} adjustments): {:.2} kWh",
        result.periods.len(),
        total_kwh
    );
    println!(
        "  {} optimal (no adjustments)           : {:.2} kWh",
        if is_hsat { "HSAT" } else { "Fixed" },
        fixed_kwh
    );

    if result.periods.len() > 1 {
        if improvement > 0.1 {
            println!("  Improvement from seasonal adjustments    : {:+.1}%", improvement);
        } else if improvement < -0.1 {
            println!("  Note: Fixed setting slightly better (rounding): {:+.1}%", improvement);
        } else {
            println!("  Seasonal adjustments provide negligible improvement");
        }
    }

    // Provide practical advice based on number of adjustments
    println!();
    if is_hsat {
        match result.periods.len() {
            1 => println!("Tip: This is HSAT with fixed tilt - no seasonal adjustments."),
            2 => println!(
                "Tip: 2 tilt adjustments/year can optimize HSAT for summer/winter sun angles."
            ),
            4 => println!("Tip: Quarterly tilt adjustments provide good HSAT optimization."),
            n if n <= 12 => println!("Tip: Monthly tilt adjustments maximize HSAT performance."),
            _ => println!("Tip: For frequent adjustments, consider a dual-axis tracker."),
        }
    } else {
        match result.periods.len() {
            1 => println!("Tip: This is a fixed installation - no adjustments needed."),
            2 => println!(
                "Tip: 2 adjustments/year (summer/winter) typically yields ~5% more energy."
            ),
            4 => {
                println!("Tip: 4 adjustments/year (quarterly) typically yields ~7-8% more energy.")
            }
            n if n <= 12 => println!("Tip: Monthly adjustments can yield up to ~10% more energy."),
            _ => println!("Tip: For daily adjustments, consider a single-axis tracker instead."),
        }
    }

    println!();
    println!("(Optimization used {} evaluations)", result.evaluations);
}
