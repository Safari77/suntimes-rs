use chrono::{DateTime, Datelike, Duration, TimeZone, Utc};
use chrono_tz::Tz;
use clap::Parser;
use parse_datetime::parse_datetime_at_date;
use solar_positioning::{
    time::DeltaT,
    types::{RefractionCorrection, SunriseResult},
};

mod cli;
mod geo;
mod optimize;
mod output;
mod solar;
mod solar_panel;
mod time;
mod uv;

use cli::{Args, DepInfo};
use geo::{SOLAR_RADIUS_DEG, horizon_dip_deg};
use solar::{SolarCalc, day_length};
use time::{parse_time_ns, resolve_timezone, system_timezone};

fn parse_relative_to_chrono<Tz>(anchor: &DateTime<Tz>, input: &str) -> Result<DateTime<Tz>, String>
where
    Tz: TimeZone + std::fmt::Display,
    Tz::Offset: std::fmt::Display,
{
    // 1. Convert Chrono Anchor -> Jiff Anchor
    let tz_name = anchor.timezone().to_string();

    let jiff_tz = jiff::tz::TimeZone::get(&tz_name)
        .map_err(|e| format!("Jiff failed to load tz '{}': {}", tz_name, e))?;

    let jiff_anchor =
        jiff::Timestamp::new(anchor.timestamp(), anchor.timestamp_subsec_nanos() as i32)
            .map_err(|e| format!("Timestamp conversion failed: {}", e))?
            .to_zoned(jiff_tz);

    let parsed_jiff = parse_datetime_at_date(jiff_anchor, input)
        .map_err(|e| format!("Failed to parse '{}': {}", input, e))?;

    // Convert Jiff Result -> Chrono
    let ts = parsed_jiff.timestamp();
    let dt_utc = Utc
        .timestamp_opt(ts.as_second(), ts.subsec_nanosecond() as u32)
        .single()
        .ok_or("Resulting timestamp is invalid in Chrono")?;

    Ok(dt_utc.with_timezone(&anchor.timezone()))
}

// ===================== MAIN =====================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().any(|arg| arg == "--show-build-info") {
        println!("Built from Git commit: {}\n", env!("APP_GIT_HASH"));
        const DEP_INFO_RAW: &str = include_str!(env!("DEPS_INFO_PATH"));
        let deps: Vec<DepInfo> = serde_json::from_str(DEP_INFO_RAW).unwrap();

        println!("Found {} dependencies.", deps.len());
        for dep in deps {
            println!("- {} v{}", dep.name, dep.version);
            if let Some(sum) = dep.checksum {
                println!("    Checksum: {}", sum);
            }
            if let Some(src) = dep.source {
                println!("    Source:   {}", src);
            }
        }
        return Ok(());
    }

    let args = Args::parse();
    // Manual override logic:
    // If --no-argos is present, it forces Argos mode off,
    // ignoring ENABLE_ARGOS env var or --argos flag.
    let use_argos = args.argos && !args.no_argos;

    let tz = if args.utc {
        Tz::UTC
    } else {
        match args.timezone.as_str() {
            "system" => system_timezone(),
            "location" => resolve_timezone(args.longitude, args.latitude),
            other => other.parse().unwrap_or(Tz::UTC),
        }
    };
    let now_dt =
        if args.utc { Utc::now().with_timezone(&Tz::UTC) } else { Utc::now().with_timezone(&tz) };

    // Anchor 'today' to the target timezone
    let anchor_time = Utc::now().with_timezone(&tz);
    let date = if let Some(s) = &args.date {
        parse_relative_to_chrono(&anchor_time, s)?
    } else {
        anchor_time
    };

    let base_alt = if args.civil {
        -6.0
    } else if args.nautical {
        -12.0
    } else if args.astro {
        -18.0
    } else {
        -SOLAR_RADIUS_DEG
    };

    // Helper to find the "Logical Midnight" (Start of Day) safely.
    // Returns the first valid time on that day (00:00, or 01:00 if a DST gap
    // skipped midnight, or `date` if neither exists).
    let find_start_of_day =
        |d: chrono::NaiveDate| -> DateTime<Tz> { time::start_of_day(tz, d, date) };

    let delta_t: f64 = DeltaT::estimate_from_date(date.year(), date.month())?;

    // Twilight is geometric; disable refraction if twilight flag is set
    let is_twilight_mode = args.civil || args.nautical || args.astro;
    let calc = match args.model {
        cli::SunModel::Noaa => {
            let dip =
                if is_twilight_mode { 0.0 } else { horizon_dip_deg(args.latitude, args.altitude) };
            // If twilight mode, force refr to None
            let refr = if is_twilight_mode { None } else { Some(RefractionCorrection::standard()) };

            SolarCalc {
                lat: args.latitude,
                lon: args.longitude,
                alt: args.altitude,
                delta_t,
                refr,
                target: base_alt - dip,
            }
        }
        cli::SunModel::Horizons => {
            let refr = if is_twilight_mode { None } else { Some(RefractionCorrection::standard()) };
            SolarCalc {
                lat: args.latitude,
                lon: args.longitude,
                alt: args.altitude,
                delta_t,
                refr,
                target: base_alt,
            }
        }
        cli::SunModel::Physical => {
            let dip =
                if is_twilight_mode { 0.0 } else { horizon_dip_deg(args.latitude, args.altitude) };

            let pressure = if (args.pressure - 1013.25).abs() > f64::EPSILON {
                args.pressure
            } else {
                1013.25 * (1.0 - 2.25577e-5 * args.altitude).powf(5.25588)
            };

            let refr = if is_twilight_mode {
                None
            } else {
                Some(RefractionCorrection::new(pressure, args.temperature)?)
            };

            SolarCalc {
                lat: args.latitude,
                lon: args.longitude,
                alt: args.altitude,
                delta_t,
                refr,
                target: base_alt - dip,
            }
        }
    };

    let transit_res = calc.get_transit(date).ok_or("Failed to compute solar transit")?;
    let transit = calc.extract_transit_time(&transit_res);

    // Get current sun position (only used for Argos or --at "now", but cheap to compute)
    let now_pos = calc.position(now_dt)?;

    // Get sun position at solar noon (max altitude)
    let transit_pos = calc.position(transit)?;

    let (sr, ss) = calc.solve_from_noon(transit)?;

    // Calculate tomorrow's day length
    let len_today = day_length(&sr, &ss);

    let mut len_tomorrow = None;
    let mut check_date = date.date_naive().succ_opt();
    // Look ahead up to 2 days to handle day skips
    for _ in 0..2 {
        if let Some(d) = check_date {
            // Find the first valid local instant on this calendar day.
            // `start_of_day_opt` also correctly handles DST fall-back (`01:00` ambiguous),
            // whereas the previous inline code used `.single()` which discarded that case.
            let tomorrow_dt = time::start_of_day_opt(tz, d);

            if let Some(t_dt) = tomorrow_dt
                && let Some(transit2_res) = calc.get_transit(t_dt)
            {
                let transit2 = calc.extract_transit_time(&transit2_res);
                let (sr2, ss2) = calc.solve_from_noon(transit2)?;
                len_tomorrow = day_length(&sr2, &ss2);
                break; // Found a valid day!
            }
            // If we reach here, the day didn't exist (Samoa skip), try the day after
            check_date = d.succ_opt();
        }
    }

    // Calculate solar panel output if configured
    let panel_config = args.solarpanel_size.map(|size| {
        let mut config =
            solar_panel::SolarPanelConfig::new(size, args.solarpanel_tilt, args.solarpanel_azimuth)
                .with_efficiency(args.solarpanel_efficiency)
                .with_linke_turbidity(args.linke_turbidity)
                .with_albedo(args.albedo);

        // Set tracking mode based on CLI flags
        if args.solarpanel_dual_axis {
            config = config.with_tracking_mode(solar_panel::TrackingMode::DualAxis);
        } else if args.solarpanel_horizontal_tracking {
            config = config.with_tracking_mode(solar_panel::TrackingMode::HorizontalAxis);
        } else if args.solarpanel_vertical_tracking {
            config = config.with_tracking_mode(solar_panel::TrackingMode::VerticalAxis);
        }

        config
    });

    let day_of_year = date.ordinal();
    // Construct the full NaiveDateTime first, then resolve timezone once.
    let target_dt: Option<DateTime<Tz>> = if let Some(at) = args.at.as_deref() {
        if at == "now" {
            Some(if args.utc {
                Utc::now().with_timezone(&Tz::UTC)
            } else {
                Utc::now().with_timezone(&tz)
            })
        } else {
            let (h, m, s, ns) = parse_time_ns(at)?;

            // 1. Construct the full naive time including nanoseconds
            let naive_dt =
                date.date_naive().and_hms_nano_opt(h, m, s, ns).ok_or("Invalid time digits")?;

            // 2. Resolve against timezone exactly once
            match tz.from_local_datetime(&naive_dt) {
                chrono::LocalResult::Single(t) => Some(t),
                chrono::LocalResult::Ambiguous(t1, t2) => {
                    // Handle Fall Back (e.g. 03:00 happens twice).
                    // t1 is usually the earlier instant (Daylight time), t2 is Standard time.
                    eprintln!(
                        "Note: Time {} is ambiguous (DST transition). Using early option: {} (vs {})",
                        at,
                        t1.format("%H:%M:%S %Z"),
                        t2.format("%H:%M:%S %Z")
                    );
                    Some(t1)
                }
                chrono::LocalResult::None => {
                    // Handle Spring Forward (e.g. 02:30 doesn't exist).
                    // Return the error through the normal `?` path so any
                    // cleanup / destructors in main's scope get to run, rather
                    // than aborting the process abruptly.
                    return Err(format!(
                        "The time {} does not exist on {} (DST gap/Spring Forward).",
                        at,
                        date.date_naive()
                    )
                    .into());
                }
            }
        }
    } else {
        None
    };

    // Determine which position to use for the "Current Conditions" solar panel output
    let (panel_output_pos, panel_output_label) = if let Some(dt) = target_dt {
        (calc.position(dt)?, format!("at {}", dt.format("%H:%M:%S")))
    } else if args.date.is_some() {
        // If user set a date but no time, default to Solar Noon for that day
        (transit_pos, format!("at solar noon ({})", transit.format("%H:%M:%S")))
    } else {
        // Default to now
        (now_pos, "now".to_string())
    };

    let current_panel_output = panel_config.as_ref().map(|config| {
        solar_panel::calculate_output(
            config,
            panel_output_pos.elevation_angle(),
            panel_output_pos.azimuth(),
            args.altitude,
            day_of_year,
        )
    });

    // 1. Current UV (at target time)
    let month = date.month();

    // Fire the Antarctic ozone-hole caveat exactly once (if applicable),
    // rather than inside uv::calculate_uv_index which runs in hot loops.
    if uv::is_ozone_hole_season(args.latitude, month) {
        eprintln!(
            "Warning: Significant ozone depletion possible at lat {:.1} during month {}. \
             Climatology may overestimate UVI.",
            args.latitude, month
        );
    }

    let uv_current = uv::calculate_uv_index(
        panel_output_pos.elevation_angle(),
        args.latitude,
        month,
        args.altitude,
    );

    // 2. Max UV (at Solar Noon)
    let uv_max =
        uv::calculate_uv_index(transit_pos.elevation_angle(), args.latitude, month, args.altitude);

    let uv_data = Some((uv_current, uv_max));

    // Calculate sun positions for the day (used by both daily energy and optimization)
    let sun_positions: Option<Vec<(f64, f64, f64)>> = {
        // 1. Determine integration window (start_dt, end_dt)
        let (start_dt, end_dt) = match transit_res {
            SunriseResult::RegularDay { .. } => match (&sr, &ss) {
                (Some((sr_t, _)), Some((ss_t, _))) => (*sr_t, *ss_t),
                _ => (date, date),
            },
            SunriseResult::AllDay { .. } => {
                let start = find_start_of_day(date.date_naive());
                (start, start + Duration::hours(24))
            }
            SunriseResult::AllNight { .. } => (date, date),
        };

        if start_dt == end_dt {
            None
        } else {
            // 2. Determine "Midnight Reference" for X-axis calculation
            let midnight_ref = find_start_of_day(date.date_naive());

            let step = Duration::minutes(1);
            let mut positions: Vec<(f64, f64, f64)> = Vec::new();
            let mut current_time = start_dt;

            // Loop until strictly before end_dt
            while current_time < end_dt {
                let pos = calc.position(current_time)?;

                // Calculate hours from midnight reference
                let duration_from_midnight = current_time.signed_duration_since(midnight_ref);
                let hours_relative = duration_from_midnight.num_seconds() as f64 / 3600.0;

                // Normalize hours (0..24 range)
                let hours_normalized =
                    if hours_relative < 0.0 { hours_relative + 24.0 } else { hours_relative };

                positions.push((hours_normalized, pos.elevation_angle(), pos.azimuth()));
                current_time += step;
            }

            // Always add the exact end point (Sunset or 24:00)
            let pos_end = calc.position(end_dt)?;
            let duration_end = end_dt.signed_duration_since(midnight_ref);
            let hours_end = duration_end.num_seconds() as f64 / 3600.0;
            let hours_end_norm = if hours_end < 0.0 { hours_end + 24.0 } else { hours_end };

            positions.push((hours_end_norm, pos_end.elevation_angle(), pos_end.azimuth()));

            Some(positions)
        }
    };

    // Calculate daily energy
    let daily_energy = panel_config.as_ref().and_then(|config| {
        if matches!(transit_res, SunriseResult::AllNight { .. }) {
            return Some(0.0);
        }
        sun_positions.as_ref().map(|positions| {
            solar_panel::calculate_daily_energy(
                config,
                args.altitude,
                day_of_year,
                positions.iter().copied(),
            )
        })
    });

    let optimization_result = if args.solarpanel_find_optimum && args.solarpanel_size.is_some() {
        // Build constraints from CLI arguments
        let constraints = optimize::OptimizationConstraints::default()
            .with_tilt_range(args.solarpanel_tilt_range)
            .with_azimuth_range(args.solarpanel_azimuth_range);

        panel_config.as_ref().and_then(|cfg| {
            sun_positions.as_ref().map(|positions| {
                if args.solarpanel_horizontal_tracking {
                    // HSAT: optimize tilt only (azimuth is tracked)
                    optimize::optimize_hsat_tilt(
                        *cfg,
                        args.altitude,
                        day_of_year,
                        positions,
                        &constraints,
                    )
                } else if args.solarpanel_vertical_tracking {
                    // VSAT: optimize azimuth only (tilt is tracked)
                    optimize::optimize_vsat_azimuth(
                        *cfg,
                        args.altitude,
                        day_of_year,
                        positions,
                        &constraints,
                    )
                } else {
                    // Fixed panel: optimize both tilt and azimuth
                    optimize::optimize_fixed_panel_constrained(
                        *cfg,
                        args.altitude,
                        day_of_year,
                        positions,
                        &constraints,
                    )
                }
            })
        })
    } else {
        None
    };

    let yearly_optimization_result = None;
    if let Some(num_adjustments) = args.solarpanel_yearly_adjustments {
        // Build constraints from CLI arguments
        let constraints = optimize::OptimizationConstraints::default()
            .with_tilt_range(args.solarpanel_tilt_range)
            .with_azimuth_range(args.solarpanel_azimuth_range);

        // Create the sun position generator closure
        // We need to move all values into the closure
        let calc_clone = calc;
        let tz_clone = tz;
        let year = date.year();

        let sun_pos_generator: optimize::SunPositionGenerator =
            Box::new(move |_lat, _lon, _alt, day_of_year| {
                let find_start_of_day =
                    |d: chrono::NaiveDate, fallback_dt: DateTime<Tz>| -> DateTime<Tz> {
                        time::start_of_day(tz_clone, d, fallback_dt)
                    };

                // Convert day_of_year to a date
                let target_date = chrono::NaiveDate::from_yo_opt(year, day_of_year)
                    .unwrap_or_else(|| chrono::NaiveDate::from_yo_opt(year, 1).unwrap());

                // Get solar transit for this day
                let day_dt = match tz_clone
                    .from_local_datetime(&target_date.and_hms_opt(12, 0, 0).unwrap())
                {
                    chrono::LocalResult::Single(t) => t,
                    chrono::LocalResult::Ambiguous(t, _) => t,
                    chrono::LocalResult::None => {
                        match tz_clone
                            .from_local_datetime(&target_date.and_hms_opt(13, 0, 0).unwrap())
                        {
                            chrono::LocalResult::Single(t) => t,
                            chrono::LocalResult::Ambiguous(t, _) => t,
                            chrono::LocalResult::None => {
                                // Fallback: use UTC for this day
                                let naive = target_date.and_hms_opt(12, 0, 0).unwrap();
                                DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc)
                                    .with_timezone(&tz_clone)
                            }
                        }
                    }
                };

                // Get transit for this day
                let transit_res = calc_clone.get_transit(day_dt);
                if transit_res.is_none() {
                    return Vec::new();
                }

                let transit = calc_clone.extract_transit_time(&transit_res.unwrap());
                let (sr, ss) = match calc_clone.solve_from_noon(transit) {
                    Ok(pair) => pair,
                    Err(_) => return Vec::new(),
                };

                // Determine start and end times
                let (start_dt, end_dt) = match (&sr, &ss) {
                    (Some((sr_t, _)), Some((ss_t, _))) => (*sr_t, *ss_t),
                    _ => {
                        // Polar day/night handling
                        let start = find_start_of_day(target_date, day_dt);
                        (start, start + Duration::hours(24))
                    }
                };

                if start_dt >= end_dt {
                    return Vec::new();
                }

                // Midnight reference for this day
                let midnight_ref = find_start_of_day(target_date, day_dt);

                let step = Duration::minutes(10); // Use 10-minute steps for efficiency
                let mut positions: Vec<(f64, f64, f64)> = Vec::new();
                let mut current_time = start_dt;

                while current_time < end_dt {
                    // On SPA failure we abort this day and return what we have so
                    // far; yearly optimization sums over sampled days, so a single
                    // truncated/empty day only locally reduces accuracy.
                    let pos = match calc_clone.position(current_time) {
                        Ok(p) => p,
                        Err(_) => return positions,
                    };
                    let duration_from_midnight = current_time.signed_duration_since(midnight_ref);
                    let hours_relative = duration_from_midnight.num_seconds() as f64 / 3600.0;
                    let hours_normalized =
                        if hours_relative < 0.0 { hours_relative + 24.0 } else { hours_relative };
                    positions.push((hours_normalized, pos.elevation_angle(), pos.azimuth()));
                    current_time += step;
                }

                // Add endpoint
                let pos_end = match calc_clone.position(end_dt) {
                    Ok(p) => p,
                    Err(_) => return positions,
                };
                let duration_end = end_dt.signed_duration_since(midnight_ref);
                let hours_end = duration_end.num_seconds() as f64 / 3600.0;
                let hours_end_norm = if hours_end < 0.0 { hours_end + 24.0 } else { hours_end };
                positions.push((hours_end_norm, pos_end.elevation_angle(), pos_end.azimuth()));

                positions
            });

        panel_config.as_ref().map(|cfg| {
            let result = if args.solarpanel_horizontal_tracking {
                // HSAT: optimize tilt adjustments with azimuth tracking
                optimize::optimize_yearly_hsat(
                    *cfg,
                    args.latitude,
                    args.longitude,
                    args.altitude,
                    year,
                    num_adjustments,
                    &constraints,
                    &sun_pos_generator,
                )
            } else {
                // Fixed panel: optimize tilt adjustments with fixed azimuth
                optimize::optimize_yearly_adjustments(
                    *cfg,
                    args.latitude,
                    args.longitude,
                    args.altitude,
                    year,
                    num_adjustments,
                    &constraints,
                    &sun_pos_generator,
                )
            };
            (result, constraints, args.solarpanel_horizontal_tracking)
        })
    } else {
        None
    };

    if use_argos {
        output::print_argos(
            &now_pos,
            &sr,
            &ss,
            &transit_res,
            transit,
            &transit_pos,
            len_today,
            len_tomorrow,
            current_panel_output.as_ref(),
            daily_energy,
            uv_data,
        );
        return Ok(());
    }

    println!("Location : lat={:.6}, lon={:.6}", args.latitude, args.longitude);
    println!("Timezone : {}", tz);
    println!("Date     : {}", date.date_naive());
    println!("Target altitude : {:.6}°", calc.target);
    println!();

    if let Some(dt) = target_dt {
        let pos = calc.position(dt)?;
        output::print_sun_position_at_time(dt, &pos, uv_current, uv_max);

        if let Some(ref config) = panel_config {
            let at_output = solar_panel::calculate_output(
                config,
                pos.elevation_angle(),
                pos.azimuth(),
                args.altitude,
                day_of_year,
            );
            output::print_solar_panel_at_time(dt, &at_output);
        }

        println!();
    }

    output::print_sun_events(
        &sr,
        &ss,
        &transit_res,
        transit,
        &transit_pos,
        len_today,
        len_tomorrow,
        &calc,
        date,
    );

    if target_dt.is_none() {
        output::print_uv_info(&panel_output_label, uv_current, uv_max);
    }

    if let Some(ref config) = panel_config
        && let Some(ref output) = current_panel_output
    {
        output::print_solar_panel_info(config, output, &panel_output_label, daily_energy);
    }

    // Print optimization results if requested
    if let Some(ref opt_result) = optimization_result {
        println!();

        if args.solarpanel_horizontal_tracking {
            println!("=== Optimal HSAT Configuration ===");
            println!("Optimized tilt for horizontal single-axis tracking on {}", date.date_naive());
            println!();
            println!("Optimal tilt    : {:6.1}°", opt_result.tilt_deg);
            println!("Daily energy    : {}", solar_panel::format_energy(opt_result.energy_wh));
        } else if args.solarpanel_vertical_tracking {
            println!("=== Optimal VSAT Configuration ===");
            println!(
                "Optimized azimuth for vertical single-axis tracking on {}",
                date.date_naive()
            );
            println!();
            println!("Optimal azimuth : {:6.1}°", opt_result.azimuth_deg);
            println!("Daily energy    : {}", solar_panel::format_energy(opt_result.energy_wh));

            // Compare with current azimuth
            if let Some(current_energy) = daily_energy {
                let improvement = opt_result.energy_wh - current_energy;
                let improvement_pct =
                    if current_energy > 0.0 { (improvement / current_energy) * 100.0 } else { 0.0 };

                if improvement > 0.1 && improvement_pct > 0.1 {
                    println!(
                        "vs current ({:.1}° az): {} ({:+.1}%)",
                        args.solarpanel_azimuth,
                        solar_panel::format_energy(current_energy),
                        improvement_pct
                    );
                }
            }
        } else {
            println!("=== Optimal Fixed Panel Orientation ===");
            println!("Optimized for maximum daily energy on {}", date.date_naive());
            println!();

            // Check if optimal is essentially the same as user-specified
            let tilt_diff = (opt_result.tilt_deg - args.solarpanel_tilt).abs();
            let azimuth_diff = {
                let diff = (opt_result.azimuth_deg - args.solarpanel_azimuth).abs();
                // Handle wrap-around (e.g., 359° vs 1°)
                diff.min(360.0 - diff)
            };

            let is_same_as_user = tilt_diff < 0.15 && azimuth_diff < 0.15;

            if is_same_as_user {
                println!(
                    "Your configuration ({:.1}° tilt, {:.1}° azimuth) is already optimal!",
                    args.solarpanel_tilt, args.solarpanel_azimuth
                );
                println!("Daily energy    : {}", solar_panel::format_energy(opt_result.energy_wh));
            } else {
                println!("Optimal tilt    : {:6.1}°", opt_result.tilt_deg);
                println!("Optimal azimuth : {:6.1}°", opt_result.azimuth_deg);
                println!("Daily energy    : {}", solar_panel::format_energy(opt_result.energy_wh));
                println!();

                // Compare with current configuration
                if let Some(current_energy) = daily_energy {
                    let improvement = opt_result.energy_wh - current_energy;
                    let improvement_pct = if current_energy > 0.0 {
                        (improvement / current_energy) * 100.0
                    } else {
                        0.0
                    };

                    if improvement.abs() < 0.1 || improvement_pct.abs() < 0.1 {
                        println!(
                            "Your configuration ({:.1}° tilt, {:.1}° az) produces the same energy.",
                            args.solarpanel_tilt, args.solarpanel_azimuth
                        );
                    } else if improvement > 0.0 {
                        println!(
                            "vs current ({:.1}° tilt, {:.1}° az): {} ({:+.1}%)",
                            args.solarpanel_tilt,
                            args.solarpanel_azimuth,
                            solar_panel::format_energy(current_energy),
                            improvement_pct
                        );
                    } else {
                        // User config is actually better (shouldn't happen with correct optimization)
                        println!(
                            "Note: Your configuration ({:.1}° tilt, {:.1}° az) produces {} - same as optimal.",
                            args.solarpanel_tilt,
                            args.solarpanel_azimuth,
                            solar_panel::format_energy(current_energy),
                        );
                    }
                }
            }

            // Calculate and show tracking comparison
            if let (Some(positions), Some(panel_size)) = (&sun_positions, args.solarpanel_size) {
                let tracking_config = solar_panel::SolarPanelConfig::new(
                    panel_size, 0.0, // Ignored for dual-axis tracking
                    0.0, // Ignored for dual-axis tracking
                )
                .with_efficiency(args.solarpanel_efficiency)
                .with_linke_turbidity(args.linke_turbidity)
                .with_albedo(args.albedo)
                .with_tracking_mode(solar_panel::TrackingMode::DualAxis);

                let tracking_energy = solar_panel::calculate_daily_energy(
                    &tracking_config,
                    args.altitude,
                    day_of_year,
                    positions.iter().copied(),
                );

                let tracking_vs_optimal = ((tracking_energy / opt_result.energy_wh) - 1.0) * 100.0;
                println!(
                    "vs dual-axis    : {} ({:+.1}%)",
                    solar_panel::format_energy(tracking_energy),
                    tracking_vs_optimal
                );
            }
        }

        println!();
        println!("(Optimization used {} evaluations)", opt_result.evaluations);
    }

    // Print yearly optimization results if requested
    if let Some((ref yearly_result, ref constraints, is_hsat)) = yearly_optimization_result {
        output::print_yearly_optimization(yearly_result, date.year(), constraints, is_hsat);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    #[test]
    fn test_conversion_preserves_timezone() {
        // 1. Setup: Create a time in New York (UTC-5)
        // 2023-01-01 12:00:00 UTC is 07:00:00 in New York
        let ny_tz: chrono_tz::Tz = "America/New_York".parse().unwrap();
        let anchor = Utc.timestamp_opt(1672574400, 0).unwrap().with_timezone(&ny_tz);
        assert_eq!(anchor.to_string(), "2023-01-01 07:00:00 EST");
        let result = parse_relative_to_chrono(&anchor, "tomorrow");
        assert!(result.is_ok());
        let val = result.unwrap();
        // The time shifted by 1 day
        assert_eq!(val.to_string(), "2023-01-02 07:00:00 EST");
        // The resulting Chrono object is still in New York time, not UTC
        assert_eq!(val.timezone(), ny_tz);
    }

    #[test]
    fn test_handles_utc_gracefully() {
        let anchor = Utc::now();
        let result = parse_relative_to_chrono(&anchor, "tomorrow");

        assert!(result.is_ok());
        let val = result.unwrap();
        // Ensure it's still UTC
        assert_eq!(val.timezone(), Utc);
        // Ensure time moved forward
        assert!(val > anchor);
    }

    #[test]
    fn test_propagates_parsing_errors() {
        let anchor = Utc::now();
        let result = parse_relative_to_chrono(&anchor, "shallnotpass");

        assert!(result.is_err());
        let err_msg = result.err().unwrap();
        assert!(err_msg.contains("Failed to parse"));
    }
}
