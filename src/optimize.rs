#![allow(clippy::too_many_arguments, clippy::type_complexity)]

use crate::solar_panel::{self, SolarPanelConfig, TrackingMode};

// Optimization constraints matching the original grid resolution
const TILT_DEGENERATE: f64 = 0.5;
const TILT_GRID_STEP: f64 = 0.5;
const AZ_GRID_STEP: f64 = 1.0;

/// Tilt used to probe the azimuth objective when the current best tilt is
/// degenerate. A flat panel collects the same energy whichever way it faces,
/// so searching azimuth at a degenerate tilt optimizes a constant function.
/// Any clearly non-degenerate tilt works; 45° keeps the panel's response to
/// azimuth strong at both high and low sun.
const AZ_PROBE_TILT: f64 = 45.0;

/// Constraints for panel optimization (tilt and azimuth ranges)
///
/// If `azimuth_min > azimuth_max` the azimuth range wraps through north,
/// e.g. (300, 60) allows NW through NE.
#[derive(Debug, Clone, Copy)]
pub struct OptimizationConstraints {
    /// Minimum allowed tilt (degrees)
    pub tilt_min: f64,
    /// Maximum allowed tilt (degrees)
    pub tilt_max: f64,
    /// Minimum allowed azimuth (degrees)
    pub azimuth_min: f64,
    /// Maximum allowed azimuth (degrees)
    pub azimuth_max: f64,
    /// True if the user explicitly set a non-trivial azimuth range.
    ///
    /// Gates the azimuth search strategy: a full 0..360 circle (default, or
    /// explicitly passed) uses a cardinal-quadrant scan to handle 0°/360°
    /// wraparound, while a user-set narrower range uses direct golden-section
    /// search inside it (the range itself may wrap, see struct docs).
    ///
    /// Tilt has no wraparound, so no corresponding flag is needed.
    azimuth_explicit: bool,
}

impl Default for OptimizationConstraints {
    fn default() -> Self {
        Self {
            tilt_min: 0.0,
            tilt_max: 90.0,
            azimuth_min: 0.0,
            azimuth_max: 360.0,
            azimuth_explicit: false,
        }
    }
}

impl OptimizationConstraints {
    pub fn with_tilt_range(mut self, range: Option<(f64, f64)>) -> Self {
        if let Some((min, max)) = range {
            // Accept the endpoints in either order: tilt has no wraparound
            // semantics and an inverted pair would later panic in f64::clamp
            // (which requires min <= max)
            let (min, max) = if min <= max { (min, max) } else { (max, min) };
            self.tilt_min = min.clamp(0.0, 90.0);
            self.tilt_max = max.clamp(0.0, 90.0);
        }
        self
    }

    pub fn with_azimuth_range(mut self, range: Option<(f64, f64)>) -> Self {
        if let Some((min, max)) = range {
            self.azimuth_min = min.clamp(0.0, 360.0);
            self.azimuth_max = max.clamp(0.0, 360.0);
            // min > max means the range wraps through north (e.g. 300..60).
            // An explicit full circle is the same as no constraint at all, so
            // keep the wraparound-aware quadrant search for it.
            self.azimuth_explicit = !(self.azimuth_min <= 0.0 && self.azimuth_max >= 360.0);
        }
        self
    }

    /// True if the azimuth range is effectively the full 0..360 circle
    /// (left at the default, or explicitly set to the full circle).
    /// Callers use this to pick the azimuth search algorithm.
    pub fn azimuth_unconstrained(&self) -> bool {
        !self.azimuth_explicit
    }

    /// Azimuth range as a continuous interval for searching; the upper bound
    /// exceeds 360 when the range wraps through north (e.g. 300..60 becomes
    /// 300..420). Angles taken from this interval must be evaluated and
    /// reported modulo 360.
    fn azimuth_interval(&self) -> (f64, f64) {
        if self.azimuth_min <= self.azimuth_max {
            (self.azimuth_min, self.azimuth_max)
        } else {
            (self.azimuth_min, self.azimuth_max + 360.0)
        }
    }

    /// Center of the (possibly wrapped) azimuth range, normalized to [0, 360)
    fn azimuth_center(&self) -> f64 {
        let (lo, hi) = self.azimuth_interval();
        ((lo + hi) / 2.0).rem_euclid(360.0)
    }

    /// True if `az` (any real angle) points inside the allowed azimuth range.
    /// Both the angle and the range are compared on the circle, so a range
    /// ending at 360 contains 0 and a range wrapping through north works too.
    pub fn contains_azimuth(&self, az: f64) -> bool {
        if !self.azimuth_explicit {
            return true;
        }
        let (lo, hi) = self.azimuth_interval();
        let az = az.rem_euclid(360.0);
        (az >= lo && az <= hi) || (az + 360.0 >= lo && az + 360.0 <= hi)
    }

    /// Nearest allowed azimuth to `az`, measured as an angle so the two range
    /// endpoints compete across the 0/360 seam.
    ///
    /// Used where an azimuth is picked by rule rather than searched (the
    /// yearly schedulers face the panel at the equator) so that a user-set
    /// azimuth range is still honoured instead of silently overridden.
    pub fn clamp_azimuth(&self, az: f64) -> f64 {
        if self.contains_azimuth(az) {
            return az.rem_euclid(360.0);
        }
        let (lo, hi) = self.azimuth_interval();
        let separation = |edge: f64| {
            let d = (az - edge).rem_euclid(360.0);
            d.min(360.0 - d)
        };
        if separation(lo) <= separation(hi) { lo.rem_euclid(360.0) } else { hi.rem_euclid(360.0) }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PanelOptimum {
    pub tilt_deg: f64,
    pub azimuth_deg: f64,
    pub energy_wh: f64,
    pub evaluations: usize,
}

/// Result for a single adjustment period in yearly optimization
#[derive(Debug, Clone)]
pub struct AdjustmentPeriod {
    /// Day of year when this adjustment starts (1-366)
    pub start_day: u32,
    /// Day of year when this adjustment ends (1-366)
    pub end_day: u32,
    /// Optimal tilt for this period
    pub tilt_deg: f64,
    /// Optimal azimuth for this period
    pub azimuth_deg: f64,
    /// Total energy produced during this period (Wh)
    pub period_energy_wh: f64,
}

/// Complete yearly optimization result
#[derive(Debug, Clone)]
pub struct YearlyOptimizationResult {
    /// List of adjustment periods with optimal settings
    pub periods: Vec<AdjustmentPeriod>,
    /// Total yearly energy (Wh)
    pub total_energy_wh: f64,
    /// Energy with fixed optimal settings (no adjustments)
    pub fixed_optimal_energy_wh: f64,
    /// Total optimization evaluations
    pub evaluations: usize,
}

/// Callback type for sun position generation
/// Takes (latitude, longitude, altitude, day_of_year) and returns sun positions
pub type SunPositionGenerator =
    Box<dyn Fn(f64, f64, f64, u32) -> Vec<(f64, f64, f64)> + Send + Sync>;

#[allow(dead_code)]
fn optimize_fixed_panel(
    base: SolarPanelConfig,
    altitude_m: f64,
    day_of_year: u32,
    sun_positions: &[(f64, f64, f64)],
) -> PanelOptimum {
    optimize_fixed_panel_constrained(
        base,
        altitude_m,
        day_of_year,
        sun_positions,
        &OptimizationConstraints::default(),
    )
}

/// Sun azimuth at the highest sampled elevation: a hemisphere-agnostic initial
/// guess for the optimal panel azimuth (~180° north of the tropics, ~0/360°
/// south of them). Falls back to 180° if the sun never rises in the sample.
fn culmination_azimuth(sun_positions: &[(f64, f64, f64)]) -> f64 {
    sun_positions
        .iter()
        .filter(|p| p.1 > 0.0)
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|p| p.2)
        .unwrap_or(180.0)
}

/// Optimize panel orientation with constraints
pub fn optimize_fixed_panel_constrained(
    base: SolarPanelConfig,
    altitude_m: f64,
    day_of_year: u32,
    sun_positions: &[(f64, f64, f64)],
    constraints: &OptimizationConstraints,
) -> PanelOptimum {
    let mut evaluations = 0;

    // `is_unconstrained` here gates the azimuth search strategy only — tilt has
    // no wraparound, so it's always GSS on the (possibly constrained) range.
    let is_unconstrained = constraints.azimuth_unconstrained();

    // Helper closure to calculate energy and count evaluations
    // We capture specific variables to avoid cloning the whole context repeatedly
    // Azimuth is normalized modulo 360 so searches may probe outside [0, 360)
    // (wraparound ranges and the ±90° quadrant window need this)
    let mut calc_energy = |tilt: f64, az: f64| -> f64 {
        evaluations += 1;
        let cfg = SolarPanelConfig { tilt_deg: tilt, azimuth_deg: az.rem_euclid(360.0), ..base };
        solar_panel::calculate_daily_energy(
            &cfg,
            altitude_m,
            day_of_year,
            sun_positions.iter().copied(),
        )
    };

    // --------------------------------------------------
    // 1) Optimize tilt within constraints
    // --------------------------------------------------
    // For the unconstrained case, start from the sun's culmination azimuth:
    // a hard-coded 180° tunes the tilt for a panel facing away from the sun
    // in the southern hemisphere, and the tilt is what the azimuth search
    // builds on. For the constrained case, use the center of the azimuth range.
    let initial_az = if is_unconstrained {
        culmination_azimuth(sun_positions)
    } else {
        constraints.azimuth_center()
    };

    // GSS on constrained tilt range
    let (raw_tilt, _) =
        golden_section_search(constraints.tilt_min, constraints.tilt_max, 0.4, |t| {
            calc_energy(t, initial_az)
        });

    // Snap to the original grid resolution (0.5 degrees)
    let mut best_tilt = (raw_tilt / TILT_GRID_STEP).round() * TILT_GRID_STEP;
    best_tilt = best_tilt.clamp(constraints.tilt_min, constraints.tilt_max);

    let mut best_az = initial_az;

    // --------------------------------------------------
    // 2) Optimize azimuth
    // --------------------------------------------------
    // A degenerate tilt must not suppress the azimuth search. It only means
    // the panel is best off flat *at the azimuth tried so far*, which for a
    // constrained range is the range center — and the center can easily be a
    // direction worth lying flat in while a tilted panel elsewhere in the same
    // range collects much more (a 60..120 range under a low sun lost ~19% of
    // the available energy this way). Azimuth does have no effect on a flat
    // panel though, so probe the azimuth objective at a non-degenerate tilt
    // and re-tune the tilt afterwards for whichever azimuth that finds.
    if constraints.tilt_max >= TILT_DEGENERATE {
        let az_probe_tilt = if best_tilt >= TILT_DEGENERATE {
            best_tilt
        } else {
            AZ_PROBE_TILT.clamp(constraints.tilt_min, constraints.tilt_max)
        };

        if is_unconstrained {
            // Original algorithm: Coarse Quadrant Search then fine GSS
            // Because azimuth wraps (0 == 360), we check cardinal directions first
            // to identify the correct sector and avoid getting stuck at a boundary.
            let mut best_sector_az = initial_az;
            let mut max_sector_e = calc_energy(az_probe_tilt, initial_az);

            for &az_check in &[0.0, 90.0, 180.0, 270.0] {
                let e = calc_energy(az_probe_tilt, az_check);
                if e > max_sector_e {
                    max_sector_e = e;
                    best_sector_az = az_check;
                }
            }

            // Fine GSS within +/- 90 degrees of the best sector
            // (may probe outside [0, 360); calc_energy normalizes)
            let search_min = best_sector_az - 90.0;
            let search_max = best_sector_az + 90.0;

            let (raw_az, _) = golden_section_search(search_min, search_max, 0.8, |a| {
                calc_energy(az_probe_tilt, a)
            });

            // Snap to grid resolution (1.0 degree)
            best_az = (raw_az / AZ_GRID_STEP).round() * AZ_GRID_STEP;

            // Normalize to [0, 360) for final output consistency
            best_az = best_az.rem_euclid(360.0);
        } else {
            // Constrained case: do GSS directly within the range. The interval
            // is continuous even for ranges wrapping through north (300..60
            // becomes 300..420) so lo <= hi always holds and clamp is safe.
            let (az_lo, az_hi) = constraints.azimuth_interval();

            if az_hi - az_lo > AZ_GRID_STEP {
                let (raw_az, _) =
                    golden_section_search(az_lo, az_hi, 0.8, |a| calc_energy(az_probe_tilt, a));

                // Snap to grid resolution (1.0 degree)
                best_az = (raw_az / AZ_GRID_STEP).round() * AZ_GRID_STEP;
                best_az = best_az.clamp(az_lo, az_hi);

                // Normalize to [0, 360) for final output consistency
                best_az = best_az.rem_euclid(360.0);
            }
        }

        // --------------------------------------------------
        // 3) Re-optimize tilt at the chosen azimuth
        // --------------------------------------------------
        // Tilt and azimuth interact; a single sequential pass leaves the tilt
        // tuned for the initial azimuth guess (a coordinate-descent trap that
        // cost ~10% energy when the initial azimuth was far from optimal).
        let (raw_tilt2, _) =
            golden_section_search(constraints.tilt_min, constraints.tilt_max, 0.4, |t| {
                calc_energy(t, best_az)
            });
        best_tilt = (raw_tilt2 / TILT_GRID_STEP).round() * TILT_GRID_STEP;
        best_tilt = best_tilt.clamp(constraints.tilt_min, constraints.tilt_max);
    } else {
        // The tilt range itself pins the panel flat → azimuth physically
        // undefined. For unconstrained: use 180° (south).
        // For constrained: use center of range.
        best_az = if is_unconstrained { 180.0 } else { constraints.azimuth_center() };
    }

    // 2b) Canonicalize azimuth for near-horizontal panels
    // Matches original logic: if tilt is very low, force canonical value
    if best_tilt <= 5.0 {
        best_az = if is_unconstrained { 180.0 } else { constraints.azimuth_center() };
    }

    // Final energy at the exact orientation we return, so energy_wh always
    // corresponds to (tilt_deg, azimuth_deg) even after canonicalization
    let best_energy = calc_energy(best_tilt, best_az);

    PanelOptimum {
        // Round to 1 decimal place for display consistency (matches original behavior)
        // Use + 0.0 to normalize -0.0 to 0.0
        tilt_deg: (best_tilt * 10.0).round() / 10.0 + 0.0,
        azimuth_deg: (best_az * 10.0).round() / 10.0 + 0.0,
        energy_wh: best_energy,
        evaluations,
    }
}

/// Optimize tilt for HSAT (Horizontal Single-Axis Tracking)
/// HSAT tracks the sun by rotating about its axis automatically, so we only
/// need to find the optimal axis tilt. The axis lean direction is taken from
/// `base.azimuth_deg` (the rest-position facing, typically 180° toward the
/// equator in the northern hemisphere, 0° in the southern).
pub fn optimize_hsat_tilt(
    base: SolarPanelConfig,
    altitude_m: f64,
    day_of_year: u32,
    sun_positions: &[(f64, f64, f64)],
    constraints: &OptimizationConstraints,
) -> PanelOptimum {
    let mut evaluations = 0;

    // For HSAT, the panel rotation tracks the sun, so we create a config with
    // HSAT mode and only optimize the axis tilt
    let mut calc_energy = |tilt: f64| -> f64 {
        evaluations += 1;
        let cfg = SolarPanelConfig {
            tilt_deg: tilt,
            azimuth_deg: base.azimuth_deg, // Axis lean direction (rest-position facing)
            tracking_mode: TrackingMode::HorizontalAxis,
            ..base
        };
        solar_panel::calculate_daily_energy(
            &cfg,
            altitude_m,
            day_of_year,
            sun_positions.iter().copied(),
        )
    };

    // GSS on constrained tilt range, kept inside the range the tracker model
    // accepts: it clamps an axis tilt of 90° down to near-vertical internally,
    // so searching past that would report a tilt the energy figure is not for
    let tilt_max = constraints.tilt_max.min(solar_panel::MAX_TRACKER_AXIS_TILT_DEG);
    let tilt_min = constraints.tilt_min.min(tilt_max);
    let (raw_tilt, _) = golden_section_search(tilt_min, tilt_max, 0.4, &mut calc_energy);

    let best_tilt = (raw_tilt / TILT_GRID_STEP).round() * TILT_GRID_STEP;
    let best_tilt = best_tilt.clamp(tilt_min, tilt_max);
    let best_energy = calc_energy(best_tilt);

    PanelOptimum {
        tilt_deg: (best_tilt * 10.0).round() / 10.0 + 0.0,
        azimuth_deg: base.azimuth_deg, // Axis lean direction used during optimization
        energy_wh: best_energy,
        evaluations,
    }
}

/// Optimize tilt for VSAT (Vertical Single-Axis Tracking)
/// VSAT tracks the sun's azimuth automatically (panel spins about a vertical
/// axis), so we only need to find the optimal fixed tilt.
pub fn optimize_vsat_tilt(
    base: SolarPanelConfig,
    altitude_m: f64,
    day_of_year: u32,
    sun_positions: &[(f64, f64, f64)],
    constraints: &OptimizationConstraints,
) -> PanelOptimum {
    let mut evaluations = 0;

    // For VSAT, azimuth tracks the sun, so we create a config with VSAT mode
    // and only optimize tilt
    let mut calc_energy = |tilt: f64| -> f64 {
        evaluations += 1;
        let cfg = SolarPanelConfig {
            tilt_deg: tilt,
            azimuth_deg: base.azimuth_deg, // Ignored for VSAT - azimuth tracks sun
            tracking_mode: TrackingMode::VerticalAxis,
            ..base
        };
        solar_panel::calculate_daily_energy(
            &cfg,
            altitude_m,
            day_of_year,
            sun_positions.iter().copied(),
        )
    };

    // GSS on constrained tilt range
    let (raw_tilt, _) =
        golden_section_search(constraints.tilt_min, constraints.tilt_max, 0.4, &mut calc_energy);

    let best_tilt = (raw_tilt / TILT_GRID_STEP).round() * TILT_GRID_STEP;
    let best_tilt = best_tilt.clamp(constraints.tilt_min, constraints.tilt_max);
    let best_energy = calc_energy(best_tilt);

    PanelOptimum {
        tilt_deg: (best_tilt * 10.0).round() / 10.0 + 0.0,
        azimuth_deg: 0.0, // Not applicable for VSAT - azimuth tracks sun
        energy_wh: best_energy,
        evaluations,
    }
}

/// Optimize yearly panel output with a given number of adjustment periods
///
/// This function divides the year into `num_adjustments` periods and finds
/// the optimal tilt for each period to maximize total yearly energy.
///
/// **Important**: In real-world seasonal adjustment systems, only tilt is adjusted.
/// Azimuth remains fixed (equator-facing, so typically south at 180° in the
/// northern hemisphere, and moved to the nearest allowed heading if the
/// constraints exclude it) because:
/// - The sun's daily E-W movement is captured by exposure time, not panel rotation
/// - Only continuous tracking systems (HSAT and dual-axis) adjust azimuth
/// - Manual seasonal adjustments only change tilt 2-4 times per year
///
/// # Arguments
/// * `base_config` - Base solar panel configuration
/// * `latitude` - Observer latitude in degrees
/// * `longitude` - Observer longitude in degrees
/// * `altitude_m` - Observer altitude in meters
/// * `year` - Year for calculations (affects leap year handling)
/// * `num_adjustments` - Number of adjustment periods (1-366)
/// * `constraints` - Tilt and azimuth constraints
/// * `sun_pos_generator` - Function to generate sun positions for a given day
///
/// # Returns
/// YearlyOptimizationResult with optimal settings for each period
pub fn optimize_yearly_adjustments(
    base_config: SolarPanelConfig,
    latitude: f64,
    longitude: f64,
    altitude_m: f64,
    year: i32,
    num_adjustments: u32,
    constraints: &OptimizationConstraints,
    sun_pos_generator: &SunPositionGenerator,
) -> YearlyOptimizationResult {
    optimize_yearly_tilt_schedule(
        base_config,
        latitude,
        longitude,
        altitude_m,
        year,
        num_adjustments,
        constraints,
        TrackingMode::Fixed,
        sun_pos_generator,
    )
}

/// Optimize yearly HSAT panel output with a given number of tilt adjustment periods
///
/// HSAT tracks azimuth automatically throughout each day, but tilt can be
/// seasonally adjusted for optimal performance.
///
/// # Arguments
/// * `base_config` - Base solar panel configuration
/// * `latitude` - Observer latitude in degrees
/// * `longitude` - Observer longitude in degrees
/// * `altitude_m` - Observer altitude in meters
/// * `year` - Year for calculations (affects leap year handling)
/// * `num_adjustments` - Number of tilt adjustment periods (1-366)
/// * `constraints` - Tilt constraints
/// * `sun_pos_generator` - Function to generate sun positions for a given day
///
/// # Returns
/// YearlyOptimizationResult with optimal tilt settings for each period
pub fn optimize_yearly_hsat(
    base_config: SolarPanelConfig,
    latitude: f64,
    longitude: f64,
    altitude_m: f64,
    year: i32,
    num_adjustments: u32,
    constraints: &OptimizationConstraints,
    sun_pos_generator: &SunPositionGenerator,
) -> YearlyOptimizationResult {
    optimize_yearly_tilt_schedule(
        base_config,
        latitude,
        longitude,
        altitude_m,
        year,
        num_adjustments,
        constraints,
        TrackingMode::HorizontalAxis,
        sun_pos_generator,
    )
}

/// Shared driver for the yearly tilt schedulers
///
/// `tracking_mode` selects what a period's tilt means: the panel tilt itself
/// for a fixed panel, or the axis tilt for HSAT (whose rotation tracks the sun
/// within each day). Azimuth is not adjusted in either case, so it is chosen
/// once up front.
fn optimize_yearly_tilt_schedule(
    base_config: SolarPanelConfig,
    latitude: f64,
    longitude: f64,
    altitude_m: f64,
    year: i32,
    num_adjustments: u32,
    constraints: &OptimizationConstraints,
    tracking_mode: TrackingMode,
    sun_pos_generator: &SunPositionGenerator,
) -> YearlyOptimizationResult {
    let is_leap = is_leap_year(year);
    let days_in_year: u32 = if is_leap { 366 } else { 365 };

    let num_adjustments = num_adjustments.clamp(1, days_in_year);

    // Calculate period boundaries
    // We divide the year into equal periods
    let periods = calculate_period_boundaries(days_in_year, num_adjustments);

    // Rest-position facing of the panel: toward the equator, so south-facing
    // (180°) in the northern hemisphere and north-facing (0°) in the southern.
    // Projected into the caller's azimuth range when one was set.
    let azimuth = constraints.clamp_azimuth(if latitude >= 0.0 { 180.0 } else { 0.0 });

    // Generate every day's sun positions once and keep them. Beyond saving the
    // generator calls, it lets the fixed-tilt baseline below be measured on
    // exactly the same sample days and weights as the per-period optima.
    let samples: Vec<PeriodSamples> = periods
        .iter()
        .map(|&(start_day, end_day)| {
            PeriodSamples::generate(
                latitude,
                longitude,
                altitude_m,
                start_day,
                end_day,
                sun_pos_generator,
            )
        })
        .collect();

    let mut total_evaluations = 0;
    let mut adjustment_periods = Vec::with_capacity(periods.len());
    let mut total_energy = 0.0;

    // For each period, find the optimal tilt (azimuth stays fixed)
    for (&(start_day, end_day), period) in periods.iter().zip(&samples) {
        let period_result = optimize_tilt_over(
            std::slice::from_ref(period),
            &base_config,
            altitude_m,
            azimuth,
            tracking_mode,
            constraints,
        );

        total_evaluations += period_result.evaluations;
        total_energy += period_result.period_energy_wh;

        adjustment_periods.push(AdjustmentPeriod {
            start_day,
            end_day,
            tilt_deg: period_result.tilt_deg,
            azimuth_deg: azimuth,
            period_energy_wh: period_result.period_energy_wh,
        });
    }

    // Calculate fixed optimal (single tilt setting for entire year).
    // Measured over all the periods' samples, so it is the same quadrature the
    // per-period totals used: a single tilt is one of the choices open to
    // every period, so the adjusted total can no longer come out *below* the
    // never-adjusted one.
    let fixed_result =
        optimize_tilt_over(&samples, &base_config, altitude_m, azimuth, tracking_mode, constraints);
    total_evaluations += fixed_result.evaluations;

    YearlyOptimizationResult {
        periods: adjustment_periods,
        total_energy_wh: total_energy,
        fixed_optimal_energy_wh: fixed_result.period_energy_wh,
        evaluations: total_evaluations,
    }
}

/// Internal structure for period optimization results
struct PeriodOptimum {
    tilt_deg: f64,
    period_energy_wh: f64,
    evaluations: usize,
}

/// Sun positions for one period's sample days, plus the weight that scales a
/// sampled total up to the full period length
struct PeriodSamples {
    /// (day_of_year, sun positions for that day)
    days: Vec<(u32, Vec<(f64, f64, f64)>)>,
    /// period length / number of sample days
    weight: f64,
}

impl PeriodSamples {
    fn generate(
        latitude: f64,
        longitude: f64,
        altitude_m: f64,
        start_day: u32,
        end_day: u32,
        sun_pos_generator: &SunPositionGenerator,
    ) -> Self {
        // Sample representative days within the period
        let sample_days = select_sample_days(start_day, end_day);
        let weight = (end_day - start_day + 1) as f64 / sample_days.len() as f64;

        // Pre-generate sun positions for all sample days
        let days = sample_days
            .into_iter()
            .map(|day| (day, sun_pos_generator(latitude, longitude, altitude_m, day)))
            .collect();

        Self { days, weight }
    }

    /// Energy over the whole period, scaled up from the sampled days
    fn energy(&self, cfg: &SolarPanelConfig, altitude_m: f64) -> f64 {
        let sampled: f64 = self
            .days
            .iter()
            .map(|(day, positions)| {
                solar_panel::calculate_daily_energy(
                    cfg,
                    altitude_m,
                    *day,
                    positions.iter().copied(),
                )
            })
            .sum();
        sampled * self.weight
    }
}

/// Find the single tilt that maximizes total energy across `samples`
///
/// Given one period this is that period's optimal tilt; given every period it
/// is the best never-adjusted tilt for the year. Azimuth is fixed, which
/// matches real-world seasonal adjustment: only tilt gets changed by hand.
fn optimize_tilt_over(
    samples: &[PeriodSamples],
    base_config: &SolarPanelConfig,
    altitude_m: f64,
    azimuth_deg: f64,
    tracking_mode: TrackingMode,
    constraints: &OptimizationConstraints,
) -> PeriodOptimum {
    let mut evaluations = 0;

    // An HSAT axis cannot stand exactly vertical, so keep the search inside
    // the range the tracker model accepts; otherwise the reported axis tilt
    // is not the one the energy figure was computed for
    let tilt_max = match tracking_mode {
        TrackingMode::HorizontalAxis => {
            constraints.tilt_max.min(solar_panel::MAX_TRACKER_AXIS_TILT_DEG)
        }
        _ => constraints.tilt_max,
    };
    let tilt_min = constraints.tilt_min.min(tilt_max);

    // Calculate total energy for a given tilt (azimuth is fixed)
    let mut calc_energy = |tilt: f64| -> f64 {
        evaluations += 1;
        let cfg = SolarPanelConfig { tilt_deg: tilt, azimuth_deg, tracking_mode, ..*base_config };
        samples.iter().map(|period| period.energy(&cfg, altitude_m)).sum()
    };

    // Optimize tilt using Golden Section Search
    let (raw_tilt, _) = golden_section_search(tilt_min, tilt_max, 0.4, &mut calc_energy);

    let best_tilt = (raw_tilt / TILT_GRID_STEP).round() * TILT_GRID_STEP;
    let best_tilt = best_tilt.clamp(tilt_min, tilt_max);

    // Calculate final energy with optimal tilt
    let period_energy = calc_energy(best_tilt);

    PeriodOptimum {
        tilt_deg: (best_tilt * 10.0).round() / 10.0 + 0.0,
        period_energy_wh: period_energy,
        evaluations,
    }
}

/// Calculate period boundaries for the year
fn calculate_period_boundaries(days_in_year: u32, num_periods: u32) -> Vec<(u32, u32)> {
    assert!(num_periods >= 1);
    let mut periods = Vec::with_capacity(num_periods as usize);
    let period_length = days_in_year as f64 / num_periods as f64;

    for i in 0..num_periods {
        let start = (i as f64 * period_length).floor() as u32 + 1;
        let end = if i == num_periods - 1 {
            days_in_year
        } else {
            ((i + 1) as f64 * period_length).floor() as u32
        };
        periods.push((start, end));
    }

    periods
}

/// Select representative sample days for a period
/// For efficiency, we don't compute every day - we sample strategically
///
/// Samples are spread evenly across the period with both endpoints included,
/// so each one stands for the same number of days. Callers average them and
/// scale by the period length, which stepping by a fixed 3 or 7 days broke:
/// that left a short final gap before the appended end day, overweighting the
/// end of every period whose length was not a multiple of the step.
fn select_sample_days(start_day: u32, end_day: u32) -> Vec<u32> {
    let period_length = end_day - start_day + 1;

    // Short period: use all days
    if period_length <= 7 {
        return (start_day..=end_day).collect();
    }

    // Medium period: roughly every 3rd day. Long period: roughly weekly.
    let step = if period_length <= 31 { 3 } else { 7 };
    let count = period_length.div_ceil(step).max(2);

    let span = (period_length - 1) as f64;
    let steps = (count - 1) as f64;
    (0..count).map(|i| start_day + (i as f64 * span / steps).round() as u32).collect()
}

/// Check if a year is a leap year
fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

/// Generic Golden Section Search for finding the maximum of a unimodal function `f`
/// within the range `[min, max]`. Requires `min <= max` (callers normalize
/// their ranges; wrapped azimuth ranges are mapped to a continuous interval
/// before reaching this function).
///
/// Returns (x_at_max, max_value)
fn golden_section_search<F>(min: f64, max: f64, tol: f64, mut f: F) -> (f64, f64)
where
    F: FnMut(f64) -> f64,
{
    debug_assert!(min <= max, "golden_section_search needs an ordered range, got {min}..{max}");

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let resphi = 2.0 - phi; // approx 0.382

    let mut a = min;
    let mut b = max;

    // Initialize internal points
    let mut c = a + resphi * (b - a);
    let mut d = b - resphi * (b - a);

    let mut fc = f(c);
    let mut fd = f(d);

    while (b - a).abs() > tol {
        if fc < fd {
            // Max is in [c, b]
            a = c;
            c = d;
            fc = fd;
            d = b - resphi * (b - a);
            fd = f(d);
        } else {
            // Max is in [a, d]
            b = d;
            d = c;
            fd = fc;
            c = a + resphi * (b - a);
            fc = f(c);
        }
    }

    let best_x = (a + b) / 2.0;
    let best_y = fc.max(fd); // Approximate max value
    (best_x, best_y)
}

#[cfg(test)]
fn generate_sun_positions(
    lat: f64,
    lon: f64,
    altitude_m: f64,
    tz: chrono_tz::Tz,
    date: chrono::NaiveDate,
) -> Vec<(f64, f64, f64)> {
    use crate::solar::SolarCalc;
    use chrono::{Duration, TimeZone};

    let solar =
        SolarCalc { lat, lon, alt: altitude_m, delta_t: 69.0, refr: None, target: -1.064699 };
    // Use noon as a stable anchor.
    let noon = tz
        .from_local_datetime(&date.and_hms_opt(12, 0, 0).unwrap())
        .earliest()
        .expect("invalid local datetime in test helper");

    // Step size must match energy integration resolution
    let step_minutes = 10;
    let step = Duration::minutes(step_minutes);

    let mut positions = Vec::new();
    let start = noon - Duration::hours(12);
    let end = noon + Duration::hours(12);
    let mut t = start;

    while t <= end {
        let pos = solar.position(t).expect("SPA failed in test helper");

        // Elapsed hours since the start of the sampling window: strictly
        // monotonic and measured in real time, so it is immune to both the
        // hour-of-day midnight wrap and DST clock jumps which would distort
        // the trapezoidal time steps in calculate_daily_energy
        let hours = (t - start).num_seconds() as f64 / 3600.0;

        positions.push((hours, pos.elevation_angle(), pos.azimuth()));

        t += step;
    }

    positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solar_panel::SolarPanelConfig;
    use chrono::Datelike;
    use chrono::NaiveDate;
    use chrono_tz::Europe::Madrid;

    #[test]
    fn test_granada_summer_flat_panel_azimuth_regression() {
        // Granada, Spain
        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        // June 23, 2025 (matches your CLI example)
        let date = NaiveDate::from_ymd_opt(2025, 6, 23).unwrap();
        let day_of_year = date.ordinal();

        let sun_positions = generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

        assert!(sun_positions.len() > 20, "Sun position sampling unexpectedly small");

        let panel_config = SolarPanelConfig {
            area_m2: 1.0,
            tilt_deg: 35.0,
            azimuth_deg: 180.0,
            efficiency: 0.20,
            linke_turbidity: 3.0,
            albedo: 0.20,
            tracking_mode: crate::solar_panel::TrackingMode::Fixed,
        };

        let optimum = optimize_fixed_panel(panel_config, altitude_m, day_of_year, &sun_positions);

        // ---- Regression assertions ----

        // Near-horizontal is expected at summer solstice
        assert!(
            (optimum.azimuth_deg - 180.0).abs() <= AZ_GRID_STEP,
            "Azimuth should be near south for near-horizontal panels, got {:.1}°",
            optimum.azimuth_deg
        );

        // Azimuth MUST be canonicalized, never random
        assert_eq!(optimum.azimuth_deg, 180.0, "Azimuth must be canonicalized for flat panels");

        // Sanity check on energy magnitude
        assert!(optimum.energy_wh > 800.0, "Energy unexpectedly low: {:.1} Wh", optimum.energy_wh);
    }

    #[test]
    fn test_granada_winter_optimal_steep_and_south() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Madrid;

        // Granada, Spain
        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        // December 23, 2025 (winter solstice period)
        let date = NaiveDate::from_ymd_opt(2025, 12, 23).unwrap();
        let day_of_year = date.ordinal();

        let sun_positions = generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

        assert!(sun_positions.len() > 20, "Sun position sampling unexpectedly small");

        let panel_config = SolarPanelConfig {
            area_m2: 1.0,
            tilt_deg: 35.0,
            azimuth_deg: 180.0,
            efficiency: 0.20,
            linke_turbidity: 3.0,
            albedo: 0.20,
            tracking_mode: crate::solar_panel::TrackingMode::Fixed,
        };

        let optimum = optimize_fixed_panel(panel_config, altitude_m, day_of_year, &sun_positions);

        // ---- Winter-specific assertions ----

        // In winter at ~37°N, optimal tilt should be steep
        assert!(
            optimum.tilt_deg >= 50.0,
            "Expected steep winter optimal tilt, got {:.1}°",
            optimum.tilt_deg
        );

        // Azimuth must be meaningful and close to south
        assert!(
            (optimum.azimuth_deg - 180.0).abs() <= AZ_GRID_STEP.max(2.0),
            "Expected south-facing winter azimuth, got {:.1}°",
            optimum.azimuth_deg
        );

        // Must not be canonicalized-flat behavior
        assert!(optimum.tilt_deg > 5.0, "Winter tilt incorrectly treated as near-horizontal");

        // Sanity check on energy
        assert!(
            optimum.energy_wh > 300.0,
            "Winter energy unexpectedly low: {:.1} Wh",
            optimum.energy_wh
        );
    }

    #[test]
    fn test_granada_annual_average_tilt() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Madrid;

        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        // Representative days across the year
        let dates = [
            NaiveDate::from_ymd_opt(2025, 1, 15).unwrap(), // winter
            NaiveDate::from_ymd_opt(2025, 3, 20).unwrap(), // equinox
            NaiveDate::from_ymd_opt(2025, 6, 23).unwrap(), // summer
            NaiveDate::from_ymd_opt(2025, 9, 22).unwrap(), // equinox
            NaiveDate::from_ymd_opt(2025, 12, 23).unwrap(), // winter
        ];

        let mut tilts = Vec::new();

        for date in dates {
            let sun_positions =
                generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

            let panel_config = SolarPanelConfig {
                area_m2: 1.0,
                tilt_deg: 35.0,
                azimuth_deg: 180.0,
                efficiency: 0.20,
                linke_turbidity: 3.0,
                albedo: 0.20,
                tracking_mode: crate::solar_panel::TrackingMode::Fixed,
            };

            let optimum =
                optimize_fixed_panel(panel_config, altitude_m, date.ordinal(), &sun_positions);

            tilts.push(optimum.tilt_deg);
        }

        let avg_tilt: f64 = tilts.iter().sum::<f64>() / tilts.len() as f64;

        // ---- Annual-average assertion ----
        assert!(
            (32.0..=45.0).contains(&avg_tilt),
            "Expected annual-average tilt near latitude (+0\u{2026}+10°), got {:.1}°",
            avg_tilt
        );
    }

    #[test]
    fn test_granada_equinox_symmetry() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Madrid;

        // Granada, Spain
        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        // Spring equinox (approx)
        let date = NaiveDate::from_ymd_opt(2025, 3, 20).unwrap();
        let day_of_year = date.ordinal();

        let sun_positions = generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

        let panel_config = SolarPanelConfig {
            area_m2: 1.0,
            tilt_deg: 35.0,
            azimuth_deg: 180.0,
            efficiency: 0.20,
            linke_turbidity: 3.0,
            albedo: 0.20,
            tracking_mode: crate::solar_panel::TrackingMode::Fixed,
        };

        let optimum = optimize_fixed_panel(panel_config, altitude_m, day_of_year, &sun_positions);

        // ---- Equinox assertions ----

        // Tilt should be moderate (neither flat nor winter-steep)
        assert!(
            (20.0..=50.0).contains(&optimum.tilt_deg),
            "Expected moderate equinox tilt, got {:.1}°",
            optimum.tilt_deg
        );

        // Azimuth should be symmetric around south
        assert!(
            (optimum.azimuth_deg - 180.0).abs() <= 2.0,
            "Expected south-facing equinox azimuth, got {:.1}°",
            optimum.azimuth_deg
        );

        // Energy sanity check
        assert!(
            optimum.energy_wh > 600.0,
            "Equinox energy unexpectedly low: {:.1} Wh",
            optimum.energy_wh
        );
    }

    #[test]
    fn test_constrained_tilt_optimization() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Madrid;

        // Granada, Spain
        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        // Winter solstice - normally would have steep tilt
        let date = NaiveDate::from_ymd_opt(2025, 12, 23).unwrap();
        let day_of_year = date.ordinal();

        let sun_positions = generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

        let panel_config = SolarPanelConfig {
            area_m2: 1.0,
            tilt_deg: 35.0,
            azimuth_deg: 180.0,
            efficiency: 0.20,
            linke_turbidity: 3.0,
            albedo: 0.20,
            tracking_mode: crate::solar_panel::TrackingMode::Fixed,
        };

        // Constrain tilt to 20-40 degrees (roof limitations)
        let constraints = OptimizationConstraints::default().with_tilt_range(Some((20.0, 40.0)));

        let optimum = optimize_fixed_panel_constrained(
            panel_config,
            altitude_m,
            day_of_year,
            &sun_positions,
            &constraints,
        );

        // Tilt should be at constraint maximum (40°) since winter optimal is steeper
        assert!(
            optimum.tilt_deg >= 20.0 && optimum.tilt_deg <= 40.0,
            "Tilt should be within constraints, got {:.1}°",
            optimum.tilt_deg
        );

        // Should be close to max allowed since winter needs steeper angle
        assert!(
            optimum.tilt_deg >= 38.0,
            "Winter tilt should be near constraint max, got {:.1}°",
            optimum.tilt_deg
        );
    }

    #[test]
    fn test_constrained_azimuth_optimization() {
        use chrono::NaiveDate;
        use chrono_tz::Europe::Madrid;

        // Granada, Spain
        let latitude = 37.17;
        let longitude = -3.62;
        let altitude_m = 619.0;

        let date = NaiveDate::from_ymd_opt(2025, 6, 23).unwrap();
        let day_of_year = date.ordinal();

        let sun_positions = generate_sun_positions(latitude, longitude, altitude_m, Madrid, date);

        let panel_config = SolarPanelConfig {
            area_m2: 1.0,
            tilt_deg: 35.0,
            azimuth_deg: 180.0,
            efficiency: 0.20,
            linke_turbidity: 3.0,
            albedo: 0.20,
            tracking_mode: crate::solar_panel::TrackingMode::Fixed,
        };

        // Constrain azimuth to east-facing only (e.g., roof constraints)
        let constraints =
            OptimizationConstraints::default().with_azimuth_range(Some((60.0, 120.0)));

        let optimum = optimize_fixed_panel_constrained(
            panel_config,
            altitude_m,
            day_of_year,
            &sun_positions,
            &constraints,
        );

        // Azimuth should be within constraints
        assert!(
            optimum.azimuth_deg >= 60.0 && optimum.azimuth_deg <= 120.0,
            "Azimuth should be within constraints, got {:.1}°",
            optimum.azimuth_deg
        );
    }

    #[test]
    fn test_period_boundaries_calculation() {
        // Test equal division of year
        let periods = calculate_period_boundaries(365, 4);
        assert_eq!(periods.len(), 4);
        assert_eq!(periods[0].0, 1); // First period starts at day 1
        assert_eq!(periods[3].1, 365); // Last period ends at day 365

        // Check no gaps
        for i in 1..periods.len() {
            assert_eq!(
                periods[i].0,
                periods[i - 1].1 + 1,
                "Gap detected between periods {} and {}",
                i - 1,
                i
            );
        }

        // Test single period
        let single = calculate_period_boundaries(365, 1);
        assert_eq!(single.len(), 1);
        assert_eq!(single[0], (1, 365));

        // Test daily adjustments
        let daily = calculate_period_boundaries(365, 365);
        assert_eq!(daily.len(), 365);
        for (i, (start, end)) in daily.iter().enumerate() {
            assert_eq!(*start, (i as u32) + 1);
            assert_eq!(*end, (i as u32) + 1);
        }
    }

    #[test]
    fn test_sample_days_selection() {
        // Short period
        let short = select_sample_days(1, 5);
        assert_eq!(short.len(), 5);
        assert_eq!(short, vec![1, 2, 3, 4, 5]);

        // Medium period (2 weeks)
        let medium = select_sample_days(1, 14);
        assert!(medium.contains(&1)); // Start
        assert!(medium.contains(&14)); // End
        assert!(medium.len() >= 5 && medium.len() <= 10);

        // Long period (90 days)
        let long = select_sample_days(1, 90);
        assert!(long.contains(&1)); // Start
        assert!(long.contains(&90)); // End
        // Should have roughly weekly samples
        assert!(long.len() >= 10 && long.len() <= 20);
    }

    #[test]
    fn test_golden_section_search_basic() {
        // Test with a simple parabola: f(x) = -(x - 3)^2 + 10
        // This function has a global maximum at x = 3.0 with value 10.0
        // We search in the range [0.0, 6.0] which contains the maximum.

        let tolerance = 1e-6;
        let (x, y) = golden_section_search(0.0, 6.0, tolerance, |x| -(x - 3.0).powi(2) + 10.0);

        // Verify the location of the maximum
        assert!((x - 3.0).abs() < 1e-8, "Expected maximum at x=3.0, got x={:.5}", x);

        // Verify the value of the maximum
        assert!((y - 10.0).abs() < 1e-8, "Expected max value y=10.0, got y={:.5}", y);
    }

    // ------------------------------------------------------------------
    // Synthetic sun-position tests (no chrono / SPA dependency)
    // ------------------------------------------------------------------

    /// Synthetic clear day: the sun rises due east at 06:00, culminates at
    /// noon at `peak_elevation` in the given direction (180 = southern sky for
    /// northern-hemisphere sites, 0/360 = northern sky for southern-hemisphere
    /// sites) and sets due west at 18:00. 15-minute resolution.
    fn synthetic_day(peak_elevation: f64, culmination_az: f64) -> Vec<(f64, f64, f64)> {
        let southern_sky = (culmination_az - 180.0).abs() < 90.0;
        (0..=48)
            .map(|i| {
                let h = 6.0 + i as f64 * 0.25;
                let frac = (h - 6.0) / 12.0; // 0 at sunrise, 1 at sunset
                let elev = peak_elevation * (std::f64::consts::PI * frac).sin();
                // Sweep east -> culmination direction -> west
                let az = if southern_sky {
                    90.0 + 180.0 * frac
                } else {
                    (90.0 - 180.0 * frac).rem_euclid(360.0)
                };
                (h, elev, az)
            })
            .collect()
    }

    /// True if `az` lies within the (possibly north-wrapping) range
    fn azimuth_in_range(az: f64, min: f64, max: f64) -> bool {
        let az = az.rem_euclid(360.0);
        if min <= max { az >= min && az <= max } else { az >= min || az <= max }
    }

    fn synthetic_base() -> SolarPanelConfig {
        SolarPanelConfig::new(10.0, 0.0, 180.0)
            .with_efficiency(0.20)
            .with_linke_turbidity(3.0)
            .with_albedo(0.2)
    }

    #[test]
    fn test_inverted_tilt_range_does_not_panic() {
        let constraints = OptimizationConstraints::default().with_tilt_range(Some((40.0, 20.0)));
        let positions = synthetic_day(50.0, 180.0);
        let result =
            optimize_fixed_panel_constrained(synthetic_base(), 0.0, 172, &positions, &constraints);
        assert!(
            result.tilt_deg >= 20.0 && result.tilt_deg <= 40.0,
            "Tilt {} outside normalized range 20..40",
            result.tilt_deg
        );
        assert!(result.energy_wh > 0.0);
    }

    #[test]
    fn test_wrapped_azimuth_range() {
        let constraints =
            OptimizationConstraints::default().with_azimuth_range(Some((300.0, 60.0)));
        // Southern-hemisphere style day: sun culminates in the northern sky,
        // so the optimum azimuth is near 0/360 - inside the wrapped range
        let positions = synthetic_day(50.0, 0.0);
        let result =
            optimize_fixed_panel_constrained(synthetic_base(), 0.0, 355, &positions, &constraints);
        assert!(
            azimuth_in_range(result.azimuth_deg, 300.0, 60.0),
            "Azimuth {} outside wrapped range 300..60",
            result.azimuth_deg
        );
        assert!(result.tilt_deg > 5.0, "Expected a meaningful tilt, got {}", result.tilt_deg);
        assert!(result.energy_wh > 0.0);
    }

    #[test]
    fn test_explicit_full_circle_is_unconstrained() {
        // Passing the full circle explicitly must behave exactly like the
        // default: same wraparound-aware search, same result
        let full = OptimizationConstraints::default().with_azimuth_range(Some((0.0, 360.0)));
        assert!(full.azimuth_unconstrained());

        let positions = synthetic_day(45.0, 0.0);
        let explicit =
            optimize_fixed_panel_constrained(synthetic_base(), 0.0, 355, &positions, &full);
        let default = optimize_fixed_panel_constrained(
            synthetic_base(),
            0.0,
            355,
            &positions,
            &OptimizationConstraints::default(),
        );
        assert_eq!(explicit.tilt_deg, default.tilt_deg);
        assert_eq!(explicit.azimuth_deg, default.azimuth_deg);
    }

    #[test]
    fn test_southern_hemisphere_near_brute_force() {
        // The sun culminates in the northern sky.
        let positions = synthetic_day(40.0, 0.0); // winter-ish southern day
        let base = synthetic_base();

        let mut brute_best = 0.0_f64;
        let mut tilt = 0.0;
        while tilt <= 90.0 {
            let mut az = 0.0;
            while az < 360.0 {
                let cfg = SolarPanelConfig { tilt_deg: tilt, azimuth_deg: az, ..base };
                let e =
                    solar_panel::calculate_daily_energy(&cfg, 0.0, 355, positions.iter().copied());
                brute_best = brute_best.max(e);
                az += 5.0;
            }
            tilt += 2.5;
        }

        let result = optimize_fixed_panel_constrained(
            base,
            0.0,
            355,
            &positions,
            &OptimizationConstraints::default(),
        );

        // Optimum azimuth should be in the northern sky (near 0/360)
        assert!(
            azimuth_in_range(result.azimuth_deg, 315.0, 45.0),
            "Azimuth {} should face the sun's northern culmination",
            result.azimuth_deg
        );
        assert!(
            result.energy_wh >= brute_best * 0.98,
            "Optimizer found {:.0} Wh, brute force found {:.0} Wh",
            result.energy_wh,
            brute_best
        );
    }

    #[test]
    fn test_reported_energy_matches_reported_orientation() {
        // energy_wh must always be the energy at the returned (tilt, azimuth),
        // including after near-horizontal azimuth canonicalization
        let base = synthetic_base();
        for (peak, culm, day) in [(75.0, 180.0, 172), (40.0, 0.0, 355), (88.0, 180.0, 172)] {
            let positions = synthetic_day(peak, culm);
            let result = optimize_fixed_panel_constrained(
                base,
                0.0,
                day,
                &positions,
                &OptimizationConstraints::default(),
            );
            let cfg = SolarPanelConfig {
                tilt_deg: result.tilt_deg,
                azimuth_deg: result.azimuth_deg,
                ..base
            };
            let recomputed =
                solar_panel::calculate_daily_energy(&cfg, 0.0, day, positions.iter().copied());
            // Allow only rounding slack from the 0.1° display rounding
            assert!(
                (result.energy_wh - recomputed).abs() <= recomputed * 0.001 + 1e-6,
                "Reported {:.2} Wh != {:.2} Wh at reported orientation ({}, {})",
                result.energy_wh,
                recomputed,
                result.tilt_deg,
                result.azimuth_deg
            );
        }
    }

    #[test]
    fn test_hsat_tilt_optimizer_synthetic() {
        let positions = synthetic_day(35.0, 180.0); // low winter sun
        let base = SolarPanelConfig { azimuth_deg: 180.0, ..synthetic_base() };
        let result =
            optimize_hsat_tilt(base, 0.0, 355, &positions, &OptimizationConstraints::default());
        assert!(
            result.tilt_deg >= 0.0 && result.tilt_deg <= 90.0,
            "Axis tilt {} out of range",
            result.tilt_deg
        );
        // With a low sun, a leaning axis should beat a flat one
        let flat_cfg = SolarPanelConfig {
            tilt_deg: 0.0,
            tracking_mode: crate::solar_panel::TrackingMode::HorizontalAxis,
            ..base
        };
        let flat_energy =
            solar_panel::calculate_daily_energy(&flat_cfg, 0.0, 355, positions.iter().copied());
        assert!(result.tilt_deg > 0.0, "Expected a tilted axis for a low sun");
        assert!(result.energy_wh >= flat_energy, "Optimized axis tilt should not lose energy");
        assert_eq!(result.azimuth_deg, 180.0, "Lean azimuth should be reported");
    }

    /// Best energy found by scanning the allowed (tilt, azimuth) box.
    /// Coarser than the reporting grid to keep the test quick; callers allow
    /// for that in their tolerance.
    fn brute_force_best(
        base: SolarPanelConfig,
        day_of_year: u32,
        positions: &[(f64, f64, f64)],
        constraints: &OptimizationConstraints,
    ) -> f64 {
        let (az_lo, az_hi) = constraints.azimuth_interval();
        let mut best = 0.0_f64;
        let mut tilt = constraints.tilt_min;
        while tilt <= constraints.tilt_max {
            let mut az = az_lo;
            while az <= az_hi {
                let cfg =
                    SolarPanelConfig { tilt_deg: tilt, azimuth_deg: az.rem_euclid(360.0), ..base };
                best = best.max(solar_panel::calculate_daily_energy(
                    &cfg,
                    0.0,
                    day_of_year,
                    positions.iter().copied(),
                ));
                az += 2.0;
            }
            tilt += 1.0;
        }
        best
    }

    #[test]
    fn test_constrained_azimuth_not_trapped_by_flat_tilt() {
        // An azimuth range whose center is a direction the panel would rather
        // lie flat in must not suppress the azimuth search: the optimum here
        // is a tilted panel at the edge of the range, and settling for the
        // flat panel at the center threw away ~19% of the available energy
        let base = synthetic_base();
        let constraints =
            OptimizationConstraints::default().with_azimuth_range(Some((60.0, 120.0)));

        for (peak, day) in [(35.0, 355), (75.0, 172), (50.0, 100)] {
            let positions = synthetic_day(peak, 180.0);
            let result = optimize_fixed_panel_constrained(base, 0.0, day, &positions, &constraints);
            let best = brute_force_best(base, day, &positions, &constraints);

            assert!(
                azimuth_in_range(result.azimuth_deg, 60.0, 120.0),
                "Azimuth {} outside the requested range",
                result.azimuth_deg
            );
            assert!(
                result.energy_wh >= best * 0.995,
                "Peak {peak}°: optimizer found {:.0} Wh, grid search found {:.0} Wh",
                result.energy_wh,
                best
            );
        }
    }

    #[test]
    fn test_azimuth_constraint_membership_and_clamping() {
        let unconstrained = OptimizationConstraints::default();
        assert!(unconstrained.contains_azimuth(0.0));
        assert!(unconstrained.contains_azimuth(275.0));
        assert_eq!(unconstrained.clamp_azimuth(275.0), 275.0);

        let east = OptimizationConstraints::default().with_azimuth_range(Some((60.0, 120.0)));
        assert!(east.contains_azimuth(90.0));
        assert!(!east.contains_azimuth(180.0));
        // 180 is 60° from the 120 edge but 120° from the 60 edge
        assert_eq!(east.clamp_azimuth(180.0), 120.0);

        // Wrapping range: 0 and 350 are inside, 180 is not
        let wrapped = OptimizationConstraints::default().with_azimuth_range(Some((300.0, 60.0)));
        assert!(wrapped.contains_azimuth(0.0));
        assert!(wrapped.contains_azimuth(350.0));
        assert!(!wrapped.contains_azimuth(180.0));
        // Clamping always lands somewhere inside the range
        assert!(wrapped.contains_azimuth(wrapped.clamp_azimuth(180.0)));
        assert!(east.contains_azimuth(east.clamp_azimuth(300.0)));

        // A range ending at 360 contains north, whichever way north is written
        let north_edge =
            OptimizationConstraints::default().with_azimuth_range(Some((300.0, 360.0)));
        assert!(north_edge.contains_azimuth(0.0));
        assert!(north_edge.contains_azimuth(360.0));
    }

    /// Sun position generator over a synthetic year: the day length is fixed
    /// but the peak elevation and culmination direction follow the season
    fn synthetic_year_generator() -> SunPositionGenerator {
        Box::new(|lat: f64, _lon: f64, _alt: f64, day: u32| {
            let decl = 23.45_f64.to_radians()
                * (2.0 * std::f64::consts::PI * (284.0 + day as f64) / 365.0).sin();
            let decl_deg = decl.to_degrees();
            let peak = (90.0 - (lat - decl_deg).abs()).max(0.0);
            let culmination = if lat >= decl_deg { 180.0 } else { 0.0 };
            synthetic_day(peak, culmination)
        })
    }

    #[test]
    fn test_yearly_adjustments_never_lose_to_fixed() {
        // A single tilt is available to every period, so adjusting can never
        // produce less energy than not adjusting. Measuring the fixed baseline
        // on its own independent sampling of the year used to break this and
        // report negative gains.
        let generator = synthetic_year_generator();
        let constraints = OptimizationConstraints::default();

        for latitude in [37.17, 60.2, -33.9] {
            for num_adjustments in [1, 2, 3, 4, 6] {
                let result = optimize_yearly_adjustments(
                    synthetic_base(),
                    latitude,
                    0.0,
                    0.0,
                    2025,
                    num_adjustments,
                    &constraints,
                    &generator,
                );

                assert_eq!(result.periods.len(), num_adjustments as usize);
                assert!(
                    result.total_energy_wh >= result.fixed_optimal_energy_wh * 0.9999,
                    "lat {latitude}, {num_adjustments} adjustments: adjusted {:.0} Wh < fixed {:.0} Wh",
                    result.total_energy_wh,
                    result.fixed_optimal_energy_wh
                );

                // One adjustment IS the fixed case
                if num_adjustments == 1 {
                    assert!(
                        (result.total_energy_wh - result.fixed_optimal_energy_wh).abs()
                            <= result.fixed_optimal_energy_wh * 1e-9,
                        "Single period should equal the fixed optimum"
                    );
                }

                // Periods must tile the year without gaps or overlaps
                assert_eq!(result.periods[0].start_day, 1);
                assert_eq!(result.periods[result.periods.len() - 1].end_day, 365);
                for pair in result.periods.windows(2) {
                    assert_eq!(pair[1].start_day, pair[0].end_day + 1);
                }
            }
        }
    }

    #[test]
    fn test_yearly_adjustments_respect_azimuth_constraints() {
        // The equator-facing default must be moved into the allowed range
        // instead of silently reporting an azimuth the caller excluded
        let generator = synthetic_year_generator();
        let constraints =
            OptimizationConstraints::default().with_azimuth_range(Some((60.0, 120.0)));

        let result = optimize_yearly_adjustments(
            synthetic_base(),
            50.0,
            0.0,
            0.0,
            2025,
            4,
            &constraints,
            &generator,
        );

        for period in &result.periods {
            assert!(
                azimuth_in_range(period.azimuth_deg, 60.0, 120.0),
                "Period azimuth {} outside the requested 60..120 range",
                period.azimuth_deg
            );
        }
        // 120 is the allowed heading closest to south
        assert_eq!(result.periods[0].azimuth_deg, 120.0);
    }

    #[test]
    fn test_yearly_hsat_leap_year_and_bounds() {
        let generator = synthetic_year_generator();
        let result = optimize_yearly_hsat(
            synthetic_base(),
            60.2,
            0.0,
            0.0,
            2024, // leap year
            4,
            &OptimizationConstraints::default(),
            &generator,
        );

        assert_eq!(result.periods[result.periods.len() - 1].end_day, 366);
        assert!(result.total_energy_wh >= result.fixed_optimal_energy_wh * 0.9999);
        for period in &result.periods {
            // The HSAT axis is never searched past what the tracker accepts
            assert!(
                period.tilt_deg <= crate::solar_panel::MAX_TRACKER_AXIS_TILT_DEG,
                "Axis tilt {} exceeds the tracker limit",
                period.tilt_deg
            );
        }
    }

    #[test]
    fn test_hsat_axis_tilt_stays_within_tracker_limit() {
        // Asking for a near-vertical axis must report a tilt the tracker
        // model actually used, not the 90° that gets clamped internally
        let positions = synthetic_day(35.0, 180.0);
        let constraints = OptimizationConstraints::default().with_tilt_range(Some((89.5, 90.0)));
        let result = optimize_hsat_tilt(synthetic_base(), 0.0, 355, &positions, &constraints);
        assert!(
            result.tilt_deg <= crate::solar_panel::MAX_TRACKER_AXIS_TILT_DEG,
            "Reported axis tilt {} is past the tracker limit",
            result.tilt_deg
        );
    }

    #[test]
    fn test_sample_days_are_evenly_spread() {
        // Equal weighting in the period energy sum is only correct if each
        // sample stands for the same number of days
        for (start, end) in [(1u32, 8u32), (1, 14), (1, 31), (1, 32), (1, 90), (92, 182)] {
            let days = select_sample_days(start, end);
            assert_eq!(days[0], start, "Sampling must include the first day");
            assert_eq!(days[days.len() - 1], end, "Sampling must include the last day");

            let gaps: Vec<u32> = days.windows(2).map(|w| w[1] - w[0]).collect();
            let min_gap = *gaps.iter().min().unwrap();
            let max_gap = *gaps.iter().max().unwrap();
            assert!(min_gap >= 1, "Sample days must be strictly increasing");
            assert!(max_gap - min_gap <= 1, "Uneven sampling for {start}..{end}: gaps {gaps:?}");
        }
    }

    #[test]
    fn test_vsat_tilt_optimizer_synthetic() {
        let positions = synthetic_day(35.0, 180.0);
        let result = optimize_vsat_tilt(
            synthetic_base(),
            0.0,
            355,
            &positions,
            &OptimizationConstraints::default(),
        );
        // The sun's mean zenith over this day is large (peak elevation 35°),
        // so a substantial fixed tilt is expected
        assert!(
            result.tilt_deg > 30.0 && result.tilt_deg <= 90.0,
            "VSAT tilt {} unreasonable for a 35° peak sun",
            result.tilt_deg
        );
        assert!(result.energy_wh > 0.0);

        // And the optimizer must respect tilt constraints
        let constrained = optimize_vsat_tilt(
            synthetic_base(),
            0.0,
            355,
            &positions,
            &OptimizationConstraints::default().with_tilt_range(Some((10.0, 25.0))),
        );
        assert!(
            constrained.tilt_deg >= 10.0 && constrained.tilt_deg <= 25.0,
            "Constrained VSAT tilt {} outside 10..25",
            constrained.tilt_deg
        );
    }
}
