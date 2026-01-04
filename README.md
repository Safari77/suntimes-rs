suntimes-rs ☀️
==============

A high-precision, Rust-based CLI tool for calculating accurate sun position, sunset and sunrise times, twilight phases, UV index, and estimating solar panel energy yields. Argos support to display info in GNOME Shell status bar.

suntimes-rs goes beyond simple sunrise/sunset times; it features a robust solar energy modeling engine capable of optimizing panel tilt and azimuth, supporting solar tracking systems (dual-axis, HSAT, VSAT), and generating seasonal adjustment schedules to maximize your energy production.

Conceptual model definitions
============================

### noaa

Purpose: match NOAA sunrise/sunset tables

Applies:
Upper limb (−0.266°)
Geometric horizon dip
Standard refraction only (1013.25 hPa, 15 °C)

Does NOT apply:
Pressure-from-altitude
Height-of-eye atmospheric correction
This is what NOAA effectively does.

### horizons

Purpose: match JPL Horizons Observer Table
Applies:
Observer altitude in SPA geometry
Standard refraction
No geometric dip correction

Does NOT apply:
Horizon dip
Pressure-from-altitude

### physical

Purpose: physically realistic observer model
Applies:
Observer altitude in geometry
Geometric horizon dip
Pressure estimated from altitude
Temperature from user

Does NOT apply:
Any NOAA/Horizons compatibility constraints

This is not directly comparable to Horizons — by design.

Argos setup
===========
Assuming you have already installed Argos.

```bash
$ cat > ~/.config/environment.d/argos.conf  << END
> ARGOS_SUNTIMES_ENABLE=true
ARGOS_SUNTIMES_LONGITUDE=5.5
ARGOS_SUNTIMES_LATITUDE=45.5
ARGOS_SUNTIMES_ALTITUDE=120
> END

$ ln -s ~/.cargo/bin/suntimes-rs ~/.config/argos/suntimes.5s.bin
```
And relogin.
You can also use bash script as the argos binary so changing parameters is easier:
```
#!/bin/sh
exec suntimes-rs --latitude ...
```

Solar panel power output offline estimation
===========================================
# Basic solar panel output
`suntimes-rs --latitude 60.17 --longitude 24.94 --solarpanel-size 10`

# Full configuration
```
suntimes-rs --latitude 60.17 --longitude 24.94 \
  --solarpanel-size 25 \
  --solarpanel-tilt 40 \
  --solarpanel-azimuth 170 \
  --solarpanel-efficiency 0.22 \
  --linke-turbidity 2.5 \
  --albedo 0.3
```

# Find optimal 4-season tilt adjustments
./suntimes-rs --latitude 60.2 --longitude 24.9 --solarpanel-size 10 \
  --solarpanel-yearly-adjustments 4

# With tilt constraints (roof limitations)
./suntimes-rs --latitude 60.2 --longitude 24.9 --solarpanel-size 10 \
  --solarpanel-yearly-adjustments 4 --solarpanel-tilt-range 25-55

# Simulate HSAT tracker output
./suntimes-rs --latitude 60.2 --longitude 24.9 --solarpanel-size 10 \
  --solarpanel-hsat

# Simulate VSAT tracker output
./suntimes-rs --latitude 60.2 --longitude 24.9 --solarpanel-size 10 \
  --solarpanel-vsat

# Compare fixed vs dual-axis tracking
./suntimes-rs --latitude 60.2 --longitude 24.9 --solarpanel-size 10 \
  --solarpanel-tracking

`--solarpanel-find-optimum`: Find optimum fixed tilt and azimuth for maximum daily energy

*Valid combinations*
| Options | Valid? | Notes |
| :--- | :---: | :--- |
| --solarpanel-horizontal-tracking | ✓ | HSAT |
| --solarpanel-vertical-tracking | ✓ | VSAT |
| --solarpanel-find-optimum | ✓ | Find optimal azimuth and tilt for fixed panel |
| --solarpanel-dual-axis | ✓ | Full tracking |
| --solarpanel-vertical-tracking --solarpanel-find-optimum | ✓ | Find optimal azimuth for VSAT |
| --solarpanel-horizontal-tracking --solarpanel-find-optimum | ✓ | Find optimal tilt for HSAT |
| --solarpanel-dual-axis --solarpanel-find-optimum | ✗ | Error (nothing to optimize) |
| --solarpanel-yearly-adjustments + any tracking | ✗ | Error (adjustments are for fixed panels) |

# Argos mode with solar panel
suntimes-rs --latitude 60.7 --longitude 22.94 --argos --solarpanel-size 10

Temporal validity and supported epochs
======================================
This program computes solar positions, sunrise/sunset, and twilight times using
high-precision astronomical models. Its numerical validity depends on several
underlying models, each with different applicability ranges. The effective
supported epoch is the intersection of these models.

### Supported calendar range
Supported years: −500 to +5000 (astronomical year numbering)
Location to timezone feature (--timezone location) is reliable only from year 1970 onwards.
Calendar system: proleptic Gregorian
Time representation: UTC or modern IANA time zones
Dates outside this range are rejected because required physical models (ΔT) are not available.

### Model components and limits
Solar geometry (SPA)
 - Uses the NREL Solar Position Algorithm (SPA)
 - Valid approximately from −2000 to +6000
 - Provides sub-arcminute accuracy for solar azimuth and altitude

This is not the limiting factor.

ΔT (Terrestrial Time − Universal Time)
ΔT is required to convert between Earth rotation time and ephemeris time
ΔT is estimated, not known exactly, outside the modern era
The ΔT model used by this program is defined for years ≥ −500
For earlier years, ΔT uncertainty becomes too large and the model refuses to extrapolate.

Notes
=====
`--timezone location` uses about 30 ms on startup.
tzf-rs loads a binary map of the entire world’s timezone boundaries. This map is stored as a Protocol Buffer.
