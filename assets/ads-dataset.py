# -*- coding: utf-8 -*-
import argparse
import xarray as xr
import numpy as np
import sys
import math


def format_day_ranges(days):
    if len(days) == 0:
        return "None"
    ranges = []
    start = prev = days[0]
    for d in days[1:]:
        if d != prev + 1:
            ranges.append(f"{start}-{prev}" if start != prev else str(start))
            start = d
        prev = d
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ", ".join(ranges)


def main():
    parser = argparse.ArgumentParser(
        description="Process ECMWF data into an optimized 3D binary blob."
    )
    parser.add_argument("files", nargs="+", help="NetCDF files to process")
    parser.add_argument(
        "-g",
        "--grid",
        type=float,
        default=5.0,
        help="Grid size in degrees (default: 5.0)",
    )
    parser.add_argument(
        "-d",
        "--days",
        type=int,
        default=1,
        help="Day granularity / chunk size (default: 1)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="atmosphere_grid.bin",
        help="Output binary file",
    )
    args = parser.parse_args()

    print(
        f"Loading files, target grid: {args.grid}°, granularity: {args.days} day(s)..."
    )
    try:
        ds = xr.open_mfdataset(args.files, combine="by_coords", join="outer")
    except Exception as e:
        print(f"Error loading files: {e}")
        sys.exit(1)

    print("Precalculating AOD_310 and Ozone DU...")
    # 1. Calculate Alpha, then immediately calculate AOD_310
    aod550 = ds["aod550"].clip(min=1e-5)
    aod865 = ds["aod865"].clip(min=1e-5)
    alpha = -np.log(aod550 / aod865) / np.log(550.0 / 865.0)

    # Precalculate AOD at 310nm. We will export this directly.
    ds["aod310"] = aod550 * (310.0 / 550.0) ** (-alpha)

    # 2. divisor for GTCO3 → DU is the column mass of 1 DU in kg/m²
    ds["ozone_du"] = ds["gtco3"] / 2.1415e-5

    print("Filtering and padding time dimension to 365 days...")
    # 3. Filter Leap Days
    ds = ds.sel(
        valid_time=~((ds.valid_time.dt.month == 2) & (ds.valid_time.dt.day == 29))
    )

    # 4. Collapse to 365 daily averages
    ds_daily = ds.groupby("valid_time.dayofyear").mean("valid_time")

    # --- NEW: Print coverage analysis for each dataset ---
    print("\n--- Data Coverage Analysis ---")
    for var_name in ["ozone_du", "aod550", "aluvd", "aluvp"]:
        if var_name in ds_daily:
            # Count how many valid (non-NaN) spatial pixels exist for each day
            valid_spatial_counts = (
                ds_daily[var_name].count(dim=["latitude", "longitude"]).values
            )
            # Keep only the days that have at least one valid pixel
            valid_days = ds_daily.dayofyear.values[valid_spatial_counts > 0]

            if len(valid_days) >= 365:
                print(f"  {var_name:>8}: Full 365-day coverage")
            else:
                print(
                    f"  {var_name:>8}: Partial ({len(valid_days)} days) -> Days {format_day_ranges(valid_days)}"
                )
    print("------------------------------\n")

    ds_daily = (
        ds_daily.reindex(dayofyear=np.arange(1, 366))
        .ffill(dim="dayofyear")
        .bfill(dim="dayofyear")
    )

    print(f"Chunking data into {args.days}-day bins...")
    # 5. Apply Day Granularity
    # This groups days (e.g., 1-7, 8-14) into single time bins
    ds_daily.coords["time_bin"] = (ds_daily.dayofyear - 1) // args.days
    ds_binned = ds_daily.groupby("time_bin").mean("dayofyear")

    print("Regridding spatial dimensions...")
    # 6. Regrid Longitude/Latitude
    ds_binned.coords["longitude"] = (ds_binned.coords["longitude"] + 180) % 360 - 180
    ds_binned = ds_binned.sortby(ds_binned.longitude)

    new_lat = np.arange(-90, 90, args.grid)
    new_lon = np.arange(-180, 180, args.grid)
    grid_target = ds_binned.interp(latitude=new_lat, longitude=new_lon, method="linear")

    # 7. Extract arrays with safe fallbacks
    ozone_arr = grid_target["ozone_du"].fillna(300.0).values.astype("<f4")
    aod310_arr = grid_target["aod310"].fillna(0.04).values.astype("<f4")

    # Check if the merged dataset contains Albedo
    if "aluvd" in grid_target:
        print("Found Diffuse UV Albedo (aluvd) in dataset.")
        albedo_arr = grid_target["aluvd"].fillna(1.02).values.astype("<f4")
    elif "aluvp" in grid_target:
        # Fallback just in case you only checked the direct box
        print("Found Direct UV Albedo (aluvp) in dataset.")
        albedo_arr = grid_target["aluvp"].fillna(1.02).values.astype("<f4")
    else:
        print("No Albedo data found. Falling back to default (1.02).")
        albedo_arr = np.full_like(ozone_arr, 1.02).astype("<f4")

    # 8. Stack into 3 variables (Ozone, AOD310, Albedo) -> 12 bytes per cell
    blob_3d = np.stack((ozone_arr, aod310_arr, albedo_arr), axis=-1)
    final_blob = np.ascontiguousarray(blob_3d)

    # 9. Build the expanded 16-byte Header
    time_steps = blob_3d.shape[0]
    header = np.array(
        [
            args.days,  # Granularity
            time_steps,  # Number of time chunks
            len(new_lat),  # Number of latitudes
            len(new_lon),  # Number of longitudes
        ],
        dtype="<u4",
    )

    with open(args.output, "wb") as f:
        f.write(header.tobytes())
        f.write(final_blob.tobytes())

    total_mb = (header.nbytes + final_blob.nbytes) / (1024 * 1024)
    print(f"\nSuccess! Packed into {total_mb:.2f} MiB.")
    print(
        f"Header Config: {args.days} Days/Chunk | {time_steps} Time Steps | {len(new_lat)} Lats | {len(new_lon)} Lons"
    )


if __name__ == "__main__":
    main()
