$ python3 ads-dataset.py -g 4 -d 4 /tmp/data_sfc_*
Loading files, target grid: 4.0°, granularity: 4 day(s)...
Precalculating AOD_310 and Ozone DU...
Filtering and padding time dimension to 365 days...

--- Data Coverage Analysis ---
  ozone_du: Full 365-day coverage
    aod550: Full 365-day coverage
     aluvd: Full 365-day coverage
     aluvp: Full 365-day coverage
------------------------------

Chunking data into 4-day bins...
Regridding spatial dimensions...
Found Diffuse UV Albedo (aluvd) in dataset.

Success! Packed into 4.26 MiB.
Header Config: 4 Days/Chunk | 92 Time Steps | 45 Lats | 90 Lons
