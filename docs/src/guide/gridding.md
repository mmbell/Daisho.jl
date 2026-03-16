# Gridding

Daisho.jl implements a beam-aware radar gridding algorithm that accounts for Earth curvature, standard atmospheric refraction, and the radar beam pattern.

## Grid Types

Several gridding modes are available:

| Mode | Function | Dimensions | Use Case |
|------|----------|------------|----------|
| Volume | `grid_radar_volume` | 3D (X, Y, Z) | Full volumetric analysis |
| Lat/Lon Volume | `grid_radar_latlon_volume` | 3D (lon, lat, Z) | Geographic coordinate grids |
| RHI | `grid_radar_rhi` | 2D (R, Z) | Range-height cross sections |
| PPI | `grid_radar_ppi` | 2D (X, Y) | Plan position indicator |
| Composite | `grid_radar_composite` | 2D (X, Y) | Maximum value composite |
| Column | `grid_radar_column` | 1D (Z) | Vertical profile above radar |

## Basic Usage

### Volume Gridding

```julia
Daisho.grid_radar_volume(
    volume, qc_dict, grid_type_dict, "output.nc", index_time,
    -50000.0, 500.0, 201,    # X: min, increment, dimension
    -50000.0, 500.0, 201,    # Y: min, increment, dimension
    0.0, 500.0, 37,          # Z: min, increment, dimension
    0.01,                     # beam_inflation
    0.5                       # power_threshold
)
```

### PPI Gridding

```julia
Daisho.grid_radar_ppi(
    volume, qc_dict, grid_type_dict, "output_ppi.nc", index_time,
    -125000.0, 500.0, 501,   # X
    -125000.0, 500.0, 501,   # Y
    0.01, 0.5                # beam_inflation, power_threshold
)
```

## Gridding Algorithm

The algorithm uses a BallTree spatial index for efficient nearest-neighbor queries:

1. Build a BallTree of all radar gate locations in projected coordinates
2. For each grid point, find nearby radar gates within the radius of influence
3. Compute beam-pattern weights (Gaussian based on angular separation)
4. Compute range weights
5. Apply weighted interpolation (linear for dBZ, weighted average for velocity, etc.)

## Key Parameters

- **`beam_inflation`**: Expands the radius of influence proportionally to distance from the radar. Accounts for beam broadening with range. Typical value: 0.01.
- **`power_threshold`**: Minimum beam pattern weight to include a gate. Gates with weight below this threshold are excluded. Typical value: 0.5.
- **`missing_key`**: The moment used to determine if a gate has valid data (e.g., "SQI").
- **`valid_key`**: The moment used to check for non-missing data (e.g., "DBZ").

## Coordinate Systems

Daisho uses Transverse Mercator projections (via CoordRefSystems.jl) centered on the radar location. The approximate inverse projection function converts projected coordinates back to lat/lon:

```julia
lat, lon = Daisho.appx_inverse_projection(ref_lat, ref_lon, [y_meters, x_meters])
```

## Output Format

Gridded data is written as CF-compliant NetCDF files with:
- Coordinate variables (X, Y, Z or lat, lon)
- Grid mapping metadata (Transverse Mercator)
- All gridded moments as 3D or 2D variables
- Time, start_time, stop_time
- Latitude/longitude grids for each horizontal point
