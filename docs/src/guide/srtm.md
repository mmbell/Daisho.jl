# SRTM Terrain Integration

Daisho.jl integrates NASA Shuttle Radar Topography Mission (SRTM) elevation data for terrain-aware quality control.

## Overview

SRTM data is used to identify radar gates that are contaminated by ground clutter. By comparing the height of each radar gate with the terrain elevation at the same horizontal location, gates that intersect the ground can be removed.

## Tile Format

SRTM data is distributed as `.hgt` files, named by their southwest corner coordinates:

- `N16W025.hgt` - tile covering 16-17°N, 24-25°W
- `S34W071.hgt` - tile covering 33-34°S, 70-71°W

## Usage

### Single Point Lookup

```julia
# From a directory of .hgt files
elevation = Daisho.terrain_height("/path/to/srtm/tiles", 16.886, -24.988)

# Returns elevation in meters, or -1.0 if no data available (e.g., ocean)
```

### Pre-loading Tiles for Performance

For repeated lookups (e.g., during QC), pre-load all tiles into memory:

```julia
# Load all tiles from a directory
tiles = Daisho.read_srtm_elevation_multi("/path/to/srtm/tiles")

# Fast lookup using pre-loaded tiles
elevation = Daisho.terrain_height(tiles, 16.886, -24.988)
```

### Terrain-Based QC

Apply terrain masking to radar data:

```julia
qc_moments = Daisho.threshold_terrain_height(
    volume, raw_dict, qc_moments, qc_dict,
    5.0,                      # terrain threshold (meters)
    "/path/to/srtm/tiles"     # tile directory
)
```

Gates where the terrain height exceeds the threshold are set to missing. A threshold of 0 meters masks all land; higher values allow coastal and low-elevation returns.

## Dependencies

SRTM reading uses [Rasters.jl](https://github.com/rafaqz/Rasters.jl) and [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl) for geospatial raster I/O.
