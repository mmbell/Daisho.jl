# API Reference

## Core Types

```@docs
Daisho.radar
```

## Radar I/O

```@docs
Daisho.read_cfradial
Daisho.get_radar_orientation
Daisho.write_qced_cfradial_sigmet
Daisho.write_qced_cfradial_singlepol
Daisho.write_qced_cfradial_dualpol
Daisho.write_qced_cfradial_P3
```

## Radar Utilities

```@docs
Daisho.initialize_moment_dictionaries
Daisho.initialize_qc_fields
Daisho.split_sweeps
Daisho.beam_height
Daisho.dB_to_linear
Daisho.linear_to_dB
Daisho.dB_to_linear!
Daisho.linear_to_dB!
```

## Quality Control

```@docs
Daisho.fix_SEAPOL_RHOHV!
Daisho.threshold_qc
Daisho.despeckle
Daisho.despeckle_azimuthal
Daisho.stddev_phidp_threshold
Daisho.remove_platform_motion!
Daisho.threshold_dbz
Daisho.threshold_height
Daisho.threshold_terrain_height
Daisho.add_azimuthal_offset
Daisho.mask_sector
Daisho.smooth_sqi
```

## Gridding

### Grid Initialization

```@docs
Daisho.initialize_regular_grid
```

### Radar Information

```@docs
Daisho.get_radar_zyx
Daisho.get_beam_info
Daisho.radar_arrays
Daisho.radar_balltree_yx
Daisho.radar_balltree_r
Daisho.appx_inverse_projection
```

### High-Level Gridding

```@docs
Daisho.grid_radar_volume
Daisho.grid_radar_latlon_volume
Daisho.grid_radar_rhi
Daisho.grid_radar_ppi
Daisho.grid_radar_composite
Daisho.grid_radar_column
```

### Core Gridding Algorithms

```@docs
Daisho.grid_volume
Daisho.grid_rhi
Daisho.grid_ppi
Daisho.grid_composite
Daisho.grid_column
```

### Gridded Data I/O

```@docs
Daisho.write_gridded_radar_volume
Daisho.write_gridded_radar_rhi
Daisho.write_gridded_radar_ppi
Daisho.write_gridded_radar_column
Daisho.read_gridded_radar
Daisho.read_gridded_ppi
Daisho.read_gridded_rhi
```

## SRTM / Terrain

```@docs
Daisho.terrain_height
Daisho.read_srtm_elevation_rasters
Daisho.read_srtm_elevation_multi
Daisho.read_srtm_elevation_dict
```
