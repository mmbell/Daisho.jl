# Quality Control

Daisho.jl provides a suite of quality control functions for radar data. The QC workflow typically proceeds in several stages.

## QC Workflow

### 1. Initialize QC Fields

Copy raw moments into QC arrays that will be progressively filtered:

```julia
qc_moments = Daisho.initialize_qc_fields(volume, raw_dict, qc_dict)
```

### 2. Apply Threshold Filters

Remove gates based on signal quality indicators:

```julia
# Remove gates with low SQI
qc_moments = Daisho.threshold_qc(volume.moments, raw_dict, qc_moments, qc_dict,
    "SQI", 0.3, "SQI", true)

# Remove gates with high spectrum width
qc_moments = Daisho.threshold_qc(volume.moments, raw_dict, qc_moments, qc_dict,
    "WIDTH", 8.0, "WIDTH", false)
```

Or apply multiple thresholds at once:

```julia
qc_moments = Daisho.threshold_qc(volume.moments, raw_dict, qc_moments, qc_dict,
    0.35, 6.0, 0.7, 8.0)  # SQI, SNR, RHOHV, WIDTH thresholds
```

### 3. Remove Speckle

Remove isolated gates (speckle) in range and azimuth:

```julia
qc_moments = Daisho.despeckle(3, qc_moments, qc_dict, n_gates, n_rays)
qc_moments = Daisho.despeckle_azimuthal(3, qc_moments, qc_dict, n_gates, n_rays)
```

### 4. PHIDP-Based Filtering

Remove gates with high differential phase variability:

```julia
qc_moments = Daisho.stddev_phidp_threshold(qc_moments, qc_dict, n_gates, n_rays, 11, 12)
```

### 5. Platform Motion Removal

For moving platforms (ships, aircraft), remove platform motion from radial velocity:

```julia
Daisho.remove_platform_motion!(volume, qc_moments, qc_dict)
```

### 6. Terrain Masking

Remove gates contaminated by ground clutter using SRTM elevation data:

```julia
qc_moments = Daisho.threshold_terrain_height(volume, raw_dict, qc_moments, qc_dict,
    5.0, "/path/to/srtm/tiles")
```

### 7. Sector Masking

Mask data in specific heading-relative azimuth sectors (e.g., ship superstructure blockage):

```julia
heading = get_heading_array(volume)
qc_moments = Daisho.mask_sector(volume, raw_dict, qc_moments, qc_dict,
    heading, 170.0, 190.0, "SQI")
```

## SEA-POL Specific

For the SEA-POL C-band radar, RHOHV requires special handling:

```julia
Daisho.fix_SEAPOL_RHOHV!(volume, raw_dict)
```

## Utility Functions

- [`Daisho.add_azimuthal_offset`](@ref) - Apply constant azimuth offset (with wrapping)
- [`Daisho.threshold_height`](@ref) - Remove gates below a height threshold
- [`Daisho.threshold_dbz`](@ref) - Advanced clutter detection based on reflectivity, velocity, and spectrum width
