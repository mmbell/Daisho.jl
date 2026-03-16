# Gridding Algorithm

## Overview

Daisho.jl implements a beam-aware radar gridding algorithm that maps irregular radar observations onto regular grids. The algorithm accounts for:

- Earth curvature and standard atmospheric refraction
- Radar beam pattern (Gaussian approximation)
- Range weighting
- Beam broadening with distance (beam inflation)

## Algorithm Steps

### 1. Spatial Indexing

A BallTree (from NearestNeighbors.jl) is constructed from the projected (Y, X) locations of all radar gates. This enables efficient radius queries during gridding.

### 2. For Each Grid Point

For each grid point ``(x_g, y_g, z_g)``:

1. **Find nearby gates**: Query the BallTree for all gates within the effective horizontal radius of influence
2. **Check vertical coverage**: Verify at least one gate is within the vertical radius of influence
3. **For each nearby valid gate**:
   a. Compute the effective elevation angle from the gate to the grid point
   b. Compute the effective azimuth angle
   c. Calculate the angular separation (spherical angle)
   d. Apply beam-pattern weighting
   e. Apply range weighting
   f. Accumulate weighted contributions

### 3. Weight Computation

#### Beam Pattern Weight

The beam pattern is modeled as a Gaussian with half-power beamwidth of 1 degree:

```math
w_{beam} = \exp(-\Delta\theta \times 79.43)
```

where ``\Delta\theta`` is the spherical angle between the beam direction and the grid point direction. The constant 79.43 corresponds to ``-\ln(0.5) / (0.5° \times \pi/180)``.

Gates with ``w_{beam}`` below the power threshold are excluded.

#### Range Weight

```math
w_{range} = \frac{r_{gridpt}}{r_{gate}}
```

where ``r_{gridpt}`` is the estimated range from the radar to the grid point, and ``r_{gate}`` is the actual range of the gate. Gates where the range difference exceeds the radius of influence are excluded.

#### Total Weight

```math
w_{total} = w_{beam} \times w_{range}
```

### 4. Interpolation Types

The interpolation method depends on the moment type:

| Grid Type | Method | Use Case |
|-----------|--------|----------|
| `:linear` | Convert to linear, weighted average, convert back to dB | Reflectivity (dBZ) |
| `:weighted` | Direct weighted average | Velocity, spectrum width |
| `:nearest` | Keep value with highest weight | Nearest-neighbor |

### 5. Beam Inflation

The effective radius of influence can expand proportionally to the distance from the radar:

```math
R_{eff} = \max(\alpha \times d, R_{base})
```

where ``\alpha`` is the `beam_inflation` parameter, ``d`` is the distance from the radar origin, and ``R_{base}`` is the base radius of influence.

This accounts for the broadening of the radar beam with range.

## Grid Flag Values

- **-32768.0**: No radar coverage (no gates within radius of influence)
- **-9999.0**: Radar coverage exists but all data has been QC'd out
- **Other values**: Valid gridded data

## Composite Gridding

The composite mode (`grid_composite`) selects the maximum reflectivity value from all gates within the radius of influence, rather than performing weighted interpolation. All other moments are taken from the same gate as the maximum reflectivity.

## Column Gridding

The column mode (`grid_column`) produces a single vertical profile above the radar location. Only vertical (elevation angle) weighting is applied.

## Parallelization

All gridding functions use `Threads.@threads` for parallelization across horizontal grid points. Set the `JULIA_NUM_THREADS` environment variable to control the number of threads.
