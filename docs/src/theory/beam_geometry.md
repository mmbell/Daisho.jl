# Beam Geometry

## Beam Height Calculation

Daisho.jl calculates radar beam height using the standard 4/3 effective Earth radius model, which accounts for atmospheric refraction under standard conditions.

### Effective Earth Radius

The effective Earth radius for standard refraction is:

```math
R_{eff} = \frac{4}{3} R_E = \frac{4}{3} \times 6371000 \text{ m} \approx 8495333 \text{ m}
```

### Beam Height Formula

The height of the radar beam center at slant range ``r`` and elevation angle ``\theta_e`` is:

```math
h = \sqrt{r^2 + R_{eff}^2 + 2 r R_{eff} \sin(\theta_e)} - R_{eff} + h_{radar}
```

where ``h_{radar}`` is the radar altitude above sea level.

This is implemented in [`Daisho.beam_height`](@ref).

### Reference

Doviak, R. J., and D. S. Zrnic, 1993: *Doppler Radar and Weather Observations*, 2nd ed., Academic Press.

## Surface Range

The surface (ground) range corresponding to a beam at slant range ``r`` and height ``h`` is:

```math
s = R_{eff} \arcsin\left(\frac{r \cos(\theta_e)}{R_{eff} + h}\right)
```

## Coordinate System

Daisho uses a Transverse Mercator projection centered on the radar location (or a reference point). Coordinates are:

- **Z**: Height above sea level (meters)
- **Y**: Northward distance from the projection origin (meters)
- **X**: Eastward distance from the projection origin (meters)

The projection is computed using CoordRefSystems.jl for accuracy.

## Approximate Inverse Projection

For converting projected (Y, X) coordinates back to (latitude, longitude), Daisho uses an approximate formula based on FCC wireless communication specifications:

```math
\Delta\text{lat} = \frac{Y}{f_{lat} \times 1000}
```
```math
\Delta\text{lon} = \frac{X}{f_{lon} \times 1000}
```

where the latitude and longitude factors are:

```math
f_{lat} = 111.13209 - 0.56605 \cos(2\phi) + 0.00012 \cos(4\phi) - 0.000002 \cos(6\phi)
```
```math
f_{lon} = 111.41513 \cos(\phi) - 0.09455 \cos(3\phi) + 0.00012 \cos(5\phi)
```

This is implemented in [`Daisho.appx_inverse_projection`](@ref).
