using Rasters, ArchGDAL, Printf

"""
    read_srtm_elevation_rasters(filename::String, lat::Float64, lon::Float64) -> Union{Float64, Nothing}

Read a NASA SRTM `.hgt` file and extract the elevation at a specified latitude and longitude
using Rasters.jl.

The SRTM file is loaded as a raster in WGS84 geographic coordinates, and the nearest-neighbor
elevation value is returned. Void data values (`-32768` or `missing`) return `nothing`.

# Arguments
- `filename::String`: Path to an SRTM `.hgt` file.
- `lat::Float64`: Latitude in decimal degrees (positive north).
- `lon::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters at the specified coordinates.
- `nothing`: If the elevation data is missing or void at the location.

# Throws
- `ErrorException`: If the coordinates are outside the raster bounds.
"""
function read_srtm_elevation_rasters(filename::String, lat::Float64, lon::Float64)
    try
        # Read the SRTM file as a raster
        # SRTM files are in geographic coordinates (WGS84)
        dem = Raster(filename)
        
        # Extract elevation at the specified coordinates
        # Rasters.jl handles the coordinate transformation automatically
        elevation = dem[X(Near(lon)),Y(Near(lat))]
        
        # Handle missing/void data (typically represented as missing or specific fill values)
        if ismissing(elevation) || elevation == -32768
            return nothing
        else
            return Float64(elevation)
        end
        
    catch e
        if isa(e, BoundsError)
            error("Coordinates ($lat, $lon) are outside the raster bounds")
        else
            rethrow(e)
        end
    end
end

"""
    read_srtm_elevation_multi(tile_directory::String, lat::Float64, lon::Float64) -> Union{Float64, Nothing}

Look up terrain elevation from SRTM tiles stored in a directory for a given latitude and longitude.

The appropriate SRTM tile filename is automatically determined from the coordinates using standard
naming conventions (e.g., `N16W025.hgt`). If the tile file does not exist, ocean is assumed and
`nothing` is returned.

# Arguments
- `tile_directory::String`: Path to directory containing SRTM `.hgt` tile files.
- `lat::Float64`: Latitude in decimal degrees (positive north).
- `lon::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters at the specified coordinates.
- `nothing`: If the tile is missing (assumed ocean) or elevation data is void.
"""
function read_srtm_elevation_multi(tile_directory::String, lat::Float64, lon::Float64)
    # Generate the expected tile name based on coordinates
    lat_char = lat >= 0 ? 'N' : 'S'
    lon_char = lon >= 0 ? 'E' : 'W'
    
    lat_deg = Int(floor(abs(lat)))
    lon_deg = Int(floor(abs(lon)))
    
    # Adjust for tile naming convention
    if lat < 0
        lat_deg = lat_deg + 1
    end
    if lon < 0  
        lon_deg = lon_deg + 1
    end
    
    tile_name = string(lat_char, lpad(lat_deg, 2, '0'), lon_char, lpad(lon_deg, 3, '0'), ".hgt")
    tile_path = joinpath(tile_directory, tile_name)
    
    if !isfile(tile_path)
        #error("SRTM tile not found: $tile_path")
        # Assume it is ocean if tile is missing
        return nothing
    end
    
    return read_srtm_elevation_rasters(tile_path, lat, lon)
end

"""
    terrain_height(tile_directory::String, latitude::Float64, longitude::Float64) -> Float64

Look up terrain elevation at a single geographic point from SRTM tiles stored on disk.

Loads the appropriate SRTM tile from the specified directory and returns the elevation.
Returns `-1.0` if no elevation data is available (e.g., over ocean) or if an error occurs.

# Arguments
- `tile_directory::String`: Path to directory containing SRTM `.hgt` tile files.
- `latitude::Float64`: Latitude in decimal degrees (positive north).
- `longitude::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters, or `-1.0` if data is unavailable or an error occurs.
"""
function terrain_height(tile_directory::String, latitude::Float64, longitude::Float64)

    try
        elevation = read_srtm_elevation_multi(tile_directory, latitude, longitude)
        if elevation === nothing
            #println("No elevation data available at ($latitude, $longitude)")
            return -1.0
        else
            #println("Elevation at ($latitude, $longitude): $(elevation) meters")
            return elevation
        end
    catch e
        println("Error: $e")
        return -1.0
    end
end

"""
    terrain_height(tiles::Dict{String, Raster}, latitude::Float64, longitude::Float64) -> Float64

Look up terrain elevation at a single geographic point using a pre-loaded dictionary of SRTM
raster tiles.

This method avoids repeated disk I/O by using tiles already loaded into memory via
[`read_srtm_elevation_multi(tile_directory)`](@ref). Returns `-1.0` if no elevation data is
available or if an error occurs.

# Arguments
- `tiles::Dict{String, Raster}`: Dictionary mapping tile name strings (e.g., `"N16W025"`) to loaded `Raster` objects.
- `latitude::Float64`: Latitude in decimal degrees (positive north).
- `longitude::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters, or `-1.0` if data is unavailable or an error occurs.
"""
function terrain_height(tiles::Dict{String, Raster}, latitude::Float64, longitude::Float64)

    try
        elevation = read_srtm_elevation_dict(tiles, latitude, longitude)
        if elevation === nothing
            #println("No elevation data available at ($latitude, $longitude)")
            return -1.0
        else
            #println("Elevation at ($latitude, $longitude): $(elevation) meters")
            return elevation
        end
    catch e
        println("Error: $e")
        return -1.0
    end
end

"""
    read_srtm_elevation_multi(tile_directory::String) -> Dict{String, Raster}

Load all SRTM `.hgt` tile files from a directory into memory as a dictionary of `Raster` objects.

Each tile is keyed by its base name without extension (e.g., `"N16W025"`). The returned
dictionary can be passed to [`terrain_height(tiles, lat, lon)`](@ref) for efficient repeated
lookups without disk I/O.

# Arguments
- `tile_directory::String`: Path to directory containing SRTM `.hgt` tile files.

# Returns
- `Dict{String, Raster}`: Dictionary mapping tile name strings to loaded `Raster` objects.
"""
function read_srtm_elevation_multi(tile_directory::String)

    # Load in all the tiles in the directory
    tile_files = filter(f -> endswith(f, ".hgt"), readdir(tile_directory, join=true))
    tiles = Dict{String, Raster}()
    for tile_file in tile_files
        tile_name = split(basename(tile_file), ".")[1]
        tiles[tile_name] = Raster(tile_file)
    end
    return tiles
end

"""
    read_srtm_elevation_rasters(raster_dem::Raster, lat::Float64, lon::Float64) -> Union{Float64, Nothing}

Extract the elevation at a specified latitude and longitude from an already-loaded SRTM `Raster` object.

This method avoids re-reading the file from disk by accepting a pre-loaded `Raster`. The
nearest-neighbor elevation value is returned. Void data values (`-32768` or `missing`) return `nothing`.

# Arguments
- `raster_dem::Raster`: A pre-loaded SRTM raster in WGS84 geographic coordinates.
- `lat::Float64`: Latitude in decimal degrees (positive north).
- `lon::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters at the specified coordinates.
- `nothing`: If the elevation data is missing or void at the location.

# Throws
- `ErrorException`: If the coordinates are outside the raster bounds.
"""
function read_srtm_elevation_rasters(raster_dem::Raster, lat::Float64, lon::Float64)

    try
        # SRTM files are in geographic coordinates (WGS84)        
        # Extract elevation at the specified coordinates
        # Rasters.jl handles the coordinate transformation automatically
        elevation = raster_dem[X(Near(lon)),Y(Near(lat))]
        
        # Handle missing/void data (typically represented as missing or specific fill values)
        if ismissing(elevation) || elevation == -32768
            return nothing
        else
            return Float64(elevation)
        end
        
    catch e
        if isa(e, BoundsError)
            error("Coordinates ($lat, $lon) are outside the raster bounds")
        else
            rethrow(e)
        end
    end
end

"""
    read_srtm_elevation_dict(tiles::Dict{String, Raster}, lat::Float64, lon::Float64) -> Union{Float64, Nothing}

Look up terrain elevation from a pre-loaded dictionary of SRTM raster tiles.

The appropriate tile name is determined from the coordinates using standard SRTM naming
conventions. If the tile is not present in the dictionary, ocean is assumed and `nothing` is
returned.

# Arguments
- `tiles::Dict{String, Raster}`: Dictionary mapping tile name strings (e.g., `"N16W025"`) to loaded `Raster` objects.
- `lat::Float64`: Latitude in decimal degrees (positive north).
- `lon::Float64`: Longitude in decimal degrees (positive east).

# Returns
- `Float64`: Elevation in meters at the specified coordinates.
- `nothing`: If the tile is not in the dictionary (assumed ocean) or elevation data is void.
"""
function read_srtm_elevation_dict(tiles::Dict{String, Raster}, lat::Float64, lon::Float64)
    # Generate the expected tile name based on coordinates
    lat_char = lat >= 0 ? 'N' : 'S'
    lon_char = lon >= 0 ? 'E' : 'W'
    
    lat_deg = Int(floor(abs(lat)))
    lon_deg = Int(floor(abs(lon)))
    
    # Adjust for tile naming convention
    if lat < 0
        lat_deg = lat_deg + 1
    end
    if lon < 0  
        lon_deg = lon_deg + 1
    end
    
    tile_name = string(lat_char, lpad(lat_deg, 2, '0'), lon_char, lpad(lon_deg, 3, '0'))
    
    if !haskey(tiles, tile_name)
        #error("SRTM tile not found: $tile_path")
        # Assume it is ocean if tile is missing
        return nothing
    end

    return read_srtm_elevation_rasters(tiles[tile_name], lat, lon)
end
# Run examples
#tile_dir = "/Users/mmbell/Science/PICCOLO/srtm/SRTM3"
#terrain_height(tile_dir, 16.886, -24.988)
# example_beam_blockage()