using Rasters, ArchGDAL, Printf

"""
Read NASA SRTM .hgt file and get elevation at specified latitude/longitude using Rasters.jl
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
Read multiple SRTM tiles and get elevation for coordinates that might span tiles
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

# Single point elevation lookup example
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
    end
end

# Get elevation using pre-loaded tiles dictionary
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
    end
end

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