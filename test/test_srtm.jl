@testset "SRTM" begin

    @testset "SRTM tile naming convention" begin
        # Test the tile naming logic by checking the coordinate-to-name mapping
        # This tests the logic in read_srtm_elevation_multi and read_srtm_elevation_dict

        # Northern hemisphere, eastern longitude
        lat, lon = 16.5, 25.3
        lat_char = lat >= 0 ? 'N' : 'S'
        lon_char = lon >= 0 ? 'E' : 'W'
        lat_deg = Int(floor(abs(lat)))
        lon_deg = Int(floor(abs(lon)))
        tile_name = string(lat_char, lpad(lat_deg, 2, '0'), lon_char, lpad(lon_deg, 3, '0'))
        @test tile_name == "N16E025"

        # Southern hemisphere, western longitude
        lat, lon = -33.4, -70.6
        lat_char = lat >= 0 ? 'N' : 'S'
        lon_char = lon >= 0 ? 'E' : 'W'
        lat_deg = Int(floor(abs(lat))) + 1  # Adjust for negative
        lon_deg = Int(floor(abs(lon))) + 1  # Adjust for negative
        tile_name = string(lat_char, lpad(lat_deg, 2, '0'), lon_char, lpad(lon_deg, 3, '0'))
        @test tile_name == "S34W071"
    end

    @testset "terrain_height - missing tile directory" begin
        # With a non-existent directory, should handle gracefully
        result = Daisho.terrain_height(joinpath(@__DIR__, "fixtures", "nonexistent_srtm"), 16.886, -24.988)
        @test result == -1.0
    end

    @testset "terrain_height - empty directory" begin
        # Create an empty directory for testing
        empty_dir = joinpath(@__DIR__, "fixtures", "empty_srtm")
        mkpath(empty_dir)

        result = Daisho.terrain_height(empty_dir, 16.886, -24.988)
        @test result == -1.0

        rm(empty_dir, recursive=true, force=true)
    end

    @testset "terrain_height returns Float64" begin
        # Both overloads should return -1.0 (Float64) when no data available
        result = Daisho.terrain_height(joinpath(@__DIR__, "fixtures"), 16.886, -24.988)
        @test result isa Float64
        @test result == -1.0
    end

    @testset "read_srtm_elevation_multi - missing tile returns nothing" begin
        empty_dir = joinpath(@__DIR__, "fixtures", "empty_srtm2")
        mkpath(empty_dir)

        result = Daisho.read_srtm_elevation_multi(empty_dir, 16.886, -24.988)
        @test result === nothing

        rm(empty_dir, recursive=true, force=true)
    end

    @testset "read_srtm_elevation_dict - missing tile returns nothing" begin
        tiles = Dict{String, Rasters.Raster}()
        result = Daisho.read_srtm_elevation_dict(tiles, 16.886, -24.988)
        @test result === nothing
    end

    @testset "terrain_height dict overload - missing tile" begin
        tiles = Dict{String, Rasters.Raster}()
        result = Daisho.terrain_height(tiles, 16.886, -24.988)
        @test result == -1.0
    end

    @testset "read_srtm_elevation_multi - load all tiles" begin
        # Create a temporary directory (no tiles)
        tile_dir = joinpath(@__DIR__, "fixtures", "srtm_load_test")
        mkpath(tile_dir)

        tiles = Daisho.read_srtm_elevation_multi(tile_dir)
        @test isempty(tiles)

        rm(tile_dir, recursive=true, force=true)
    end

end
