@testset "Radar" begin

    @testset "radar struct construction" begin
        vol = make_synthetic_radar()
        @test vol.scan_name == "TEST_VOL"
        @test length(vol.azimuth) == 10
        @test length(vol.range) == 20
        @test size(vol.moments) == (200, 3)
        @test length(vol.time) == 10
        @test vol.latitude[1] ≈ 16.886f0
        @test vol.longitude[1] ≈ -24.988f0
    end

    @testset "radar struct with custom parameters" begin
        vol = make_synthetic_radar(n_rays=5, n_gates=10, n_moments=4, lat=40.0, lon=-105.0)
        @test length(vol.azimuth) == 5
        @test length(vol.range) == 10
        @test size(vol.moments) == (50, 4)
        @test vol.latitude[1] ≈ 40.0f0
    end

    @testset "initialize_moment_dictionaries" begin
        raw_names = ["DBZ", "VEL", "WIDTH", "SQI"]
        qc_names = ["DBZ", "VEL", "WIDTH"]
        grid_types = [:linear, :weighted, :weighted]

        raw_dict, qc_dict, grid_dict = Daisho.initialize_moment_dictionaries(raw_names, qc_names, grid_types)

        @test raw_dict["DBZ"] == 1
        @test raw_dict["VEL"] == 2
        @test raw_dict["WIDTH"] == 3
        @test raw_dict["SQI"] == 4
        @test qc_dict["DBZ"] == 1
        @test qc_dict["VEL"] == 2
        @test qc_dict["WIDTH"] == 3
        @test grid_dict[1] == :linear
        @test grid_dict[2] == :weighted
        @test grid_dict[3] == :weighted
    end

    @testset "initialize_moment_dictionaries - empty input" begin
        raw_dict, qc_dict, grid_dict = Daisho.initialize_moment_dictionaries([], [], [])
        @test isempty(raw_dict)
        @test isempty(qc_dict)
        @test isempty(grid_dict)
    end

    @testset "initialize_qc_fields" begin
        vol = make_synthetic_radar(n_rays=5, n_gates=10, n_moments=3)
        # Set some known values
        vol.moments[:, 1] .= 20.0  # DBZ
        vol.moments[:, 2] .= 5.0   # VEL
        vol.moments[:, 3] .= 2.0   # WIDTH

        raw_dict = make_moment_dict(["DBZ", "VEL", "WIDTH"])
        qc_dict = make_moment_dict(["DBZ", "VEL"])

        qc_moments = Daisho.initialize_qc_fields(vol, raw_dict, qc_dict)

        @test size(qc_moments) == (50, 2)
        @test all(qc_moments[:, qc_dict["DBZ"]] .== 20.0)
        @test all(qc_moments[:, qc_dict["VEL"]] .== 5.0)
    end

    @testset "split_sweeps" begin
        vol = make_synthetic_radar(n_rays=10, n_gates=5, n_moments=2, n_sweeps=2)
        sweeps = Daisho.split_sweeps(vol)

        @test length(sweeps) == 2
        @test length(sweeps[1].azimuth) == 5
        @test length(sweeps[2].azimuth) == 5
        @test sweeps[1].scan_name == vol.scan_name
    end

    @testset "split_sweeps - single sweep" begin
        vol = make_synthetic_radar(n_rays=10, n_gates=5, n_moments=2, n_sweeps=1)
        sweeps = Daisho.split_sweeps(vol)

        @test length(sweeps) == 1
        @test length(sweeps[1].azimuth) == 10
    end

    @testset "beam_height" begin
        # At 0 elevation, close range, height should be approximately radar_height
        h = Daisho.beam_height(0.0, 0.0, 50.0)
        @test h ≈ 50.0

        # At 0 elevation, 100 km range, beam height should be above ground due to Earth curvature
        h = Daisho.beam_height(100000.0, 0.0, 0.0)
        @test h > 0.0

        # At 90 degrees elevation, beam height ≈ slant_range + radar_height
        h = Daisho.beam_height(10000.0, 90.0, 50.0)
        @test h ≈ 10050.0 atol=10.0

        # Typical case: 1 degree elevation, 50 km range
        h = Daisho.beam_height(50000.0, 1.0, 0.0)
        @test h > 800.0  # Should be roughly 870m + curvature
        @test h < 1200.0

        # Negative elevation (possible for airborne radar)
        h = Daisho.beam_height(10000.0, -5.0, 5000.0)
        @test h < 5000.0  # Height should be less than radar height
    end

    @testset "dB_to_linear" begin
        # 0 dB = 1.0
        result = Daisho.dB_to_linear([0.0])
        @test result[1] ≈ 1.0

        # 10 dB = 10.0
        result = Daisho.dB_to_linear([10.0])
        @test result[1] ≈ 10.0

        # -10 dB = 0.1
        result = Daisho.dB_to_linear([-10.0])
        @test result[1] ≈ 0.1

        # 20 dB = 100.0
        result = Daisho.dB_to_linear([20.0])
        @test result[1] ≈ 100.0

        # Multiple values
        result = Daisho.dB_to_linear([0.0, 10.0, 20.0])
        @test result ≈ [1.0, 10.0, 100.0]
    end

    @testset "linear_to_dB" begin
        result = Daisho.linear_to_dB([1.0])
        @test result[1] ≈ 0.0

        result = Daisho.linear_to_dB([10.0])
        @test result[1] ≈ 10.0

        result = Daisho.linear_to_dB([100.0])
        @test result[1] ≈ 20.0
    end

    @testset "dB_to_linear / linear_to_dB round-trip" begin
        original = [0.0, 5.0, 10.0, 20.0, 30.0, -5.0]
        result = Daisho.linear_to_dB(Daisho.dB_to_linear(original))
        @test result ≈ original atol=1e-10
    end

    @testset "dB_to_linear!" begin
        moment = [0.0, 10.0, 20.0]
        Daisho.dB_to_linear!(moment)
        @test moment ≈ [1.0, 10.0, 100.0]
    end

    @testset "linear_to_dB!" begin
        moment = [1.0, 10.0, 100.0]
        Daisho.linear_to_dB!(moment)
        @test moment ≈ [0.0, 10.0, 20.0]
    end

    @testset "CfRadial I/O round-trip" begin
        filepath = joinpath(@__DIR__, "fixtures", "test_cfradial.nc")
        moment_names = ["DBZ", "VEL", "WIDTH"]
        moment_dict = make_moment_dict(moment_names)

        # Create a synthetic CfRadial file
        create_synthetic_cfradial(filepath, n_sweeps=2, n_rays=10, n_gates=20,
            moment_names=moment_names)

        # Read it back
        vol = Daisho.read_cfradial(filepath, moment_dict)

        @test vol.scan_name == "TEST_VOL"
        @test length(vol.azimuth) == 10
        @test length(vol.range) == 20
        @test size(vol.moments) == (200, 3)
        @test vol.latitude[1] ≈ 16.886 atol=0.01
        @test vol.longitude[1] ≈ -24.988 atol=0.01

        # Clean up
        rm(filepath, force=true)
    end

    @testset "get_radar_orientation - no heading/pitch/roll" begin
        # Create a CfRadial file without heading/pitch/roll
        filepath = joinpath(@__DIR__, "fixtures", "test_orientation.nc")
        create_synthetic_cfradial(filepath, n_sweeps=1, n_rays=5, n_gates=10)

        result = Daisho.get_radar_orientation(filepath)
        @test size(result, 1) == 5
        @test size(result, 2) == 3
        # All should be NaN since no heading/pitch/roll in synthetic file
        @test all(isnan.(result[:, 1]))
        @test all(isnan.(result[:, 2]))
        @test all(isnan.(result[:, 3]))

        rm(filepath, force=true)
    end

end
