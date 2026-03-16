@testset "NetCDF Parameters" begin

    @testset "variable_attrib_dict" begin
        # Check that all expected moment names are present
        expected_keys = ["DBZ", "ZDR", "RHOHV", "KDP", "WIDTH", "VEL", "PHIDP", "SNR", "SQI"]
        for key in expected_keys
            @test haskey(Daisho.variable_attrib_dict, key)
        end

        # Level 1 variants
        l1_keys = ["DBZ_L1", "ZDR_L1", "RHOHV_L1", "KDP_L1", "WIDTH_L1", "VEL_L1"]
        for key in l1_keys
            @test haskey(Daisho.variable_attrib_dict, key)
        end

        # Special fields
        @test haskey(Daisho.variable_attrib_dict, "RATE_CSU_BLENDED")
        @test haskey(Daisho.variable_attrib_dict, "PID")
        @test haskey(Daisho.variable_attrib_dict, "PID_FOR_QC")
        @test haskey(Daisho.variable_attrib_dict, "HID_CSU")
        @test haskey(Daisho.variable_attrib_dict, "UNKNOWN")
    end

    @testset "common_attrib" begin
        @test haskey(Daisho.common_attrib, "_FillValue")
        @test haskey(Daisho.common_attrib, "missing_value")
        @test haskey(Daisho.common_attrib, "coordinates")
        @test haskey(Daisho.common_attrib, "grid_mapping")
        @test haskey(Daisho.common_attrib, "valid_min")
        @test haskey(Daisho.common_attrib, "valid_max")
        @test Daisho.common_attrib["_FillValue"] == Float32(-32768.0)
    end

    @testset "attribute structure" begin
        # Each variable attrib dict should have standard_name, long_name, units
        for (key, attrib) in Daisho.variable_attrib_dict
            @test haskey(attrib, "standard_name")
            @test haskey(attrib, "long_name")
            @test haskey(attrib, "units")
        end
    end

    @testset "DBZ attributes" begin
        @test Daisho.dbz_attrib["standard_name"] == "DBZ"
        @test Daisho.dbz_attrib["units"] == "dBZ"
    end

    @testset "VEL attributes" begin
        @test Daisho.vel_attrib["standard_name"] == "VEL"
        @test Daisho.vel_attrib["units"] == "m/s"
    end

end
