# Radar utilities

Base.@kwdef struct radar
    scan_name::String
    azimuth::Array{Union{Missing, Float32}}
    elevation::Array{Union{Missing, Float32}}
    range::Array{Union{Missing, Float32}}
    time::Array{DateTime}
    latitude::Array{Union{Missing, Float32}}
    longitude::Array{Union{Missing, Float32}}
    altitude::Array{Union{Missing, Float32}}
    fixed_angles::Array{Union{Missing, Float32}} 
    swpstart::Array{Union{Missing, Float32}}
    swpend::Array{Union{Missing, Float32}}
    moments::Array{Union{Missing, Float64}}
end

function initialize_moment_dictionaries(raw_moment_names, qc_moment_names)

    # This function maps names to numbers to facilitate array operations
    
    # Raw moments are recorded by the radar and have yet to be processed
    raw_moment_dict = Dict()
    for i in eachindex(raw_moment_names)
        raw_moment_dict[raw_moment_names[i]] = i
    end

    # QCed moments are post-processed names
    qc_moment_dict = Dict()
    for i in eachindex(qc_moment_names)
        qc_moment_dict[qc_moment_names[i]] = i
    end

    return raw_moment_dict, qc_moment_dict
end

function initialize_qc_fields(volume, raw_moment_dict, qc_moment_dict)

    # Initialize QC fields by copying over the originals
    n_moments = length(qc_moment_dict)
    qc_moments = Array{Union{Missing, Float32}}(undef,length(volume.range)*length(volume.azimuth),n_moments)

    # This loop assumes that all the QC moments originate as raw ones
    # This may or may not be a good assumption in the future, but need to create these fields somehow
    for key in keys(qc_moment_dict)
        qc_moments[:, qc_moment_dict[key]] .= volume.moments[:,raw_moment_dict[key]]
    end
    
    return qc_moments
    
end

function split_sweeps(radar_volume)

    # Extract sweeps from the volume and put it in an array of radar structures
    swp_indices = [ range(Int(radar_volume.swpstart[i] + 1),Int(radar_volume.swpend[i] + 1)) for i in eachindex(radar_volume.swpstart) ]

    sweeps = Array{radar}(undef,length(radar_volume.swpstart))
    for i in eachindex(swp_indices)

        swp = swp_indices[i]
        scan_name = radar_volume.scan_name
        azimuth = view(radar_volume.azimuth, swp)
        elevation = view(radar_volume.elevation, swp)
        range = radar_volume.range
        time = view(radar_volume.time, swp)
        latitude = view(radar_volume.latitude, swp)
        longitude = view(radar_volume.longitude, swp)
        altitude = view(radar_volume.altitude, swp)
        fixed_angles = view(radar_volume.fixed_angles, i)
        swpstart = view(radar_volume.swpstart, i)
        swpend = view(radar_volume.swpend, i)
        
        m_start = Int(length(radar_volume.range)*(radar_volume.swpstart[i]) + 1)
        m_end = Int(length(radar_volume.range)*(radar_volume.swpend[i] + 1))
        moments = view(radar_volume.moments, m_start:m_end, :)
        sweeps[i] = radar(scan_name, azimuth, elevation, range, time,
            latitude, longitude, altitude, fixed_angles, swpstart, swpend, moments)
    end

    return sweeps
end

function beam_height(slant_range, elevation, radar_height)

    # This function calculates the beam height using 4/3 Earth standard refraction model
    # Slant range is in meters
    # Elevation is in degrees
    Reff = 4.0 * 6371000.0 / 3.0
    elrad = deg2rad(elevation)
    h = sqrt(slant_range^2 + Reff^2 + 2*slant_range*Reff*sin(elrad)) - Reff + radar_height
    return h

end

function dB_to_linear(moment)

    # Convert Z or ZDR or other moment in dB to linear units
    return 10.0 .^ (moment[:] ./ 10.0)

end

function linear_to_dB(moment)

    # Convert Z or ZDR or other moment in linear units to dB
    return 10.0 .* log10.(moment[:])

end

function dB_to_linear!(moment)

    # Convert Z or ZDR or other moment in dB to linear units
    moment[:] .= 10.0 .^ (moment[:] ./ 10.0)

end

function linear_to_dB!(moment)

    # Convert Z or ZDR or other moment in linear units to dB
    moment[:] .= 10.0 .* log10.(moment[:])

end

function read_cfradial(file, moment_dict)

    inputds = NCDataset(file);

    scan_name = inputds.attrib["scan_name"]
    
    az = inputds["azimuth"]
    azdata = az[:]
    
    el = inputds["elevation"]
    eldata = el[:]
    
    range = inputds["range"]
    rangedata = range[:]
    
    time = inputds["time"]
    timedata = time[:]

    latitude = inputds["latitude"]
    latdata = latitude[:]
    
    longitude = inputds["longitude"]
    londata = longitude[:]

    altitude = inputds["altitude"]
    altdata = altitude[:]
    
    fixed_angle = inputds["fixed_angle"]
    angledata = fixed_angle[:]

    sweep_start = inputds["sweep_start_ray_index"]
    swpstart = sweep_start[:]
    
    sweep_end = inputds["sweep_end_ray_index"]
    swpend = sweep_end[:]

    # Store radar data
    n_moments = length(moment_dict)

    # This could mess things up in the future since some CfRadial files have time, range
    # But variable gated volumes have n_points dimension
    # Need some checks on this
    # n_points = inputds.dim["n_points"]
    # If number of gates vary then dimension is smaller
    # if n_points
    #    moment_dim = n_points
    # end
    
    radardata = Array{Union{Missing, Float32}}(undef,length(range)*length(time),n_moments)
    for key in keys(moment_dict)
        radardata[:, moment_dict[key]] = inputds[key][:]
    end

    close(inputds)
    return radar(scan_name, azdata, eldata, rangedata, timedata, latdata, londata, altdata, angledata, swpstart, swpend, radardata)

end

function write_qced_cfradial(file, qc_file, qc_moments, qc_moment_dict)

    inputds = NCDataset(file);
    
    # Write the metadata to a new file
    println("Writing $qc_file...")
    
    ds = NCDataset(qc_file,"c", attrib = OrderedDict(
        "Conventions"               => inputds.attrib["Conventions"],
        "Sub_conventions"           => inputds.attrib["Sub_conventions"],
        "version"                   => inputds.attrib["version"],
        "title"                     => inputds.attrib["title"],
        "institution"               => inputds.attrib["institution"],
        "references"                => inputds.attrib["references"],
        "source"                    => inputds.attrib["source"],
        "history"                   => inputds.attrib["history"],
        "comment"                   => inputds.attrib["comment"],
        "original_format"           => inputds.attrib["original_format"],
        "driver"                    => inputds.attrib["driver"],
        "created"                   => inputds.attrib["created"],
        "start_datetime"            => inputds.attrib["start_datetime"],
        "time_coverage_start"       => inputds.attrib["time_coverage_start"],
        "start_time"                => inputds.attrib["start_time"],
        "end_datetime"              => inputds.attrib["end_datetime"],
        "time_coverage_end"         => inputds.attrib["time_coverage_end"],
        "end_time"                  => inputds.attrib["end_time"],
        "instrument_name"           => inputds.attrib["instrument_name"],
        "site_name"                 => inputds.attrib["site_name"],
        "scan_name"                 => inputds.attrib["scan_name"],
        "scan_id"                   => Int32(inputds.attrib["scan_id"]),
        "platform_is_mobile"        => inputds.attrib["platform_is_mobile"],
        "n_gates_vary"              => inputds.attrib["n_gates_vary"],
        "ray_times_increase"        => inputds.attrib["ray_times_increase"],
    ))
    
    # Dimensions
    
    ds.dim["time"] = inputds.dim["time"]
    ds.dim["range"] = inputds.dim["range"]
    ds.dim["sweep"] = inputds.dim["sweep"]
    ds.dim["string_length_8"] = inputds.dim["string_length_8"]
    ds.dim["string_length_32"] = inputds.dim["string_length_32"]
    ds.dim["status_xml_length"] = inputds.dim["status_xml_length"]
    ds.dim["r_calib"] = inputds.dim["r_calib"]
    ds.dim["frequency"] = inputds.dim["frequency"]
    
    # Declare variables
    
    ncvolume_number = defVar(ds,"volume_number", Int32, (), attrib = OrderedDict(
        "long_name"                 => "data_volume_index_number",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncplatform_type = defVar(ds,"platform_type", Char, ("string_length_32",), attrib = OrderedDict(
        "long_name"                 => "platform_type",
        "options"                   => "fixed, vehicle, ship, aircraft_fore, aircraft_aft, aircraft_tail, aircraft_belly, aircraft_roof, aircraft_nose, satellite_orbit, satellite_geostat",
    ))
    
    ncprimary_axis = defVar(ds,"primary_axis", Char, ("string_length_32",), attrib = OrderedDict(
        "long_name"                 => "primary_axis_of_rotation",
        "options"                   => "axis_z, axis_y, axis_x, axis_z_prime, axis_y_prime, axis_x_prime",
    ))
    
    ncstatus_xml = defVar(ds,"status_xml", Char, ("status_xml_length",), attrib = OrderedDict(
        "long_name"                 => "status_of_instrument",
    ))
    
    ncinstrument_type = defVar(ds,"instrument_type", Char, ("string_length_32",), attrib = OrderedDict(
        "long_name"                 => "type_of_instrument",
        "options"                   => "radar, lidar",
        "meta_group"                => "instrument_parameters",
    ))
    
    ncradar_antenna_gain_h = defVar(ds,"radar_antenna_gain_h", Float32, (), attrib = OrderedDict(
        "long_name"                 => "nominal_radar_antenna_gain_h_channel",
        "units"                     => "db",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncradar_antenna_gain_v = defVar(ds,"radar_antenna_gain_v", Float32, (), attrib = OrderedDict(
        "long_name"                 => "nominal_radar_antenna_gain_v_channel",
        "units"                     => "db",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncradar_beam_width_h = defVar(ds,"radar_beam_width_h", Float32, (), attrib = OrderedDict(
        "long_name"                 => "half_power_radar_beam_width_h_channel",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncradar_beam_width_v = defVar(ds,"radar_beam_width_v", Float32, (), attrib = OrderedDict(
        "long_name"                 => "half_power_radar_beam_width_v_channel",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncradar_rx_bandwidth = defVar(ds,"radar_rx_bandwidth", Float32, (), attrib = OrderedDict(
        "long_name"                 => "radar_receiver_bandwidth",
        "units"                     => "s-1",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    nctime_coverage_start = defVar(ds,"time_coverage_start", Char, ("string_length_32",), attrib = OrderedDict(
        "long_name"                 => "data_volume_start_time_utc",
        "comment"                   => "ray times are relative to start time in secs",
    ))
    
    nctime_coverage_end = defVar(ds,"time_coverage_end", Char, ("string_length_32",), attrib = OrderedDict(
        "long_name"                 => "data_volume_end_time_utc",
    ))
    
    ncfrequency = defVar(ds,"frequency", Float32, ("frequency",), attrib = OrderedDict(
        "long_name"                 => "transmission_frequency",
        "units"                     => "s-1",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncgrid_mapping = defVar(ds,"grid_mapping", Int32, (), attrib = OrderedDict(
        "grid_mapping_name"         => "radar_lidar_radial_scan",
    ))
    
    ncsweep_number = defVar(ds,"sweep_number", Int32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "sweep_index_number_0_based",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncsweep_mode = defVar(ds,"sweep_mode", Char, ("string_length_32", "sweep"), attrib = OrderedDict(
        "long_name"                 => "scan_mode_for_sweep",
        "options"                   => "sector, coplane, rhi, vertical_pointing, idle, azimuth_surveillance, elevation_surveillance, sunscan, pointing, calibration, manual_ppi, manual_rhi, sunscan_rhi, doppler_beam_swinging, complex_trajectory, electronic_steering",
    ))
    
    ncpolarization_mode = defVar(ds,"polarization_mode", Char, ("string_length_32", "sweep"), attrib = OrderedDict(
        "long_name"                 => "polarization_mode_for_sweep",
        "options"                   => "horizontal, vertical, hv_alt, hv_sim, circular",
        "meta_group"                => "radar_parameters",
    ))
    
    ncprt_mode = defVar(ds,"prt_mode", Char, ("string_length_32", "sweep"), attrib = OrderedDict(
        "long_name"                 => "transmit_pulse_mode",
        "options"                   => "fixed, staggered, dual",
        "meta_group"                => "radar_parameters",
    ))
    
    ncfollow_mode = defVar(ds,"follow_mode", Char, ("string_length_32", "sweep"), attrib = OrderedDict(
        "long_name"                 => "follow_mode_for_scan_strategy",
        "options"                   => "none, sun, vehicle, aircraft, target, manual",
        "meta_group"                => "instrument_parameters",
    ))
    
    ncfixed_angle = defVar(ds,"fixed_angle", Float32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "ray_target_fixed_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    nctarget_scan_rate = defVar(ds,"target_scan_rate", Float32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "target_scan_rate_for_sweep",
        "units"                     => "degrees per second",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncsweep_start_ray_index = defVar(ds,"sweep_start_ray_index", Int32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "index_of_first_ray_in_sweep",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncsweep_end_ray_index = defVar(ds,"sweep_end_ray_index", Int32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "index_of_last_ray_in_sweep",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncrays_are_indexed = defVar(ds,"rays_are_indexed", Char, ("string_length_8", "sweep"), attrib = OrderedDict(
        "long_name"                 => "flag_for_indexed_rays",
    ))
    
    ncray_angle_res = defVar(ds,"ray_angle_res", Float32, ("sweep",), attrib = OrderedDict(
        "long_name"                 => "angular_resolution_between_rays",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_time = defVar(ds,"r_calib_time", Char, ("string_length_32", "r_calib"), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_time_utc",
        "meta_group"                => "radar_calibration",
    ))
    
    ncr_calib_pulse_width = defVar(ds,"r_calib_pulse_width", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_pulse_width",
        "units"                     => "seconds",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_xmit_power_h = defVar(ds,"r_calib_xmit_power_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_xmit_power_h_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_xmit_power_v = defVar(ds,"r_calib_xmit_power_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_xmit_power_v_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_two_way_waveguide_loss_h = defVar(ds,"r_calib_two_way_waveguide_loss_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_two_way_waveguide_loss_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_two_way_waveguide_loss_v = defVar(ds,"r_calib_two_way_waveguide_loss_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_two_way_waveguide_loss_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_two_way_radome_loss_h = defVar(ds,"r_calib_two_way_radome_loss_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_two_way_radome_loss_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_two_way_radome_loss_v = defVar(ds,"r_calib_two_way_radome_loss_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_two_way_radome_loss_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_mismatch_loss = defVar(ds,"r_calib_receiver_mismatch_loss", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_receiver_mismatch_loss",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_k_squared_water = defVar(ds,"r_calib_k_squared_water", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_k_squared_water",
        "units"                     => "",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_radar_constant_h = defVar(ds,"r_calib_radar_constant_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_constant_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_radar_constant_v = defVar(ds,"r_calib_radar_constant_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_constant_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_antenna_gain_h = defVar(ds,"r_calib_antenna_gain_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_antenna_gain_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_antenna_gain_v = defVar(ds,"r_calib_antenna_gain_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_antenna_gain_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_hc = defVar(ds,"r_calib_noise_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_noise_h_co_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_vc = defVar(ds,"r_calib_noise_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_noise_v_co_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_hx = defVar(ds,"r_calib_noise_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_noise_h_cross_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_vx = defVar(ds,"r_calib_noise_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_noise_v_cross_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_i0_dbm_hc = defVar(ds,"r_calib_i0_dbm_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_i0_dbm_hc",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_i0_dbm_vc = defVar(ds,"r_calib_i0_dbm_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_i0_dbm_vc",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_i0_dbm_hx = defVar(ds,"r_calib_i0_dbm_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_i0_dbm_hx",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_i0_dbm_vx = defVar(ds,"r_calib_i0_dbm_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_i0_dbm_vx",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_gain_hc = defVar(ds,"r_calib_receiver_gain_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_gain_h_co_polar_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_gain_vc = defVar(ds,"r_calib_receiver_gain_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_gain_v_co_polar_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_gain_hx = defVar(ds,"r_calib_receiver_gain_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_gain_h_cross_polar_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_gain_vx = defVar(ds,"r_calib_receiver_gain_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_gain_v_cross_polar_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_slope_hc = defVar(ds,"r_calib_receiver_slope_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_slope_h_co_polar_channel",
        "units"                     => "",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_slope_vc = defVar(ds,"r_calib_receiver_slope_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_slope_v_co_polar_channel",
        "units"                     => "",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_slope_hx = defVar(ds,"r_calib_receiver_slope_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_slope_h_cross_polar_channel",
        "units"                     => "",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_receiver_slope_vx = defVar(ds,"r_calib_receiver_slope_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_receiver_slope_v_cross_polar_channel",
        "units"                     => "",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_dynamic_range_db_hc = defVar(ds,"r_calib_dynamic_range_db_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_dynamic_range_db_hc",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_dynamic_range_db_vc = defVar(ds,"r_calib_dynamic_range_db_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_dynamic_range_db_vc",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_dynamic_range_db_hx = defVar(ds,"r_calib_dynamic_range_db_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_dynamic_range_db_hx",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_dynamic_range_db_vx = defVar(ds,"r_calib_dynamic_range_db_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "r_calib_dynamic_range_db_vx",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_base_dbz_1km_hc = defVar(ds,"r_calib_base_dbz_1km_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_reflectivity_at_1km_at_zero_snr_h_co_polar_channel",
        "units"                     => "dBZ",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_base_dbz_1km_vc = defVar(ds,"r_calib_base_dbz_1km_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_reflectivity_at_1km_at_zero_snr_v_co_polar_channel",
        "units"                     => "dBZ",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_base_dbz_1km_hx = defVar(ds,"r_calib_base_dbz_1km_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_reflectivity_at_1km_at_zero_snr_h_cross_polar_channel",
        "units"                     => "dBZ",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_base_dbz_1km_vx = defVar(ds,"r_calib_base_dbz_1km_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_reflectivity_at_1km_at_zero_snr_v_cross_polar_channel",
        "units"                     => "dBZ",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_sun_power_hc = defVar(ds,"r_calib_sun_power_hc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_sun_power_h_co_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_sun_power_vc = defVar(ds,"r_calib_sun_power_vc", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_sun_power_v_co_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_sun_power_hx = defVar(ds,"r_calib_sun_power_hx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_sun_power_h_cross_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_sun_power_vx = defVar(ds,"r_calib_sun_power_vx", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_sun_power_v_cross_polar_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_source_power_h = defVar(ds,"r_calib_noise_source_power_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_noise_source_power_h_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_noise_source_power_v = defVar(ds,"r_calib_noise_source_power_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_noise_source_power_v_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_power_measure_loss_h = defVar(ds,"r_calib_power_measure_loss_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_power_measurement_loss_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_power_measure_loss_v = defVar(ds,"r_calib_power_measure_loss_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_power_measurement_loss_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_coupler_forward_loss_h = defVar(ds,"r_calib_coupler_forward_loss_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_coupler_forward_loss_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_coupler_forward_loss_v = defVar(ds,"r_calib_coupler_forward_loss_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_coupler_forward_loss_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_dbz_correction = defVar(ds,"r_calib_dbz_correction", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_dbz_correction",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_zdr_correction = defVar(ds,"r_calib_zdr_correction", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_zdr_correction",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_ldr_correction_h = defVar(ds,"r_calib_ldr_correction_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_ldr_correction_h_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_ldr_correction_v = defVar(ds,"r_calib_ldr_correction_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_ldr_correction_v_channel",
        "units"                     => "db",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_system_phidp = defVar(ds,"r_calib_system_phidp", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "calibrated_radar_system_phidp",
        "units"                     => "degrees",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_test_power_h = defVar(ds,"r_calib_test_power_h", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_test_power_h_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncr_calib_test_power_v = defVar(ds,"r_calib_test_power_v", Float32, ("r_calib",), attrib = OrderedDict(
        "long_name"                 => "radar_calibration_test_power_v_channel",
        "units"                     => "dBm",
        "meta_group"                => "radar_calibration",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "time",
        "long_name"                 => "time in seconds since volume start",
        "calendar"                  => "gregorian",
        "units"                     => "seconds since 2024-08-27T21:40:08Z",
        "comment"                   => "times are relative to the volume start_time",
    ))
    
    ncrange = defVar(ds,"range", Float32, ("range",), attrib = OrderedDict(
        "long_name"                 => "range_to_center_of_measurement_volume",
        "standard_name"             => "projection_range_coordinate",
        "units"                     => "meters",
        "axis"                      => "radial_range_coordinate",
        "spacing_is_constant"       => "true",
        "meters_to_center_of_first_gate" => 400.0000059604645,
        "meters_between_gates"      => 100.00000149011612,
    ))
    
    ncray_n_gates = defVar(ds,"ray_n_gates", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "number_of_gates",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncray_start_index = defVar(ds,"ray_start_index", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "array_index_to_start_of_ray",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncray_start_range = defVar(ds,"ray_start_range", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "start_range_for_ray",
        "units"                     => "meters",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncray_gate_spacing = defVar(ds,"ray_gate_spacing", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "gate_spacing_for_ray",
        "units"                     => "meters",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncazimuth = defVar(ds,"azimuth", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "ray_azimuth_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncelevation = defVar(ds,"elevation", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "ray_elevation_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
        "positive"                  => "up",
    ))
    
    ncpulse_width = defVar(ds,"pulse_width", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "transmitter_pulse_width",
        "units"                     => "seconds",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncprt = defVar(ds,"prt", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "pulse_repetition_time",
        "units"                     => "seconds",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncprt_ratio = defVar(ds,"prt_ratio", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "pulse_repetition_frequency_ratio",
        "units"                     => "seconds",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncnyquist_velocity = defVar(ds,"nyquist_velocity", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "unambiguous_doppler_velocity",
        "units"                     => "meters per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncunambiguous_range = defVar(ds,"unambiguous_range", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "unambiguous_range",
        "units"                     => "meters",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncantenna_transition = defVar(ds,"antenna_transition", Int8, ("time",), attrib = OrderedDict(
        "long_name"                 => "antenna_is_in_transition_between_sweeps",
        "units"                     => "",
        "_FillValue"                => Int8(-128),
        "comment"                   => "1 if antenna is in transition, 0 otherwise",
    ))
    
    ncgeorefs_applied = defVar(ds,"georefs_applied", Int8, ("time",), attrib = OrderedDict(
        "long_name"                 => "georefs_have_been_applied_to_ray",
        "units"                     => "",
        "_FillValue"                => Int8(-128),
        "comment"                   => "1 if georefs have been applied, 0 otherwise",
    ))
    
    ncn_samples = defVar(ds,"n_samples", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "number_of_samples_used_to_compute_moments",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncr_calib_index = defVar(ds,"r_calib_index", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "calibration_data_array_index_per_ray",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
        "meta_group"                => "radar_calibration",
        "comment"                   => "This is the index for the calibration which applies to this ray",
    ))
    
    ncmeasured_transmit_power_h = defVar(ds,"measured_transmit_power_h", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "measured_radar_transmit_power_h_channel",
        "units"                     => "dBm",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncmeasured_transmit_power_v = defVar(ds,"measured_transmit_power_v", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "measured_radar_transmit_power_v_channel",
        "units"                     => "dBm",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "radar_parameters",
    ))
    
    ncscan_rate = defVar(ds,"scan_rate", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "antenna_angle_scan_rate",
        "units"                     => "degrees per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "instrument_parameters",
    ))
    
    ncgeoref_time = defVar(ds,"georef_time", Float64, ("time",), attrib = OrderedDict(
        "long_name"                 => "georef time in seconds since volume start",
        "units"                     => "seconds",
        "_FillValue"                => -9999.0,
    ))
    
    ncgeoref_unit_num = defVar(ds,"georef_unit_num", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "georef hardware unit number",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    ncgeoref_unit_id = defVar(ds,"georef_unit_id", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "georef hardware id or serial number",
        "units"                     => "",
        "_FillValue"                => Int32(-9999),
    ))
    
    nclatitude = defVar(ds,"latitude", Float64, ("time",), attrib = OrderedDict(
        "long_name"                 => "latitude",
        "units"                     => "degrees_north",
        "_FillValue"                => -9999.0,
    ))
    
    nclongitude = defVar(ds,"longitude", Float64, ("time",), attrib = OrderedDict(
        "long_name"                 => "longitude",
        "units"                     => "degrees_east",
        "_FillValue"                => -9999.0,
    ))
    
    ncaltitude = defVar(ds,"altitude", Float64, ("time",), attrib = OrderedDict(
        "long_name"                 => "altitude",
        "units"                     => "meters",
        "_FillValue"                => -9999.0,
        "positive"                  => "up",
    ))
    
    ncaltitude_agl = defVar(ds,"altitude_agl", Float64, ("time",), attrib = OrderedDict(
        "long_name"                 => "altitude_above_ground_level",
        "units"                     => "meters",
        "_FillValue"                => -9999.0,
        "positive"                  => "up",
    ))
    
    ncheading = defVar(ds,"heading", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_heading_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncroll = defVar(ds,"roll", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_roll_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    ncpitch = defVar(ds,"pitch", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_pitch_angle",
        "units"                     => "degrees",
        "_FillValue"                => Float32(-9999.0),
    ))
    
    nceastward_velocity = defVar(ds,"eastward_velocity", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_eastward_velocity",
        "units"                     => "meters per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "platform_velocity",
    ))
    
    ncnorthward_velocity = defVar(ds,"northward_velocity", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_northward_velocity",
        "units"                     => "meters per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "platform_velocity",
    ))
    
    ncvertical_velocity = defVar(ds,"vertical_velocity", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_vertical_velocity",
        "units"                     => "meters per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "platform_velocity",
    ))
    
    ncheading_change_rate = defVar(ds,"heading_change_rate", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_heading_angle_rate_of_change",
        "units"                     => "degrees per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "platform_velocity",
    ))
    
    ncpitch_change_rate = defVar(ds,"pitch_change_rate", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "platform_pitch_angle_rate_of_change",
        "units"                     => "degrees per second",
        "_FillValue"                => Float32(-9999.0),
        "meta_group"                => "platform_velocity",
    ))
    
    ncDBZ_TOT = defVar(ds,"DBZ_TOT_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "DBZ_TOT",
        "standard_name"             => "DBZ_TOT",
        "units"                     => "dBZ",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncDBZ = defVar(ds,"DBZ_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "DBZ",
        "standard_name"             => "DBZ",
        "units"                     => "dBZ",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncVEL = defVar(ds,"VEL_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "VEL",
        "standard_name"             => "VEL",
        "units"                     => "m/s",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncWIDTH = defVar(ds,"WIDTH_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "WIDTH",
        "standard_name"             => "WIDTH",
        "units"                     => "m/s",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(327.68),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncZDR = defVar(ds,"ZDR_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "ZDR",
        "standard_name"             => "ZDR",
        "units"                     => "dB",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncKDP = defVar(ds,"KDP_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "KDP",
        "standard_name"             => "KDP",
        "units"                     => "deg/km",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncRHOHV = defVar(ds,"RHOHV_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "RHOHV",
        "standard_name"             => "RHOHV",
        "units"                     => "",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(1.5259488e-5),
        "add_offset"                => Float32(0.5000076),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncSQI = defVar(ds,"SQI", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "SQI",
        "standard_name"             => "SQI",
        "units"                     => "",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(1.5259488e-5),
        "add_offset"                => Float32(0.5000076),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncPHIDP = defVar(ds,"PHIDP_L1", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "PHIDP",
        "standard_name"             => "PHIDP",
        "units"                     => "deg",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.0054933317),
        "add_offset"                => Float32(180.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncSNR = defVar(ds,"SNR", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "SNR",
        "standard_name"             => "SNR",
        "units"                     => "dB",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncUNKNOWN_ID_82 = defVar(ds,"UNKNOWN_ID_82", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "UNKNOWN_ID_82",
        "standard_name"             => "UNKNOWN_ID_82",
        "units"                     => "",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(1.0),
        "add_offset"                => Float32(128.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncDBZ_QC = defVar(ds,"DBZ", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "DBZ_QC",
        "standard_name"             => "DBZ_QC",
        "units"                     => "dBZ",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncVEL_QC = defVar(ds,"VEL", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "VEL_QC",
        "standard_name"             => "VEL_QC",
        "units"                     => "m/s",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))

    ncWIDTH_QC = defVar(ds,"WIDTH", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "WIDTH_QC",
        "standard_name"             => "WIDTH_QC",
        "units"                     => "m/s",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(327.68),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncZDR_QC = defVar(ds,"ZDR", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "ZDR_QC",
        "standard_name"             => "ZDR_QC",
        "units"                     => "dB",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncKDP_QC = defVar(ds,"KDP", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "KDP_QC",
        "standard_name"             => "KDP_QC",
        "units"                     => "deg/km",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.01),
        "add_offset"                => Float32(0.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncPHIDP_QC = defVar(ds,"PHIDP", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "PHIDP_QC",
        "standard_name"             => "PHIDP_QC",
        "units"                     => "deg",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(0.0054933317),
        "add_offset"                => Float32(180.0),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    ncRHOHV_QC = defVar(ds,"RHOHV", Int16, ("range", "time"), attrib = OrderedDict(
        "long_name"                 => "RHOHV_QC",
        "standard_name"             => "RHOHV_QC",
        "units"                     => "",
        "sampling_ratio"            => Float32(1.0),
        "_FillValue"                => Int16(-32768),
        "scale_factor"              => Float32(1.5259488e-5),
        "add_offset"                => Float32(0.5000076),
        "grid_mapping"              => "grid_mapping",
        "coordinates"               => "elevation azimuth range heading roll pitch rotation tilt",
    ))
    
    # Define variables
    ncvolume_number[:] = inputds["volume_number"][:]
    ncplatform_type[:] = inputds["platform_type"][:]
    ncprimary_axis[:] = inputds["primary_axis"][:]
    ncstatus_xml[:] = inputds["status_xml"][:]
    ncinstrument_type[:] = inputds["instrument_type"][:]
    ncradar_antenna_gain_h[:] = inputds["radar_antenna_gain_h"][:]
    ncradar_antenna_gain_v[:] = inputds["radar_antenna_gain_v"][:]
    ncradar_beam_width_h[:] = inputds["radar_beam_width_h"][:]
    ncradar_beam_width_v[:] = inputds["radar_beam_width_v"][:]
    ncradar_rx_bandwidth[:] = inputds["radar_rx_bandwidth"][:]
    nctime_coverage_start[:] = inputds["time_coverage_start"][:]
    nctime_coverage_end[:] = inputds["time_coverage_end"][:]
    ncfrequency[:] = inputds["frequency"][:]
    ncgrid_mapping[:] = inputds["grid_mapping"][:]
    ncsweep_number[:] = inputds["sweep_number"][:]
    ncsweep_mode[:] = inputds["sweep_mode"][:]
    ncpolarization_mode[:] = inputds["polarization_mode"][:]
    ncprt_mode[:] = inputds["prt_mode"][:]
    ncfollow_mode[:] = inputds["follow_mode"][:]
    ncfixed_angle[:] = inputds["fixed_angle"][:]
    nctarget_scan_rate[:] = inputds["target_scan_rate"][:]
    ncsweep_start_ray_index[:] = inputds["sweep_start_ray_index"][:]
    ncsweep_end_ray_index[:] = inputds["sweep_end_ray_index"][:]
    ncrays_are_indexed[:] = inputds["rays_are_indexed"][:]
    ncray_angle_res[:] = inputds["ray_angle_res"][:]
    ncr_calib_time[:] = inputds["r_calib_time"][:]
    ncr_calib_pulse_width[:] = inputds["r_calib_pulse_width"][:]
    ncr_calib_xmit_power_h[:] = inputds["r_calib_xmit_power_h"][:]
    ncr_calib_xmit_power_v[:] = inputds["r_calib_xmit_power_v"][:]
    ncr_calib_two_way_waveguide_loss_h[:] = inputds["r_calib_two_way_waveguide_loss_h"][:]
    ncr_calib_two_way_waveguide_loss_v[:] = inputds["r_calib_two_way_waveguide_loss_v"][:]
    ncr_calib_two_way_radome_loss_h[:] = inputds["r_calib_two_way_radome_loss_h"][:]
    ncr_calib_two_way_radome_loss_v[:] = inputds["r_calib_two_way_radome_loss_v"][:]
    ncr_calib_receiver_mismatch_loss[:] = inputds["r_calib_receiver_mismatch_loss"][:]
    ncr_calib_k_squared_water[:] = inputds["r_calib_k_squared_water"][:]
    ncr_calib_radar_constant_h[:] = inputds["r_calib_radar_constant_h"][:]
    ncr_calib_radar_constant_v[:] = inputds["r_calib_radar_constant_v"][:]
    ncr_calib_antenna_gain_h[:] = inputds["r_calib_antenna_gain_h"][:]
    ncr_calib_antenna_gain_v[:] = inputds["r_calib_antenna_gain_v"][:]
    ncr_calib_noise_hc[:] = inputds["r_calib_noise_hc"][:]
    ncr_calib_noise_vc[:] = inputds["r_calib_noise_vc"][:]
    ncr_calib_noise_hx[:] = inputds["r_calib_noise_hx"][:]
    ncr_calib_noise_vx[:] = inputds["r_calib_noise_vx"][:]
    ncr_calib_i0_dbm_hc[:] = inputds["r_calib_i0_dbm_hc"][:]
    ncr_calib_i0_dbm_vc[:] = inputds["r_calib_i0_dbm_vc"][:]
    ncr_calib_i0_dbm_hx[:] = inputds["r_calib_i0_dbm_hx"][:]
    ncr_calib_i0_dbm_vx[:] = inputds["r_calib_i0_dbm_vx"][:]
    ncr_calib_receiver_gain_hc[:] = inputds["r_calib_receiver_gain_hc"][:]
    ncr_calib_receiver_gain_vc[:] = inputds["r_calib_receiver_gain_vc"][:]
    ncr_calib_receiver_gain_hx[:] = inputds["r_calib_receiver_gain_hx"][:]
    ncr_calib_receiver_gain_vx[:] = inputds["r_calib_receiver_gain_vx"][:]
    ncr_calib_receiver_slope_hc[:] = inputds["r_calib_receiver_slope_hc"][:]
    ncr_calib_receiver_slope_vc[:] = inputds["r_calib_receiver_slope_vc"][:]
    ncr_calib_receiver_slope_hx[:] = inputds["r_calib_receiver_slope_hx"][:]
    ncr_calib_receiver_slope_vx[:] = inputds["r_calib_receiver_slope_vx"][:]
    ncr_calib_dynamic_range_db_hc[:] = inputds["r_calib_dynamic_range_db_hc"][:]
    ncr_calib_dynamic_range_db_vc[:] = inputds["r_calib_dynamic_range_db_vc"][:]
    ncr_calib_dynamic_range_db_hx[:] = inputds["r_calib_dynamic_range_db_hx"][:]
    ncr_calib_dynamic_range_db_vx[:] = inputds["r_calib_dynamic_range_db_vx"][:]
    ncr_calib_base_dbz_1km_hc[:] = inputds["r_calib_base_dbz_1km_hc"][:]
    ncr_calib_base_dbz_1km_vc[:] = inputds["r_calib_base_dbz_1km_vc"][:]
    ncr_calib_base_dbz_1km_hx[:] = inputds["r_calib_base_dbz_1km_hx"][:]
    ncr_calib_base_dbz_1km_vx[:] = inputds["r_calib_base_dbz_1km_vx"][:]
    ncr_calib_sun_power_hc[:] = inputds["r_calib_sun_power_hc"][:]
    ncr_calib_sun_power_vc[:] = inputds["r_calib_sun_power_vc"][:]
    ncr_calib_sun_power_hx[:] = inputds["r_calib_sun_power_hx"][:]
    ncr_calib_sun_power_vx[:] = inputds["r_calib_sun_power_vx"][:]
    ncr_calib_noise_source_power_h[:] = inputds["r_calib_noise_source_power_h"][:]
    ncr_calib_noise_source_power_v[:] = inputds["r_calib_noise_source_power_v"][:]
    ncr_calib_power_measure_loss_h[:] = inputds["r_calib_power_measure_loss_h"][:]
    ncr_calib_power_measure_loss_v[:] = inputds["r_calib_power_measure_loss_v"][:]
    ncr_calib_coupler_forward_loss_h[:] = inputds["r_calib_coupler_forward_loss_h"][:]
    ncr_calib_coupler_forward_loss_v[:] = inputds["r_calib_coupler_forward_loss_v"][:]
    ncr_calib_dbz_correction[:] = inputds["r_calib_dbz_correction"][:]
    ncr_calib_zdr_correction[:] = inputds["r_calib_zdr_correction"][:]
    ncr_calib_ldr_correction_h[:] = inputds["r_calib_ldr_correction_h"][:]
    ncr_calib_ldr_correction_v[:] = inputds["r_calib_ldr_correction_v"][:]
    ncr_calib_system_phidp[:] = inputds["r_calib_system_phidp"][:]
    ncr_calib_test_power_h[:] = inputds["r_calib_test_power_h"][:]
    ncr_calib_test_power_v[:] = inputds["r_calib_test_power_v"][:]
    nctime[:] = inputds["time"][:]
    ncrange[:] = inputds["range"][:]
    ncray_start_range[:] = inputds["ray_start_range"][:]
    ncray_gate_spacing[:] = inputds["ray_gate_spacing"][:]
    ncazimuth[:] = inputds["azimuth"][:]
    ncelevation[:] = inputds["elevation"][:]
    ncpulse_width[:] = inputds["pulse_width"][:]
    ncprt[:] = inputds["prt"][:]
    ncprt_ratio[:] = inputds["prt_ratio"][:]
    ncnyquist_velocity[:] = inputds["nyquist_velocity"][:]
    ncunambiguous_range[:] = inputds["unambiguous_range"][:]
    ncantenna_transition[:] = inputds["antenna_transition"][:]
    ncgeorefs_applied[:] = inputds["georefs_applied"][:]
    ncn_samples[:] = inputds["n_samples"][:]
    ncr_calib_index[:] = inputds["r_calib_index"][:]
    ncmeasured_transmit_power_h[:] = inputds["measured_transmit_power_h"][:]
    ncmeasured_transmit_power_v[:] = inputds["measured_transmit_power_v"][:]
    ncscan_rate[:] = inputds["scan_rate"][:]
    nclatitude[:] = inputds["latitude"][:]
    nclongitude[:] = inputds["longitude"][:]
    ncaltitude[:] = inputds["altitude"][:]
    ncaltitude_agl[:] = inputds["altitude_agl"][:]
    ncheading[:] = inputds["heading"][:]
    ncroll[:] = inputds["roll"][:]
    ncpitch[:] = inputds["pitch"][:]
    #nceastward_velocity[:] = inputds["eastward_velocity"][:]
    #ncnorthward_velocity[:] = inputds["northward_velocity"][:]
    #ncvertical_velocity[:] = inputds["vertical_velocity"][:]
    #ncheading_change_rate[:] = inputds["heading_change_rate"][:]
    #ncpitch_change_rate[:] = inputds["pitch_change_rate"][:]
    ncDBZ_TOT[:] = inputds["DBZ_TOT"][:]
    ncDBZ[:] = inputds["DBZ"][:]
    ncVEL[:] = inputds["VEL"][:]
    ncWIDTH[:] = inputds["WIDTH"][:]
    ncZDR[:] = inputds["ZDR"][:]
    ncKDP[:] = inputds["KDP"][:]
    ncRHOHV[:] = inputds["RHOHV"][:]
    ncSQI[:] = inputds["SQI"][:]
    ncPHIDP[:] = inputds["PHIDP"][:]
    ncSNR[:] = inputds["SNR"][:]
    #ncUNKNOWN_ID_82[:] = inputds["UNKNOWN_ID_82"][:]
    
    # QCed variables
    ncDBZ_QC[:] = qc_moments[:, qc_moment_dict["DBZ"]]
    ncZDR_QC[:] = qc_moments[:, qc_moment_dict["ZDR"]]
    ncVEL_QC[:] = qc_moments[:, qc_moment_dict["VEL"]]
    ncKDP_QC[:] = qc_moments[:, qc_moment_dict["KDP"]]
    ncRHOHV_QC[:] = qc_moments[:, qc_moment_dict["RHOHV"]]
    ncPHIDP_QC[:] = qc_moments[:, qc_moment_dict["PHIDP"]]
    ncWIDTH_QC[:] = qc_moments[:, qc_moment_dict["WIDTH"]]
    #ncSMOOTH_SQI[:] = qc_moments[qc_moment_dict["SQI"],:]
    
    close(ds)

end

