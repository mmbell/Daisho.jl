# Quality Control functions

function fix_SEAPOL_RHOHV!(volume, moment_dict)

    # SEAPOL RHOHV has artificially high RHOHV due to noise subtraction
    # have to get original RHOHV from this variable
    # per Brenda J: LROSE adds bias of 128 to field, then shifts by -32768
    # need to add 32640 to field to return to original values (10.19.22)
    # diff from Brody's original SEA-POL code (LROSE conversion difference)
    volume.moments[:,moment_dict["RHOHV"]] .= volume.moments[:,moment_dict["UNKNOWN_ID_82"]]
    volume.moments[:,moment_dict["RHOHV"]] .+= 32640
    volume.moments[:,moment_dict["RHOHV"]] .= (volume.moments[:,moment_dict["RHOHV"]].-1) ./65533.0

end

function threshold_qc(raw_moments, raw_moment_dict, qc_moments, qc_moment_dict,
    sqi_threshold, snr_threshold, rhohv_threshold, spec_width_threshold)

    # Probably faster/better to do a logical indexing mask by gate 
    # rather than looping through all of these keys to do the same thing
    for key in keys(qc_moment_dict)
        
        if key == "SQI" || key == "PHIDP"
            # Don't QC this field
            continue
        end

        # Less than for SQI, SNR, RHOHV
        qc_moments[:,qc_moment_dict[key]] = ifelse.(isless.(raw_moments[:,raw_moment_dict["SQI"]],sqi_threshold),
            missing, qc_moments[:,qc_moment_dict[key]])
        qc_moments[:,qc_moment_dict[key]] = ifelse.(isless.(raw_moments[:,raw_moment_dict["SNR"]],snr_threshold),
            missing, qc_moments[:,qc_moment_dict[key]])
        qc_moments[:,qc_moment_dict[key]] = ifelse.(isless.(raw_moments[:,raw_moment_dict["RHOHV"]],rhohv_threshold),
            missing, qc_moments[:,qc_moment_dict[key]])

        # Greater than for spectrum width
        qc_moments[:,qc_moment_dict[key]] = ifelse.(isless.(raw_moments[:,raw_moment_dict["WIDTH"]],spec_width_threshold),
            qc_moments[:,qc_moment_dict[key]], missing)

    end

    return qc_moments

end

function smooth_sqi(sqi)

    # Smooth the SQI with a Gaussian kernel
    # This function helps to remove second trip echo
    
    
end

function despeckle(speckle, moments, moment_dict, n_gates, n_rays)

    # Despeckle a volume of data
    0 < speckle || error("Speckle definition must be greater than 0")
    n_moments = length(moment_dict)
    for m in 1:n_moments
        moment_data = reshape(moments[:,m],n_gates,n_rays)
        for i in 1:n_rays
            raydata = moment_data[:,i]
            n = 1
            while n < n_gates
                if !ismissing(raydata[n])
                    # Found a good gate, go forward to look for missing
                    feature_size = 0
                    s = n
                    while (n < n_gates) && !ismissing(raydata[n])
                        n += 1
                        feature_size += 1
                    end
                    
                    if feature_size > speckle
                        # Feature is larger than a speckle
                        continue
                    end
                    
                    while (s < n_gates) && (s <= n)
                        moment_data[s,i] = missing
                        s += 1
                    end
                else
                    n += 1
                end
            end
        end
        moments[:,m] .= moment_data[:]
    end

    return moments
    
end

function stddev_phidp_threshold(moments, moment_dict, n_gates, n_rays, window = 11, threshold = 12)

    # Use a standard deviation of phidp to threshold
    0 < window || error("Window definition must be greater than 0")
    window % 2 == 1 || error("Window must be an odd number")
    window = Int((window - 1)/2)
    phidp = ifelse.(isless.(moments[:,moment_dict["PHIDP"]],0.0), moments[:,moment_dict["PHIDP"]].+360.0, moments[:,moment_dict["PHIDP"]])
    phidp = reshape(phidp,n_gates,n_rays)
    phidp_mask = falses(size(phidp))
    for i in 1:n_rays
        raydata = phidp[:,i]
        for n = 1:n_gates
            stddev = std(raydata[max(n-window,1):min(n+window,n_gates)])
            if isless(threshold, stddev)
                phidp_mask[n,i] = true
            end
        end
    end

    # Apply the mask
    for m in 1:length(moment_dict)
        moments[phidp_mask[:],m] .= missing
    end

    return moments

end

function remove_platform_motion!(volume, moments, moment_dict)

    # Get rid of platform motion
    n_gates = length(volume.range)
    n_rays =  length(volume.azimuth)
    vel = reshape(moments[:,moment_dict["VEL"]],n_gates,n_rays)

    # Get the platform motion
    for i in 1:n_rays
        az = volume.azimuth[i] * pi/180.0
        el = volume.elevation[i] * pi/180.0
        ew_platform = volume.ew_platform[i]
        ns_platform = volume.ns_platform[i]
        w_platform = volume.w_platform[i]
        nyquist_velocity = volume.nyquist_velocity[i]
        platform_vr = ew_platform * cos(el) * sin(az) + ns_platform * cos(el) * cos(az) + w_platform * sin(el)
        vel[:,i] .= vel[:,i] .+ platform_vr
        println("Ray: ", i, " Platform velocity: ", platform_vr)
        # Put in the Nyquist interval
        vel[:,i] = ifelse.(isless.(vel[:,i], -nyquist_velocity), vel[:,i] .+ 2.0 .* nyquist_velocity, vel[:,i])
        vel[:,i] = ifelse.(isless.(nyquist_velocity, vel[:,i]), vel[:,i] .- 2.0 .* nyquist_velocity, vel[:,i])
    end

    # Put the velocity back in the volume
    moments[:,moment_dict["VEL"]] .= vel[:]
    return moments

end

function threshold_dbz(volume, raw_moment_dict, qc_moments, qc_moment_dict,
    dbz_threshold, vel_threshold, sw_threshold)

    count = 0
    bad_count = 0
    real_count = 0
    raw_moments = volume.moments
    beam_info = get_beam_info(volume)
    # Check for dBZ values greater than the threshold
    for i in 1:size(raw_moments,1)
        if ismissing(raw_moments[i,raw_moment_dict["DBZ_QC_FINAL"]])
            continue
        end
        height = beam_info[i,4]
        count += 1
        if (raw_moments[i,raw_moment_dict["DBZ_QC_FINAL"]] >= dbz_threshold)

            if (abs(raw_moments[i,raw_moment_dict["VEL_QC_FINAL"]]) <= vel_threshold
                && raw_moments[i,raw_moment_dict["WIDTH"]] <= sw_threshold
                && height < 500.0)
            #if height < 500.0
                #raw_moments[i,raw_moment_dict["KDP"]] > 0.0 || raw_moments[i,raw_moment_dict["HID_CSU"]] < 7
                bad_count += 1
                #println("Likely clutter at $i  exceeding $(dbz_threshold) at height $(height) m")
                #println(beam_info[i,:])
                for key in keys(qc_moment_dict)
                    if key == "SQI_FOR_MASK" || key == "PID_FOR_QC"
                        # Don't QC this field
                        continue
                    end
                    #println("$key: $(raw_moments[i,raw_moment_dict[key]])")
                    qc_moments[i,:] .= missing
                end
            else
                real_count += 1
                println("Maybe real at $i  exceeding $(dbz_threshold) at height $(height) m")
                println(beam_info[i,:])
                for key in keys(raw_moment_dict)
                    println("$key: $(raw_moments[i,raw_moment_dict[key]])")
                 end
            end
        end
    end

    println("Total gates: ", count, ", Bad gates: ", bad_count, ", Real gates: ", real_count)
    if count > 0
        println(" Bad Percentage: ", round(bad_count/count * 100, digits=2), "%", ", Real Percentage: ", round(real_count/count * 100, digits=2), "%")
    end

#    for key in keys(qc_moment_dict)
#        
#        if key == "SQI" || key == "PHIDP"
#            # Don't QC this field
#            continue
#        end
#
#        qc_moments[:,qc_moment_dict[key]] = ifelse.(isless.(raw_moments[:,raw_moment_dict["DBZ"]],dbz_threshold),
#            qc_moments[:,qc_moment_dict[key]], missing)
#        println("Thresholding ", key, " with dBZ threshold: ", dbz_threshold)
#
#    end

    return qc_moments

end

function threshold_height(volume, raw_moment_dict, qc_moments, qc_moment_dict,
    height_threshold)

    raw_moments = volume.moments
    beam_info = get_beam_info(volume)
    # Check for height values lower than the threshold
    for i in 1:size(raw_moments,1)
        if ismissing(raw_moments[i,raw_moment_dict["DBZ"]])
            continue
        end
        height = beam_info[i,4]
        if (height < height_threshold)
            for key in keys(qc_moment_dict)
                if key == "SQI_FOR_MASK" || key == "PID_FOR_QC"
                    # Don't QC this field
                    continue
                end
                qc_moments[i,qc_moment_dict[key]] = missing
            end
        end
    end
    return qc_moments

end

function threshold_terrain_height(volume, raw_moment_dict, qc_moments, qc_moment_dict, terrain_threshold, tile_dir)

    # Moment data
    raw_moments = volume.moments

    # Set the reference to the first location in the volume, but could be a parameter
    reference_latitude = volume.latitude[1]
    reference_longitude = volume.longitude[1]

    # Convert the relevant radar information to arrays
    TM = CoordRefSystems.shift(TransverseMercator{1.0,reference_latitude,WGS84Latest}, lonₒ= reference_longitude)
    grid_origin, radar_zyx, beams = radar_arrays(reference_latitude, reference_longitude, volume, TM)

    # Load the SRTM tiles
    tiles = read_srtm_elevation_multi(tile_dir)

    Threads.@threads for i in 1:size(beams,1)

        if ismissing(raw_moments[i,raw_moment_dict["SQI_FOR_MASK"]])
            # In the blanking sector or out of range
            continue
        end
        dbz = raw_moments[i,raw_moment_dict["DBZ"]]
        # Get the horizontal locations of every gate in Y, X dimension
        surface_range = Reff * asin(beams[i,3] * cos(beams[i,2]) / (Reff + beams[i,4]))
        y = grid_origin.y + (radar_zyx[i][2] + surface_range * cos(beams[i,1])) * 1u"m"
        x = grid_origin.x + (radar_zyx[i][3] + surface_range * sin(beams[i,1])) * 1u"m"
        #cartTM = convert(TM,Cartesian{WGS84Latest}(x, y))
        #latlon = convert(LatLon,cartTM)
        lat, lon = appx_inverse_projection(reference_latitude, reference_longitude, [ustrip.(y), ustrip.(x)])
        #dtm_height = terrain_height(tile_dir, ustrip.(latlon.lat), ustrip.(latlon.lon))
        dtm_height = terrain_height(tiles, lat, lon)
        az = 180.0 * beams[i,1] / pi
        el = 180.0 * beams[i,2] / pi
        range = beams[i,3]/1000.0
        #if abs(az - 289.0) < 1.0 && abs(range - 17.0) < 0.5
        #    println("Ray $i Az: $az El:$el, Range: $range over water! x: $(x), y: $(y), lat $(lat), lon $(lon), height: $(dtm_height), dBZ: $(dbz)")
        #end
        if (dtm_height > terrain_threshold)
            #println("Gate $i over land at lat $(latlon.lat), lon $(latlon.lon), height $dtm_height m")
            for key in keys(qc_moment_dict)
                # Remove everything including SQI_FOR_MASK and PID_FOR_QC
                qc_moments[i,qc_moment_dict[key]] = missing
            end
        elseif !ismissing(dbz) && (dbz >= 50.0)
            #println("High dBZ gate Az: $az El:$el over water! x: $(x), y: $(y), height: $(dtm_height), dBZ: $(dbz)")
        end
    end

    return qc_moments
end

function add_azimuthal_offset(volume, az_offset)

    # Add constant offset to azimuth
    volume.azimuth[:] .= volume.azimuth[:] .+ az_offset
    return volume.azimuth[:]

end

function mask_sector(volume, raw_moment_dict, qc_moments, qc_moment_dict, heading, az_min, az_max)

    # Mask data between two heading relative azimuths
    # Create an array with all the relevant beam info
    beams = [ (heading[j], volume.azimuth[j])
            for i in eachindex(volume.range), j in eachindex(volume.elevation) ]
    beams = [ beams[i][j] for i in eachindex(beams), j in 1:2]
    for i in 1:size(qc_moments,1)
        if ismissing(qc_moments[i,qc_moment_dict["SQI_FOR_MASK"]])
            # Already in the blanking sector or out of range
            continue
        end
        az = beams[i,2] - beams[i,1]
        if (az < 0.0)
            az += 360.0
        end
        if (az >= az_min) && (az <= az_max)
            #println("Gate $i over land at lat $(latlon.lat), lon $(latlon.lon), height $dtm_height m")
            for key in keys(qc_moment_dict)
                # Remove everything including SQI_FOR_MASK and PID_FOR_QC
                qc_moments[i,qc_moment_dict[key]] = missing
            end
        end
    end

    return qc_moments
end
