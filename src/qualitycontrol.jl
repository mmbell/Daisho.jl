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
