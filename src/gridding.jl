# Radar gridding functions

function initialize_regular_grid(xmin, xincr, xdim, ymin, yincr, ydim, zmin, zincr, zdim)

    # Define and allocate a 3d regular grid
    regular_3d_grid = Array{Float64}(undef,zdim,ydim,xdim,3)
    for i in CartesianIndices(size(regular_3d_grid)[1:3])
        regular_3d_grid[i,1] = zincr * (i[1]-1) + zmin
        regular_3d_grid[i,2] = yincr * (i[2]-1) + ymin
        regular_3d_grid[i,3] = xincr * (i[3]-1) + xmin
    end
    return regular_3d_grid
end

function initialize_regular_grid(xmin, xincr, xdim, ymin, yincr, ydim)

    # Define and allocate a 2d regular grid
    regular_2d_grid = Array{Float64}(undef,ydim,xdim,2)
    for i in CartesianIndices(size(regular_2d_grid)[1:2])
        regular_2d_grid[i,1] = yincr * (i[1]-1) + ymin
        regular_2d_grid[i,2] = xincr * (i[2]-1) + xmin
    end
    return regular_2d_grid
end

function get_radar_zyx(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, radar_volume::radar, projection)

    # Radar locations mapped to transverse mercator with Z,Y,X dimensions
    radar_loc = convert.(projection,LatLon.(radar_volume.latitude, radar_volume.longitude))
    beam_origin = [ [radar_volume.altitude[i], Float64(ustrip(radar_loc[i].y)), Float64(ustrip(radar_loc[i].x))]  for i in eachindex(radar_loc)]
    radar_zyx = [ zyx for j in radar_volume.range, zyx in beam_origin ]
    return radar_zyx

end

function get_beam_info(radar_volume::radar)

    # Create an array with all the relevant beam info (azimuth, elevation, range, height)
    beams = [ (deg2rad(radar_volume.azimuth[j]), deg2rad(radar_volume.elevation[j]), radar_volume.range[i],
            beam_height(radar_volume.range[i], radar_volume.elevation[j], radar_volume.altitude[j])) 
            for i in eachindex(radar_volume.range), j in eachindex(radar_volume.elevation) ]
    beams = [ beams[i][j] for i in eachindex(beams), j in 1:4]
    return beams

end

function radar_arrays(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, radar_volume::radar, projection)

    # Grid origin
    grid_origin = convert(projection,LatLon(reference_latitude, reference_longitude))
    
    # Radar locations mapped to transverse mercator with Z,Y,X dimensions
    radar_zyx = get_radar_zyx(reference_latitude, reference_longitude, radar_volume, projection)

    # Create an array with all the relevant beam info (azimuth, elevation, range, height)
    beams = get_beam_info(radar_volume)

    return grid_origin, radar_zyx, beams
    
end

function radar_balltree_yx(radar_volume::radar, radar_zyx::AbstractArray, beams::AbstractArray)

    # Create a balltree that has horizontal locations of every gate in Y, X dimension
    gate_yx = zeros(Float64, 2, length(radar_volume.azimuth)*length(radar_volume.range))
    for i in 1:size(beams,1)
        surface_range = Reff * asin(beams[i,3] * cos(beams[i,2]) / (Reff + beams[i,4]))
        gate_yx[1,i] = radar_zyx[i][2] + surface_range * cos(beams[i,1])
        gate_yx[2,i] = radar_zyx[i][3] + surface_range * sin(beams[i,1])
    end
    balltree = BallTree(gate_yx)
    return balltree

end

function radar_balltree_r(radar_volume::radar, radar_zyx::AbstractArray, beams::AbstractArray)

    # Create a balltree that has horizontal locations of every gate in R dimension
    gate_r = zeros(Float64, 1, length(radar_volume.azimuth)*length(radar_volume.range))
    for i in 1:size(beams,1)
        surface_range = Reff * asin(beams[i,3] * cos(beams[i,2]) / (Reff + beams[i,4]))
        y = radar_zyx[i][2] + surface_range * cos(beams[i,1])
        x = radar_zyx[i][3] + surface_range * sin(beams[i,1])
        gate_r[i] = sqrt(x^2 + y^2)
    end
    balltree = BallTree(gate_r)
    return balltree

end

function appx_inverse_projection(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, yx_point::AbstractArray)

    # Approximate lat/lon from SAMURAI formula
    # This formula originated from HRD code, but is originally from FCC wireless communication specs evidently 
    latrad = reference_latitude * pi/180.0
    fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
        + 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad)
    fac_lon = 111.41513 * cos(latrad)
        - 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad)
    lon = reference_longitude + yx_point[2]/(fac_lon*1000.0)
    lat = reference_latitude + yx_point[1]/(fac_lat*1000.0)
    return lat, lon

end

function grid_radar_volume(radar_volume, moment_dict, grid_type_dict, output_file,
        xmin, xincr, xdim, ymin, yincr, ydim, zmin, zincr, zdim, beam_inflation, power_threshold,
        missing_key="SQI", valid_key="DBZ")

    # Set the reference to the first location in the volume, but could be a parameter
    reference_latitude = radar_volume.latitude[1]
    reference_longitude = radar_volume.longitude[1]
    
    # Initialize the gridpoints
    # This array is slightly different than Springsteel spectral arrays, need to reconcile later
    gridpoints = initialize_regular_grid(xmin, xincr, xdim, ymin, yincr, ydim, zmin, zincr, zdim)

    h_roi = xincr * 0.75
    v_roi = zincr * 0.75
    
    radar_grid, latlon_grid = grid_volume(reference_latitude, reference_longitude, gridpoints, 
        radar_volume, moment_dict, grid_type_dict, h_roi, v_roi, beam_inflation, power_threshold,
        missing_key="SQI", valid_key="DBZ")

    write_gridded_radar_volume(output_file, radar_volume.time[1],
        radar_volume.time[end], gridpoints, radar_grid, latlon_grid, moment_dict,
        reference_latitude, reference_longitude)
    
end

function grid_radar_rhi(radar_volume, moment_dict, grid_type_dict, output_file,
        rmin, rincr, rdim, rhi_zmin, rhi_zincr, rhi_zdim, beam_inflation, power_threshold,
        missing_key::String="SQI", valid_key::String="DBZ")

    # Set the reference to the first location in the volume, but could be a parameter
    reference_latitude = radar_volume.latitude[1]
    reference_longitude = radar_volume.longitude[1]
    
    # Initialize the gridpoints
    # This array is slightly different than Springsteel spectral arrays, need to reconcile later
    gridpoints = initialize_regular_grid(rmin, rincr, rdim, rhi_zmin, rhi_zincr, rhi_zdim)

    h_roi = rincr * 0.75
    v_roi = rhi_zincr * 0.75
    
    radar_grid, latlon_grid = grid_rhi(reference_latitude, reference_longitude, gridpoints, 
        radar_volume, moment_dict, grid_type_dict, h_roi, v_roi, beam_inflation, power_threshold, missing_key, valid_key)

    write_gridded_radar_rhi(output_file, radar_volume.time[1],
        radar_volume.time[end], gridpoints, radar_grid, latlon_grid, moment_dict,
        reference_latitude, reference_longitude)
    
end

function grid_radar_ppi(radar_volume, moment_dict, grid_type_dict, output_file,
        xmin, xincr, xdim, ymin, yincr, ydim, beam_inflation, power_threshold,
        missing_key="SQI", valid_key="DBZ", heading=-9999.0)

    # Set the reference to the first location in the volume, but could be a parameter
    reference_latitude = radar_volume.latitude[1]
    reference_longitude = radar_volume.longitude[1]
    
    # Initialize the gridpoints
    # This array is slightly different than Springsteel spectral arrays, need to reconcile later
    gridpoints = initialize_regular_grid(xmin, xincr, xdim, ymin, yincr, ydim)

    h_roi = xincr * 0.75
    
    radar_grid, latlon_grid = grid_ppi(reference_latitude, reference_longitude, gridpoints, 
        radar_volume, moment_dict, grid_type_dict, h_roi, beam_inflation, power_threshold, missing_key, valid_key)

    write_gridded_radar_ppi(output_file, radar_volume.time[1],
        radar_volume.time[end], gridpoints, radar_grid, latlon_grid, moment_dict,
        reference_latitude, reference_longitude, heading)
    
end

function grid_radar_composite(radar_volume, moment_dict, grid_type_dict, output_file,
        xmin, xincr, xdim, ymin, yincr, ydim, beam_inflation,
        missing_key="SQI", valid_key="DBZ", mean_heading=-9999.0)

    # Set the reference to the first location in the volume, but could be a parameter
    reference_latitude = radar_volume.latitude[1]
    reference_longitude = radar_volume.longitude[1]
    
    # Initialize the gridpoints
    # This array is slightly different than Springsteel spectral arrays, need to reconcile later
    gridpoints = initialize_regular_grid(xmin, xincr, xdim, ymin, yincr, ydim)

    h_roi = xincr * 0.75
    
    radar_grid, latlon_grid = grid_composite(reference_latitude, reference_longitude, gridpoints, 
        radar_volume, moment_dict, grid_type_dict, h_roi, beam_inflation, missing_key, valid_key)

    write_gridded_radar_ppi(output_file, radar_volume.time[1],
        radar_volume.time[end], gridpoints, radar_grid, latlon_grid, moment_dict,
        reference_latitude, reference_longitude, mean_heading)
    
end

function grid_volume(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, gridpoints::AbstractArray, 
        radar_volume::radar, moment_dict::Dict, grid_type_dict::Dict, horizontal_roi::Float64, vertical_roi::Float64, beam_inflation::Float64,
        power_threshold::Float64, missing_key::String="SQI", valid_key::String="DBZ")

    # Convert the relevant radar information to arrays
    TM = CoordRefSystems.shift(TransverseMercator{1.0,reference_latitude,WGS84Latest}, lonₒ= reference_longitude)
    grid_origin, radar_zyx, beams = radar_arrays(reference_latitude, reference_longitude, radar_volume, TM)
    
    # Create a balltree that has horizontal locations of every gate in Y, X dimension
    balltree = radar_balltree_yx(radar_volume, radar_zyx, beams)

    # Allocate the grid for the radar moments and the weights for each radar gate
    n_moments = length(moment_dict)
    radar_grid = fill(-32768.0,n_moments,size(gridpoints,1),size(gridpoints,2),size(gridpoints,3))
    weights = zeros(Float64,n_moments,size(gridpoints,1),size(gridpoints,2),size(gridpoints,3))

    # Allocate a regular latlon grid for the map
    latlon_grid = Array{Float64}(undef,size(gridpoints,2),size(gridpoints,3),2)

    # Loop through the horizontal indices then do each column
    Threads.@threads for i in CartesianIndices(size(gridpoints)[2:3])

        # Calculate the lat, lon of the gridpoint
        yx_point = gridpoints[1,i,2:3]

        # Using CoordRefSystems pure Julia transform, the Proj.jl wrapper was significantly slower for unknown reasons
        cartTM = convert(TM,Cartesian{WGS84Latest}(grid_origin.x + yx_point[2]u"m", grid_origin.y + yx_point[1]u"m"))
        latlon = convert(LatLon,cartTM)
        latlon_grid[i,1] = ustrip(latlon.lat)
        latlon_grid[i,2] = ustrip(latlon.lon)
        
        # Find the points within the radius of influence of the horizontal gridpoint
        eff_h_radius_influence = horizontal_roi
        eff_v_radius_influence = vertical_roi
        if beam_inflation > 0.0
            origin_dist = euclidean(yx_point, [0.0, 0.0])
            eff_h_radius_influence = max(beam_inflation * origin_dist, horizontal_roi)
            eff_v_radius_influence = max(beam_inflation * origin_dist, vertical_roi)
        end
        gates = inrange(balltree, yx_point, eff_h_radius_influence)

        if !isempty(gates)

            # Found some gates that are within range horizontally
            for j in 1:size(gridpoints,1)

                # First check to see whether at least one gate is in range vertically
                grid_z = gridpoints[j,i,1]
                min_dist, min_idx = findmin(z -> abs(z - grid_z), beams[gates,4])
                if min_dist > eff_v_radius_influence
                    # There is no gate in range of this grid box
                    continue
                else
                    if !ismissing(radar_volume.moments[gates[min_idx], moment_dict[missing_key]])
                        # There is at least one gate in range so set the flags to -9999
                        for m in 1:n_moments
                            radar_grid[m,j,i] = -9999.0
                        end
                    end
                end

                # Loops through the nearby gates with valid data
                valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[valid_key]])))
                for gate in gates[valid_gates]
                    
                    # Calculate the beam intercept to the gridpoint
                    dz = gridpoints[j,i,1] - radar_zyx[gate][1]
                    r = beams[gate,3]

                    # Calculate the effective elevation of the gridpoint taking into account
                    # earth curvature and standard refraction
                    # Sine of the height angle using equation from Doviak and Zrnic (1993)
                    sine_h = ((dz + Reff)^2 - r^2 - Reff^2) / (2*r*Reff)
                    
                    gridpt_el = missing
                    if abs(sine_h) < 1.0
                        gridpt_el = asin(sine_h)
                    else
                        # Point is not accessible by the radar?
                        # I think this is failing because range exceeds height at origin
                        # A better check is needed, but skip for now
                        continue
                    end

                    # Calculate the effective azimuth
                    dx = yx_point[2] - radar_zyx[gate][3]
                    dy = yx_point[1] - radar_zyx[gate][2]
                    gridpt_az = (pi/2.0) - atan(dy, dx)
                    if gridpt_az < 0
                        gridpt_az += 2*pi
                    end

                    # Calculate the spherical angle difference using Haversine formula
                    angle_diff = spherical_angle([beams[gate,1], beams[gate,2]], 
                        [gridpt_az, gridpt_el])

                    # Use the half-power beamwidth to define the beam
                    # 79.43 = ln(0.5) / (0.5 deg * pi / 180.0)
                    # Center of beam = 1.0 angle_weight
                    angle_weight = exp(-angle_diff * 79.43)
                    if angle_weight < power_threshold
                        angle_weight = 0
                    end
    
                    # Range weighting based on center of gridbox = 1.0 range_weight
                    gridpt_r = sin(sqrt(dx^2 + dy^2)/Reff) * (Reff + dz) / cos(gridpt_el)
                    range_weight = gridpt_r / r

                    # Multiply weights so that center of beam is 1.0
                    total_weight = range_weight * angle_weight

                    # If there is non-zero weight, add it to the grid box
                    if total_weight > 0.0
                        for m in 1:n_moments
                            if weights[m,j,i] == 0.0
                                # Initialize the radar grid box with 0 since there is a possibility that the beam hit it
                                radar_grid[m,j,i] = 0.0
                            end

                            if !ismissing(radar_volume.moments[gate,m])
                                if grid_type_dict[m] == :linear
                                    linear_z = 10.0 ^ (radar_volume.moments[gate,m] / 10.0)
                                    radar_grid[m,j,i] += total_weight * linear_z
                                    weights[m,j,i] += total_weight
                                elseif grid_type_dict[m] == :nearest
                                    if total_weight > weights[m,j,i]
                                        radar_grid[m,j,i] = radar_volume.moments[gate,m]
                                        weights[m,j,i] = total_weight
                                    end
                                else
                                    radar_grid[m,j,i] += total_weight * radar_volume.moments[gate,m]
                                    weights[m,j,i] += total_weight
                                end
                            end
                        end
                    end

                end # End of gate loop

                # Divide by the total weight for that gridbox
                for m in 1:n_moments
                    if weights[m,j,i] > 0.0 && grid_type_dict[m] != :nearest
                        radar_grid[m,j,i] /= weights[m,j,i]
                        if grid_type_dict[m] == :linear
                            if radar_grid[m,j,i] > 0.0
                                radar_grid[m,j,i] = 10.0 * log10(radar_grid[m,j,i])
                            else
                                radar_grid[m,j,i] = -9999.0
                            end
                        end
                    elseif weights[m,j,i] == 0.0
                        if radar_grid[m,j,i] == 0.0
                            # Use a different flag to indicate data has been QCed out
                            radar_grid[m,j,i] = -9999.0
                        end
                    end
                end # End of moment loop

            end # End of height loop
        end # End of gate empty test
    end # End of horizontal loop

    # Return the gridded radar array
    return radar_grid, latlon_grid
end

function grid_rhi(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, gridpoints::AbstractArray, 
        radar_volume::radar, moment_dict::Dict, grid_type_dict::Dict, horizontal_roi::Float64, vertical_roi::Float64, beam_inflation::Float64,
        power_threshold::Float64, missing_key::String="SQI", valid_key::String="DBZ")

    # Convert the relevant radar information to arrays
    TM = CoordRefSystems.shift(TransverseMercator{1.0,reference_latitude,WGS84Latest}, lonₒ= reference_longitude)
    grid_origin, radar_zyx, beams = radar_arrays(reference_latitude, reference_longitude, radar_volume, TM)

    # Create a balltree that has horizontal locations of every gate in Y, X dimension
    balltree = radar_balltree_r(radar_volume, radar_zyx, beams)

    # Allocate the grid for the radar moments and the weights for each radar gate
    n_moments = length(moment_dict)
    radar_grid = fill(-32768.0,n_moments,size(gridpoints,1),size(gridpoints,2))
    weights = zeros(Float64,n_moments,size(gridpoints,1),size(gridpoints,2))

    # Allocate a regular latlon grid for the map
    latlon_grid = Array{Float64}(undef,size(gridpoints,2),2)

    # Loop through the horizontal indices then do each column
    Threads.@threads for i in 1:size(gridpoints)[2]

        # Calculate the lat, lon of the gridpoint
        yx_point = gridpoints[1,i,1:2]

        # Using CoordRefSystems pure Julia transform, the Proj.jl wrapper was significantly slower for unknown reasons
        cartTM = convert(TM,Cartesian{WGS84Latest}(grid_origin.x + yx_point[2]u"m", grid_origin.y + yx_point[1]u"m"))
        latlon = convert(LatLon,cartTM)
        latlon_grid[i,1] = ustrip(latlon.lat)
        latlon_grid[i,2] = ustrip(latlon.lon)

        # Find the points within the radius of influence of the horizontal gridpoint
        eff_h_radius_influence = horizontal_roi
        eff_v_radius_influence = vertical_roi
        origin_dist = euclidean(yx_point, [0.0, 0.0])
        if beam_inflation > 0.0
            # No beam inflation in horizontal in this case since we are gridding in range
            eff_v_radius_influence = max(beam_inflation * origin_dist, vertical_roi)
        end
        
        gates = inrange(balltree, [ origin_dist ], eff_h_radius_influence)
        if !isempty(gates)

            # Found some gates that are within range horizontally
            for j in 1:size(gridpoints,1)

                # First check to see whether at least one gate is in range vertically
                grid_z = gridpoints[j,i,1]
                min_dist, min_idx = findmin(z -> abs(z - grid_z), beams[gates,4])
                if min_dist > eff_v_radius_influence
                    # There is no gate in range of this grid box
                    continue
                else
                    if !ismissing(radar_volume.moments[gates[min_idx], moment_dict[missing_key]])
                        # There is at least one gate in range so set the flags to -9999
                        for m in 1:n_moments
                            radar_grid[m,j,i] = -9999.0
                        end
                    end
                end

                # Loops through the nearby gates with valid data
                valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[valid_key]])))
                for gate in gates[valid_gates]
                    
                    # Calculate the beam intercept to the gridpoint
                    dz = gridpoints[j,i,1] - radar_zyx[gate][1]
                    r = beams[gate,3]

                    # Calculate the effective elevation of the gridpoint taking into account
                    # earth curvature and standard refraction
                    # Sine of the height angle using full equation
                    sine_h = ((dz + Reff)^2 - r^2 - Reff^2) / (2*r*Reff)
        
                    # Use approximation from Gao et al. 2006
                    #sine_h = (z/r) - (r/(2*Reff))
                    
                    gridpt_el = missing
                    if abs(sine_h) < 1.0
                        gridpt_el = asin(sine_h)
                    else
                        # Point is not accessible by the radar?
                        # I think this is failing because range exceeds height at origin
                        # A better check is needed, but skip for now
                        continue
                    end

                    # Don't need azimuth, just dx and dy to get the surface range
                    dx = yx_point[2] - radar_zyx[gate][3]
                    dy = yx_point[1] - radar_zyx[gate][2]

                    # Calculate the spherical angle difference using Haversine formula
                    angle_diff = spherical_angle([beams[gate,1], beams[gate,2]], 
                        [beams[gate,1], gridpt_el])

                    # Use the half-power beamwidth to define the beam
                    # 79.43 = ln(0.5) / (0.5 deg * pi / 180.0)
                    # Center of beam = 1.0 angle_weight
                    angle_weight = exp(-angle_diff * 79.43)
                    if angle_weight < power_threshold
                        angle_weight = 0
                    end
    
                    # Range weighting based on center of gridbox = 1.0 range_weight
                    gridpt_r = sin(sqrt(dx^2 + dy^2)/Reff) * (Reff + dz) / cos(gridpt_el)
                    range_weight = gridpt_r / r

                    # Multiply weights so that center of beam is 1.0
                    total_weight = range_weight * angle_weight

                    # If there is non-zero weight, add it to the grid box
                    if total_weight > 0.0
                        for m in 1:n_moments
                            if weights[m,j,i] == 0.0
                                # Initialize the radar grid box with 0 since there is a possibility that the beam hit it
                                radar_grid[m,j,i] = 0.0
                            end

                            if !ismissing(radar_volume.moments[gate,m])
                                if grid_type_dict[m] == :linear
                                    linear_z = 10.0 ^ (radar_volume.moments[gate,m] / 10.0)
                                    radar_grid[m,j,i] += total_weight * linear_z
                                    weights[m,j,i] += total_weight
                                elseif grid_type_dict[m] == :nearest
                                    if total_weight > weights[m,j,i]
                                        radar_grid[m,j,i] = radar_volume.moments[gate,m]
                                        weights[m,j,i] = total_weight
                                    end
                                else
                                    radar_grid[m,j,i] += total_weight * radar_volume.moments[gate,m]
                                    weights[m,j,i] += total_weight
                                end
                            end
                        end
                    end

                end # End of gate loop

                # Divide by the total weight for that gridbox
                for m in 1:n_moments
                    if weights[m,j,i] > 0.0 && grid_type_dict[m] != :nearest
                        radar_grid[m,j,i] /= weights[m,j,i]
                        if grid_type_dict[m] == :linear
                            if radar_grid[m,j,i] > 0.0
                                radar_grid[m,j,i] = 10.0 * log10(radar_grid[m,j,i])
                            else
                                radar_grid[m,j,i] = -9999.0
                            end
                        end
                    elseif weights[m,j,i] == 0.0
                        if radar_grid[m,j,i] == 0.0
                            # Use a different flag to indicate data has been QCed out
                            radar_grid[m,j,i] = -9999.0
                        end
                    end
                end # End of moment loop

            end # End of height loop
        end # End of gate empty test
    end # End of horizontal loop

    # Return the gridded radar array
    return radar_grid, latlon_grid
end

function grid_ppi(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, gridpoints::AbstractArray, 
        radar_volume::radar, moment_dict::Dict, grid_type_dict::Dict, horizontal_roi::Float64, beam_inflation::Float64, 
        power_threshold::Float64, missing_key::String="SQI", valid_key::String="DBZ")

    # Convert the relevant radar information to arrays
    TM = CoordRefSystems.shift(TransverseMercator{1.0,reference_latitude,WGS84Latest}, lonₒ= reference_longitude)
    grid_origin, radar_zyx, beams = radar_arrays(reference_latitude, reference_longitude, radar_volume, TM)
    
    # Create a balltree that has horizontal locations of every gate in Y, X dimension
    balltree = radar_balltree_yx(radar_volume, radar_zyx, beams)

    # Allocate the grid for the radar moments and the weights for each radar gate
    n_moments = length(moment_dict)
    radar_grid = fill(-32768.0,n_moments,size(gridpoints,1),size(gridpoints,2))
    weights = zeros(Float64,n_moments,size(gridpoints,1),size(gridpoints,2))

    # Allocate a regular latlon grid for the map
    latlon_grid = Array{Float64}(undef,size(gridpoints,1),size(gridpoints,2),2)

    # Loop through the horizontal indices then do each column
    Threads.@threads for i in CartesianIndices(size(gridpoints)[1:2])

        # Calculate the lat, lon of the gridpoint
        yx_point = gridpoints[i,1:2]

        # Using CoordRefSystems pure Julia transform, the Proj.jl wrapper was significantly slower for unknown reasons
        cartTM = convert(TM,Cartesian{WGS84Latest}(grid_origin.x + yx_point[2]u"m", grid_origin.y + yx_point[1]u"m"))
        latlon = convert(LatLon,cartTM)
        latlon_grid[i,1] = ustrip(latlon.lat)
        latlon_grid[i,2] = ustrip(latlon.lon)
        
        # Find the points within the radius of influence of the horizontal gridpoint
        eff_h_radius_influence = horizontal_roi
        if beam_inflation > 0.0
            origin_dist = euclidean(yx_point, [0.0, 0.0])
            eff_h_radius_influence = max(beam_inflation * origin_dist, horizontal_roi)
        end
        gates = inrange(balltree, yx_point, eff_h_radius_influence)

        if !isempty(gates)

            # Found some gates that are within range horizontally
            valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[missing_key]])))
            if !isempty(valid_gates)
                # There is at least one gate in range so set the flags to -9999
                for m in 1:n_moments
                    radar_grid[m,i] = -9999.0
                end
            else
                continue
            end

            # Loops through the nearby gates with valid data
            valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[valid_key]])))
            for gate in gates[valid_gates]

                # Calculate the effective azimuth
                dx = yx_point[2] - radar_zyx[gate][3]
                dy = yx_point[1] - radar_zyx[gate][2]
                gridpt_az = (pi/2.0) - atan(dy, dx)
                if gridpt_az < 0
                    gridpt_az += 2*pi
                end

                # Calculate the spherical angle difference using Haversine formula
                angle_diff = spherical_angle([beams[gate,1], beams[gate,2]], 
                    [gridpt_az, beams[gate,2]])

                # Use the half-power beamwidth to define the beam
                # 79.43 = ln(0.5) / (0.5 deg * pi / 180.0)
                # Center of beam = 1.0 angle_weight
                angle_weight = exp(-angle_diff * 79.43)
                if angle_weight < power_threshold
                    angle_weight = 0
                end

                # Range weighting based on center of gridbox = 1.0 range_weight
                r = beams[gate,3]
                gridpt_r = sqrt(dx^2 + dy^2)
                range_weight = gridpt_r / r

                # Multiply weights so that center of beam is 1.0
                total_weight = range_weight * angle_weight

                # If there is non-zero weight, add it to the grid box
                if total_weight > 0.0
                    for m in 1:n_moments
                        if weights[m,i] == 0.0
                            # Initialize the radar grid box with 0 since there is a possibility that the beam hit it
                            radar_grid[m,i] = 0.0
                        end

                        if !ismissing(radar_volume.moments[gate,m])
                            if grid_type_dict[m] == :linear
                                linear_z = 10.0 ^ (radar_volume.moments[gate,m] / 10.0)
                                radar_grid[m,i] += total_weight * linear_z
                                weights[m,i] += total_weight
                            elseif grid_type_dict[m] == :nearest
                                if total_weight > weights[m,i]
                                    radar_grid[m,i] = radar_volume.moments[gate,m]
                                    weights[m,i] = total_weight
                                end
                            else
                                radar_grid[m,i] += total_weight * radar_volume.moments[gate,m]
                                weights[m,i] += total_weight
                            end
                        end
                    end
                end

            end # End of gate loop

            # Divide by the total weight for that gridbox
            for m in 1:n_moments
                if weights[m,i] > 0.0 && grid_type_dict[m] != :nearest
                    radar_grid[m,i] /= weights[m,i]
                    if grid_type_dict[m] == :linear
                        if radar_grid[m,i] > 0.0
                            radar_grid[m,i] = 10.0 * log10(radar_grid[m,i])
                        else
                            radar_grid[m,i] = -9999.0
                        end
                    end
                elseif weights[m,i] == 0.0
                    if radar_grid[m,i] == 0.0
                        # Use a different flag to indicate data has been QCed out
                        radar_grid[m,i] = -9999.0
                    end
                end
            end # End of moment loop
        end # End of gate empty test
    end # End of horizontal loop

    # Return the gridded radar array
    return radar_grid, latlon_grid
end

function grid_composite(reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, gridpoints::AbstractArray, 
        radar_volume::radar, moment_dict::Dict, grid_type_dict::Dict, horizontal_roi::Float64, beam_inflation::Float64,
        missing_key::String="SQI", valid_key::String="DBZ")
    
    # Convert the relevant radar information to arrays
    TM = CoordRefSystems.shift(TransverseMercator{1.0,reference_latitude,WGS84Latest}, lonₒ= reference_longitude)
    grid_origin, radar_zyx, beams = radar_arrays(reference_latitude, reference_longitude, radar_volume, TM)
    
    # Create a balltree that has horizontal locations of every gate in Y, X dimension
    balltree = radar_balltree_yx(radar_volume, radar_zyx, beams)

    # Allocate the grid for the radar moments and the weights for each radar gate
    n_moments = length(moment_dict)
    radar_grid = fill(-32768.0,n_moments,size(gridpoints,1),size(gridpoints,2))
    weights = zeros(Float64,n_moments,size(gridpoints,1),size(gridpoints,2))

    # Allocate a regular latlon grid for the map
    latlon_grid = Array{Float64}(undef,size(gridpoints,1),size(gridpoints,2),2)

    # Loop through the horizontal indices then do each column
    Threads.@threads for i in CartesianIndices(size(gridpoints)[1:2])

        # Calculate the lat, lon of the gridpoint
        yx_point = gridpoints[i,1:2]

        # Using CoordRefSystems pure Julia transform, the Proj.jl wrapper was significantly slower for unknown reasons
        cartTM = convert(TM,Cartesian{WGS84Latest}(grid_origin.x + yx_point[2]u"m", grid_origin.y + yx_point[1]u"m"))
        latlon = convert(LatLon,cartTM)
        latlon_grid[i,1] = ustrip(latlon.lat)
        latlon_grid[i,2] = ustrip(latlon.lon)
        
        # Find the points within the radius of influence of the horizontal gridpoint
        eff_h_radius_influence = horizontal_roi
        if beam_inflation > 0.0
            origin_dist = euclidean(yx_point, [0.0, 0.0])
            eff_h_radius_influence = max(beam_inflation * origin_dist, horizontal_roi) #0.01745 = 1 deg beam
        end
        gates = inrange(balltree, yx_point, eff_h_radius_influence)

        if !isempty(gates)

            # Found some gates that are within range horizontally
            valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[missing_key]])))
            if !isempty(valid_gates)
                # There is at least one gate in range so set the flags to -9999
                for m in 1:n_moments
                    radar_grid[m,i] = -9999.0
                end
            else
                continue
            end

            # Loops through the nearby gates with valid data
            valid_gates = collect(keys(skipmissing(radar_volume.moments[gates,moment_dict[valid_key]])))
            if !isempty(valid_gates)
                dbzmax, max_idx = findmax(gate -> radar_volume.moments[gate,moment_dict[valid_key]], gates[valid_gates])
    
                # Divide by the total weight for that gridbox
                for m in 1:n_moments
                    if !ismissing(radar_volume.moments[gates[valid_gates][max_idx],m])
                        radar_grid[m,i] = radar_volume.moments[gates[valid_gates][max_idx],m]
                    end
                end # End of moment loop
            end # End of valid gates test
        end # End of gate empty test
    end # End of horizontal loop

    # Return the gridded radar array
    return radar_grid, latlon_grid
end

function write_gridded_radar_volume(file, start_time, stop_time, gridpoints, radar_grid, latlon_grid, moment_dict, reference_latitude::AbstractFloat, reference_longitude::AbstractFloat)

    # Delete any pre-existing file
    rm(file, force=true)
    
    ds = NCDataset(file,"c", attrib = OrderedDict(
        "Conventions"               => "CF-1.6",
        "history"                   => "Created by Michael M. Bell using custom software",
        "institution"               => "CSU",
        "source"                    => "SEAPOL",
        "title"                     => "PRELIMINARY Gridded Radar Data",
        "comment"                   => "PRELIMINARY In-field Analysis. Please use with caution!",
    ))
    
    # Dimensions
    # Could concatenate multiple volumes here
    #numswps = length(swpstart)
    xdim = size(radar_grid,4)
    ydim = size(radar_grid,3)
    zdim = size(radar_grid,2)
    ds.dim["time"] = 1
    ds.dim["X"] = xdim
    ds.dim["Y"] = ydim
    ds.dim["Z"] = zdim
    
    # Declare variables
    
    nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "time",
        "long_name"                 => "Data time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "axis"                      => "T",
        "comment"                   => "",
    ))

    ncstarttime = defVar(ds,"start_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "start_time",
        "long_name"                 => "Data start time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))

    ncstoptime = defVar(ds,"stop_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "stop_time",
        "long_name"                 => "Data stop time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))
    
    ncx = defVar(ds,"X", Float32, ("X",), attrib = OrderedDict(
        "standard_name"             => "projection_x_coordinate",
        "units"                     => "km",
        "axis"                      => "X",
    ))

    ncy = defVar(ds,"Y", Float32, ("Y",), attrib = OrderedDict(
        "standard_name"             => "projection_y_coordinate",
        "units"                     => "km",
        "axis"                      => "Y",
    ))

    ncz = defVar(ds,"Z", Float32, ("Z",), attrib = OrderedDict(
        "standard_name"             => "altitude",
        "long_name"                 => "constant altitude levels",
        "units"                     => "km",
        "positive"                  => "up",
        "axis"                      => "Z",
    ))

    nclat = defVar(ds,"latitude", Float32, ("X", "Y", "time"), attrib = OrderedDict(
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))
    
    nclon = defVar(ds,"longitude", Float32, ("X", "Y", "time"), attrib = OrderedDict(
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    ncgrid_mapping = defVar(ds,"grid_mapping", Int32, (), attrib = OrderedDict(
        "grid_mapping_name"         => "tranverse_mercator",
        "scale_factor_at_central_meridian" => 1.0,
        "longitude_of_central_meridian" => reference_longitude,
        "latitude_of_projection_origin" => reference_latitude,
        "reference_ellipsoid_name"  => "GRS80",
        "false_easting"             => 0.0,
        "false_northing"            => 0.0,        
    ))

    # Using start time for now, but eventually need to use some reference time
    nctime[:] = datetime2unix.(stop_time)
    ncstarttime[:] = datetime2unix.(start_time)
    ncstoptime[:] = datetime2unix.(stop_time)
    ncz[:] = gridpoints[:,1,1,1] 
    ncy[:] = gridpoints[1,:,1,2] 
    ncx[:] = gridpoints[1,1,:,3] 
    nclat[:] = latlon_grid[:,:,1]' 
    nclon[:] = latlon_grid[:,:,2]' 
    ncgrid_mapping[:] = -32768.0

    # Define variables
    perm = (1, 4, 3, 2)
    # moment, z, y, x -> moment, x, y, z
    ncgrid = permutedims(radar_grid,perm)
    
    # Loop through the moments
    for key in keys(moment_dict)
        var_attrib = merge(common_attrib, variable_attrib_dict[key])
        ncvar = defVar(ds, key, Float32, ("X", "Y", "Z", "time"), attrib = var_attrib)
        ncvar[:] = ncgrid[moment_dict[key],:,:,:]
    end

    close(ds)
end

function write_gridded_radar_rhi(file, start_time, stop_time, gridpoints, radar_grid, latlon_grid, moment_dict, reference_latitude::AbstractFloat, reference_longitude::AbstractFloat)

    # Delete any pre-existing file
    rm(file, force=true)
    
    ds = NCDataset(file,"c", attrib = OrderedDict(
        "Conventions"               => "CF-1.6",
        "history"                   => "Created by Michael M. Bell using custom software",
        "institution"               => "CSU",
        "source"                    => "SEAPOL",
        "title"                     => "PRELIMINARY Gridded Radar Data",
        "comment"                   => "PRELIMINARY In-field Analysis. Please use with caution!",
    ))
    
    # Dimensions
    # Could concatenate multiple volumes here
    #numswps = length(swpstart)
    rdim = size(radar_grid,3)
    zdim = size(radar_grid,2)
    ds.dim["time"] = 1
    ds.dim["R"] = rdim
    ds.dim["Z"] = zdim
    
    # Declare variables
    
    nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "time",
        "long_name"                 => "Data time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "axis"                      => "T",
        "comment"                   => "",
    ))

    ncstarttime = defVar(ds,"start_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "start_time",
        "long_name"                 => "Data start time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))

    ncstoptime = defVar(ds,"stop_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "stop_time",
        "long_name"                 => "Data stop time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))
    
    ncr = defVar(ds,"R", Float32, ("R",), attrib = OrderedDict(
        "standard_name"             => "projection_range_coordinate",
        "units"                     => "km",
        "axis"                      => "R",
    ))

    ncz = defVar(ds,"Z", Float32, ("Z",), attrib = OrderedDict(
        "standard_name"             => "altitude",
        "long_name"                 => "constant altitude levels",
        "units"                     => "km",
        "positive"                  => "up",
        "axis"                      => "Z",
    ))

    nclat = defVar(ds,"latitude", Float32, ("R", "time"), attrib = OrderedDict(
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))
    
    nclon = defVar(ds,"longitude", Float32, ("R", "time"), attrib = OrderedDict(
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    ncgrid_mapping = defVar(ds,"grid_mapping", Int32, (), attrib = OrderedDict(
        "grid_mapping_name"         => "tranverse_mercator",
        "scale_factor_at_central_meridian" => 1.0,
        "longitude_of_central_meridian" => reference_longitude,
        "latitude_of_projection_origin" => reference_latitude,
        "reference_ellipsoid_name"  => "GRS80",
        "false_easting"             => 0.0,
        "false_northing"            => 0.0,        
    ))

    # Using start time for now, but eventually need to use some reference time
    nctime[:] = datetime2unix.(stop_time)
    ncstarttime[:] = datetime2unix.(start_time)
    ncstoptime[:] = datetime2unix.(stop_time)
    ncz[:] = gridpoints[:,1,1] 
    ncr[:] = gridpoints[1,:,2]
    nclat[:] = latlon_grid[:,1]' 
    nclon[:] = latlon_grid[:,2]'
    ncgrid_mapping[:] = -32768.0

    # Define variables
    perm = (1, 3, 2)
    # moment, z, r -> moment, r, z
    ncgrid = permutedims(radar_grid,perm)
    
    # Loop through the moments
    for key in keys(moment_dict)
        var_attrib = common_attrib
        if haskey(variable_attrib_dict,key)
            var_attrib = merge(common_attrib, variable_attrib_dict[key])
        else
            var_attrib = merge(common_attrib, variable_attrib_dict["UNKNOWN"])
        end
        ncvar = defVar(ds, key, Float32, ("R", "Z", "time"), attrib = var_attrib)
        ncvar[:] = ncgrid[moment_dict[key],:,:]
    end

    close(ds)
end

function write_gridded_radar_ppi(file, start_time, stop_time, gridpoints, radar_grid, latlon_grid, moment_dict, reference_latitude::AbstractFloat, reference_longitude::AbstractFloat, mean_heading::AbstractFloat)

    # Delete any pre-existing file
    rm(file, force=true)
    
    ds = NCDataset(file,"c", attrib = OrderedDict(
        "Conventions"               => "CF-1.6",
        "history"                   => "Created by Michael M. Bell using custom software",
        "institution"               => "CSU",
        "source"                    => "SEAPOL",
        "title"                     => "Gridded Radar Data v1.0",
        "comment"                   => "Level 4 Gridded Radar Data"
    ))
    
    # Dimensions
    # Could concatenate multiple volumes here
    #numswps = length(swpstart)
    xdim = size(radar_grid,3)
    ydim = size(radar_grid,2)
    ds.dim["time"] = 1
    ds.dim["X"] = xdim
    ds.dim["Y"] = ydim
    
    # Declare variables
    
    nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "time",
        "long_name"                 => "Data time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "axis"                      => "T",
        "comment"                   => "",
    ))

    ncstarttime = defVar(ds,"start_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "start_time",
        "long_name"                 => "Data start time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))

    ncstoptime = defVar(ds,"stop_time", Float64, ("time",), attrib = OrderedDict(
        "standard_name"             => "stop_time",
        "long_name"                 => "Data stop time",
        "units"                     => "seconds since 1970-01-01T00:00:00Z",
        "comment"                   => "",
    ))
    
    ncx = defVar(ds,"X", Float32, ("X",), attrib = OrderedDict(
        "standard_name"             => "projection_x_coordinate",
        "units"                     => "km",
        "axis"                      => "X",
    ))

    ncy = defVar(ds,"Y", Float32, ("Y",), attrib = OrderedDict(
        "standard_name"             => "projection_y_coordinate",
        "units"                     => "km",
        "axis"                      => "Y",
    ))


    nclat = defVar(ds,"latitude", Float32, ("X", "Y", "time"), attrib = OrderedDict(
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))
    
    nclon = defVar(ds,"longitude", Float32, ("X", "Y", "time"), attrib = OrderedDict(
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    ncgrid_mapping = defVar(ds,"grid_mapping", Int32, (), attrib = OrderedDict(
        "grid_mapping_name"         => "tranverse_mercator",
        "scale_factor_at_central_meridian" => 1.0,
        "longitude_of_central_meridian" => reference_longitude,
        "latitude_of_projection_origin" => reference_latitude,
        "reference_ellipsoid_name"  => "GRS80",
        "false_easting"             => 0.0,
        "false_northing"            => 0.0,        
    ))

    ncheading = defVar(ds, "heading", Float32, ("time",), attrib = OrderedDict(
        "standard_name"             => "heading",
        "units"                     => "degrees",
    ))

    # Using start time for now, but eventually need to use some reference time
    nctime[:] = datetime2unix.(stop_time)
    ncstarttime[:] = datetime2unix.(start_time)
    ncstoptime[:] = datetime2unix.(stop_time)
    ncy[:] = gridpoints[:,1,1] 
    ncx[:] = gridpoints[1,:,2]
    nclat[:] = latlon_grid[:,:,1]' 
    nclon[:] = latlon_grid[:,:,2]' 
    ncgrid_mapping[:] = -32768.0
    ncheading[:] = mean_heading

    # Define variables
    perm = (1, 3, 2)
    # moment, y, x -> moment, x, y
    ncgrid = permutedims(radar_grid,perm)
    
    # Loop through the moments
    for key in keys(moment_dict)
        var_attrib = common_attrib
        if haskey(variable_attrib_dict,key)
            var_attrib = merge(common_attrib, variable_attrib_dict[key])
        else
            var_attrib = merge(common_attrib, variable_attrib_dict["UNKNOWN"])
        end
        ncvar = defVar(ds, key, Float32, ("X", "Y", "time"), attrib = var_attrib)
        ncvar[:] = ncgrid[moment_dict[key],:,:]
    end

    close(ds)
end

function read_latlon_gridded_radar(file, moment_dict)

    inputds = Dataset(file);

    x0 = inputds["x0"]
    y0 = inputds["y0"]
    z0 = inputds["z0"]
    start_time = inputds["start_time"]
    stop_time = inputds["stop_time"]
    
    # Store radar data
    n_moments = length(moment_dict)
    n_points = length(x0)*length(y0)*length(z0)
    radardata = Array{Union{Missing, Float32}}(undef,n_moments,n_points)
    for key in keys(moment_dict)
        radardata[moment_dict[key],:] = inputds[key][:]
    end
    
    return x0, y0, z0, start_time, stop_time, radardata
end

function read_cartesian_gridded_radar(file, moment_dict)

    inputds = Dataset(file);

    x0 = inputds["x0"]
    y0 = inputds["y0"]
    z0 = inputds["z0"]
    lon0 = inputds["lon0"]    
    lat0 = inputds["lat0"]
    start_time = inputds["start_time"]
    stop_time = inputds["stop_time"]
    
    # Store radar data
    n_moments = length(moment_dict)
    n_points = length(x0)*length(y0)*length(z0)
    radardata = Array{Union{Missing, Float32}}(undef,n_moments,n_points)
    for key in keys(moment_dict)
        radardata[moment_dict[key],:] = inputds[key][:]
    end

    return x0, y0, z0, lat0, lon0, start_time, stop_time, radardata
end

function read_gridded_radar(file, moment_dict)

    inputds = Dataset(file);

    x = inputds["X"]
    y = inputds["Y"]
    z = inputds["Z"]
    lon = inputds["longitude"]    
    lat = inputds["latitude"]
    start_time = inputds["start_time"]
    stop_time = inputds["stop_time"]
    
    # Store radar data
    n_moments = length(moment_dict)
    n_points = length(x)*length(y)*length(z)
    radardata = Array{Union{Missing, Float32}}(undef,n_moments,n_points)
    for key in keys(moment_dict)
        radardata[moment_dict[key],:] = inputds[key][:]
    end

    return x, y, z, lat, lon, start_time, stop_time, radardata
end

function read_gridded_ppi(file, moment_dict)

    inputds = Dataset(file);

    x = inputds["X"]
    y = inputds["Y"]
    lon = inputds["longitude"]    
    lat = inputds["latitude"]
    start_time = inputds["start_time"]
    stop_time = inputds["stop_time"]
    
    # Store radar data
    n_moments = length(moment_dict)
    n_points = length(x)*length(y)
    radardata = Array{Union{Missing, Float32}}(undef,n_moments,n_points)
    for key in keys(moment_dict)
        radardata[moment_dict[key],:] = inputds[key][:]
    end

    return x, y, lat, lon, start_time, stop_time, radardata
end

function read_gridded_rhi(file, moment_dict)

    inputds = Dataset(file);

    R = inputds["R"]
    Z = inputds["Z"]
    lon = inputds["longitude"]    
    lat = inputds["latitude"]
    start_time = inputds["start_time"]
    stop_time = inputds["stop_time"]
    
    # Store radar data
    n_moments = length(moment_dict)
    n_points = length(R)*length(Z)
    radardata = Array{Union{Missing, Float32}}(undef,n_moments,n_points)
    for key in keys(moment_dict)
        radardata[moment_dict[key],:] = inputds[key][:]
    end

    return R, Z, lat, lon, start_time, stop_time, radardata
end