"""
liquid_water_path(vec, config)

Returns liquid water path [kg/m^2]

# Inputs:
- `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
- `config` :: dictionary containing settings
"""
function liquid_water_path(vec::Vector{FT}, config) where {FT <: Real}

    @assert issubset(Set(["rho", "ql"]), Set(config["observations"]["data_names"]))

    (; n_cases, n_heights, n_times) = get_numbers_from_config(config)

    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_ql = findall(x -> x == "ql", config["observations"]["data_names"])[1]
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_heights[ind_ρ]
    n_single_case = sum(n_heights) * n_times

    lwp = []
    for i in 1:n_cases
        _v = get_case_i_vec(vec, i, n_single_case)
        _fields = get_single_case_fields(_v, n_heights, n_times[i])
        ρ = _fields[ind_ρ]
        ql = _fields[ind_ql]
        push!(lwp, sum(ql .* ρ, dims = 1) .* dz)
    end
    return lwp
end


"""
rainwater_path(vec, config)

Returns rain water path [kg/m^2]

# Inputs:
- `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
- `config` :: dictionary containing settings
"""
function rainwater_path(vec::Vector{FT}, config) where {FT <: Real}

    @assert issubset(Set(["rho", "qr"]), Set(config["observations"]["data_names"]))

    (; n_cases, n_heights, n_times) = get_numbers_from_config(config)

    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_qr = findall(x -> x == "qr", config["observations"]["data_names"])[1]
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_heights[ind_ρ]
    n_single_case = sum(n_heights) * n_times

    rwp = []
    for i in 1:n_cases
        _v = get_case_i_vec(vec, i, n_single_case)
        _fields = get_single_case_fields(_v, n_heights, n_times[i])
        ρ = _fields[ind_ρ]
        qr = _fields[ind_qr]
        push!(rwp, sum(qr .* ρ, dims = 1) .* dz)
    end
    return rwp
end

"""
rainrate(vec, config; height)

Returns rain rate [mm/hr]

# Inputs:
    - `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
    - `config` :: dictionary containing settings
    - `height` :: height at which rain rate is computed. By default, the rain rate is computed at 
        the ground level `height = 0.0`
"""
function rainrate(vec::Vector{FT}, config; height::FT = FT(0)) where {FT <: Real}

    @assert issubset(Set(["rainrate"]), Set(config["observations"]["data_names"]))

    (; n_cases, n_heights, n_times) = get_numbers_from_config(config)

    ind_rr = findall(x -> x == "rainrate", config["observations"]["data_names"])[1]
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_heights[ind_rr]
    n_single_case = sum(n_heights) * n_times
    ind_z1 = min(n_heights[ind_rr] - 1, max(1, ceil(Int, (height - config["model"]["z_min"]) / dz - FT(0.5))))
    ind_z2 = ind_z1 + 1

    a1 = ((ind_z2 - 0.5) * dz - height) / dz
    a2 = (height - (ind_z1 - 0.5) * dz) / dz

    rainrate = []
    for i in 1:n_cases
        _v = get_case_i_vec(vec, i, n_single_case)
        _fields = get_single_case_fields(_v, n_heights, n_times[i])
        _rainrate = a1 .* _fields[ind_rr][ind_z1, :] + a2 .* _fields[ind_rr][ind_z2, :]
        _rainrate[findall(x -> x < 0, _rainrate)] .= FT(0)
        push!(rainrate, _rainrate)
    end

    return rainrate
end
