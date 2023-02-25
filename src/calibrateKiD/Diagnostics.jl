"""
lwp(vec, config)

Returns liquid water path [kg/m^2]

# Inputs:
- `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
- `config` :: dictionary containing settings
"""
function lwp(vec::Vector{FT}, config) where {FT <: Real}

    @assert issubset(Set(["rho", "rl"]), Set(config["observations"]["data_names"]))

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    lwp = zeros(n_t, n_c)
    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_rl = findall(x -> x == "rl", config["observations"]["data_names"])[1]

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[((ind_ρ - 1) * n_z + 1):(ind_ρ * n_z), :]
        rl = m_[((ind_rl - 1) * n_z + 1):(ind_rl * n_z), :]
        lwp[:, i] = sum(rl .* ρ, dims = 1) .* dz
    end
    return lwp
end


"""
rwp(vec, config)

Returns rain water path [kg/m^2]

# Inputs:
- `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
- `config` :: dictionary containing settings
"""
function rwp(vec::Vector{FT}, config) where {FT <: Real}

    @assert issubset(Set(["rho", "rr"]), Set(config["observations"]["data_names"]))

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    rwp = zeros(n_t, n_c)
    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_rr = findall(x -> x == "rr", config["observations"]["data_names"])[1]

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[((ind_ρ - 1) * n_z + 1):(ind_ρ * n_z), :]
        rr = m_[((ind_rr - 1) * n_z + 1):(ind_rr * n_z), :]
        rwp[:, i] = sum(rr .* ρ, dims = 1) .* dz
    end
    return rwp
end

"""
rainrate(vec, config; height, isref)

Returns rain rate [mm/hr]

# Inputs:
    - `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
    - `config` :: dictionary containing settings
    - `height` :: height at which rain rate is computed. By default, the rain rate is computed at 
        the ground level `height = 0.0`
    - `isref` :: option defining whether terminal velocity data should be read from vector or should
        be computed. Velocity is computed by default.
"""
function rainrate(vec::Vector{FT}, config; height::FT = FT(0), isref = false) where {FT <: Real}

    @assert issubset(Set(["rho", "rr", "rain averaged terminal velocity"]), Set(config["observations"]["data_names"]))

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    rainrate = zeros(n_t, n_c)
    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_rr = findall(x -> x == "rr", config["observations"]["data_names"])[1]
    ind_vt = findall(x -> x == "rain averaged terminal velocity", config["observations"]["data_names"])[1]
    ind_z = min(n_z, max(1, ceil(Int, (height - config["model"]["z_min"]) / dz)))

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[((ind_ρ - 1) * n_z + 1):(ind_ρ * n_z), :]
        rr = m_[((ind_rr - 1) * n_z + 1):(ind_rr * n_z), :]

        ρ_z = ρ[ind_z, :]
        rr_z = rr[ind_z, :]
        if isref
            vt = m_[((ind_vt - 1) * n_z + 1):(ind_vt * n_z), :]
            vt_z = vt[ind_z, :]
        else
            ϕ_names = collect(keys(config["prior"]["parameters"]))
            ϕ_values = collect([v.mean for v in values(config["prior"]["parameters"])])
            params_calib = create_param_dict(ϕ_values, ϕ_names)
            params = create_parameter_set(Float64, config["model"], params_calib)
            microphys_params = KP.microphysics_params(params)
            vt_z = CM1.terminal_velocity.(microphys_params, CMT.RainType(), ρ_z, rr_z)
        end
        rainrate[:, i] = rr_z .* ρ_z .* vt_z .* 3600
    end
    return rainrate
end
