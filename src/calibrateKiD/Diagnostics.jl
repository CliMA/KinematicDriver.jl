"""
lwp(vec, config)

Returns liquid water path [kg/m^2]

# Inputs:
- `vec` :: vector containing field information; 
       must have the same structure as vectors in calibrations (height, variable, time, case)
- `config` :: dictionary containing settings
"""
function lwp(vec::Vector{FT}, config) where {FT <: Real}

    @assert issubset(Set(["rho", "ql"]), Set(config["observations"]["data_names"]))

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    lwp = zeros(n_t, n_c)
    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_ql = findall(x -> x == "ql", config["observations"]["data_names"])[1]

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[((ind_ρ - 1) * n_z + 1):(ind_ρ * n_z), :]
        ql = m_[((ind_ql - 1) * n_z + 1):(ind_ql * n_z), :]
        lwp[:, i] = sum(ql .* ρ, dims = 1) .* dz
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

    @assert issubset(Set(["rho", "qr"]), Set(config["observations"]["data_names"]))

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    rwp = zeros(n_t, n_c)
    ind_ρ = findall(x -> x == "rho", config["observations"]["data_names"])[1]
    ind_qr = findall(x -> x == "qr", config["observations"]["data_names"])[1]

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[((ind_ρ - 1) * n_z + 1):(ind_ρ * n_z), :]
        qr = m_[((ind_qr - 1) * n_z + 1):(ind_qr * n_z), :]
        rwp[:, i] = sum(qr .* ρ, dims = 1) .* dz
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

    (n_c, n_v, n_z, n_t) = get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    rainrate = zeros(n_t, n_c)
    ind_rr = findall(x -> x == "rainrate", config["observations"]["data_names"])[1]
    ind_z1 = min(n_z - 1, max(1, ceil(Int, (height - config["model"]["z_min"]) / dz - FT(0.5))))
    ind_z2 = ind_z1 + 1

    a1 = ((ind_z2 - 0.5) * dz - height) / dz
    a2 = (height - (ind_z1 - 0.5) * dz) / dz

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        rainrate_ = m_[((ind_rr - 1) * n_z + 1):(ind_rr * n_z), :]
        rainrate[:, i] = a1 .* rainrate_[ind_z1, :] + a2 .* rainrate_[ind_z2, :]
    end
    rainrate[findall(x -> x < 0, rainrate)] .= Float64(0)
    return rainrate
end
