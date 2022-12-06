
function get_obs!(config::Dict)

    FT = config["observations"]["data_type"]
    _n_heights =
        config["model"]["filter"]["apply"] ? config["model"]["filter"]["nz_filtered"] : config["model"]["n_elem"]
    _dz = (config["model"]["z_max"] - config["model"]["z_min"]) / _n_heights
    _heights::Array{FT} =
        collect(range(config["model"]["z_min"] + _dz / 2, config["model"]["z_max"] - _dz / 2, _n_heights))
    _times::Array{FT} = collect(config["model"]["t_calib"])
    _variables::Array{String} = config["observations"]["data_names"]

    if config["observations"]["data_source"] == "file"
        _cases::Vector{NamedTuple{(:w1, :p0, :Nd, :dir), Tuple{Float64, Float64, Float64, String}}} =
            config["observations"]["cases"]
        _y::Matrix{FT} =
            get_obs_matrix(_cases, _variables, _heights, _times; apply_filter = config["model"]["filter"]["apply"])
    elseif config["observations"]["data_source"] == "perfect_model"
        _y = get_validation_samples(config)
    else
        @error("Invalid data source!")
    end

    # Normalize variables by their maximum (times given numbers) for each case separately 
    _ynorm = config["observations"]["ynorm"]
    _n_times = config["model"]["filter"]["apply"] ? length(_times) - 1 : length(_times)
    normalize_y!(_y, _ynorm, length(_variables), _n_heights, _n_times)
    config["observations"]["ynorm"] = _ynorm

    _Γy::Matrix{FT} = cov(_y, dims = 2, corrected = false)

    return Observations.Observation(_y, _Γy, _variables)

end

# fill norm_vec with maximum of mean profiles (for all variables) in space and time
function normalize_y!(
    y::Matrix{FT},
    ynorm::Matrix{FT},
    n_variables::Int,
    n_heights::Int,
    n_times::Int,
) where {FT <: Real}

    _n_cases = size(ynorm)[1]
    _n_single_sim = n_variables * n_heights * n_times
    @assert size(y)[1] == (_n_cases * _n_single_sim)

    for i in 1:_n_cases
        _y_single_case = y[((i - 1) * _n_single_sim + 1):(i * _n_single_sim), :]

        _y_mean_vector_single_case = mean(_y_single_case, dims = 2)
        _y_mean_matrix::Matrix{FT} = reshape(_y_mean_vector_single_case, n_variables * n_heights, n_times)
        _var_max = [maximum(abs.(_y_mean_matrix[((j - 1) * n_heights + 1):(j * n_heights), :])) for j in 1:n_variables]
        _var_max = ifelse.(_var_max .> eps(FT), _var_max, 1.0)
        _norm_vec = ynorm[i, :] .* _var_max

        # normalise part of y belonging to case i by _norm_vec
        _norm_vec_extended_single_time = reshape((ones(n_heights, n_variables) .* _norm_vec'), :, 1)
        _norm_vec_extended = reshape(ones(n_heights * n_variables, n_times) .* _norm_vec_extended_single_time, :, 1)
        _y_single_case = _y_single_case ./ _norm_vec_extended
        y[((i - 1) * _n_single_sim + 1):(i * _n_single_sim), :] = _y_single_case
        ynorm[i, :] = _norm_vec
    end

end

function get_obs_matrix(
    cases::Vector{NamedTuple{(:w1, :p0, :Nd, :dir), Tuple{Float64, Float64, Float64, String}}},
    variables::Array{String},
    heights::Array{FT},
    times::Array{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _data_matrix = Matrix{FT}
    for (i, case) in enumerate(cases)
        _dir::String = case.dir
        _data_matrix_single_case::Matrix{FT} =
            get_obs_matrix(_dir, variables, heights, times; apply_filter = apply_filter)

        if i == 1
            _data_matrix = _data_matrix_single_case
        else
            cols = minimum([size(_data_matrix)[2], size(_data_matrix_single_case)[2]])
            _data_matrix = [_data_matrix[:, 1:cols]; _data_matrix_single_case[:, 1:cols]]
        end
    end

    return _data_matrix

end

function get_obs_matrix(
    dir::String,
    variables::Array{String},
    heights::Array{FT},
    times::Array{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _n_times = apply_filter ? length(times) - 1 : length(times)
    _n_single_data_vector_size::Int = length(heights) * _n_times * length(variables)

    _data_matrix = zeros(FT, _n_single_data_vector_size, 0)
    foreach(readdir(dir)) do file
        _filename = dir * file
        _file_data_vector = get_single_obs_vector(_filename, variables, heights, times; apply_filter = apply_filter)
        _data_matrix = hcat(_data_matrix, _file_data_vector)
    end

    return _data_matrix

end

function get_single_obs_vector(
    filename::String,
    variables::Array{String},
    heights::Array{FT},
    times::Array{FT};
    apply_filter::Bool = false,
) where {FT <: Real}
    _data_dict = Dict()
    _data_dict = get_single_obs_field(filename, variables, heights, times; apply_filter = apply_filter)

    _n_times = apply_filter ? length(times) - 1 : length(times)
    _outputs = FT[]
    for i in 1:_n_times
        _single_time_output = FT[]
        for var in variables
            _single_time_output = [_single_time_output; _data_dict[var][:, i]]
        end
        _outputs = [_outputs; _single_time_output]
    end

    return _outputs
end

function get_single_obs_field(
    filename::String,
    variables::Array{String},
    heights::Array{FT},
    times::Array{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _data_pysdm = read_pysdm_data(filename).variables
    _z_data = _data_pysdm["height"]
    _t_data = _data_pysdm["time"]

    _output = Dict()

    for var in variables
        if var == "qt"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = _r_tot ./ (1 .+ _r_tot)
        elseif var == "qv"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = _data_pysdm["qv"] ./ (1 .+ _r_tot)
        elseif var == "ql"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = _data_pysdm["ql"] .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qr"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = _data_pysdm["qr"] .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qlr"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = (_data_pysdm["ql"] .+ _data_pysdm["qr"]) .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qtr"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = (_r_tot .+ _data_pysdm["qr"] .* 1e-3) ./ (1 .+ _r_tot)
        elseif var == "rho"
            _r_tot = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
            _data = _data_pysdm["rhod"] .* (1 .+ _r_tot)
        elseif var == "rt"
            _data = _data_pysdm["qv"] .+ _data_pysdm["ql"] .* 1e-3
        elseif var == "rv"
            _data = _data_pysdm["qv"]
        elseif var == "rl"
            _data = _data_pysdm["ql"] .* 1e-3
        elseif var == "rr"
            _data = _data_pysdm["qr"] .* 1e-3
        elseif var == "rlr"
            _data = (_data_pysdm["ql"] .+ _data_pysdm["qr"]) .* 1e-3
        elseif var == "rtr"
            _data = _data_pysdm["qv"] .+ (_data_pysdm["ql"] .+ _data_pysdm["qr"]) .* 1e-3
        else
            _data = _data_pysdm[var]
        end

        _f = linear_interpolation((_t_data, _z_data), _data, extrapolation_bc = Line())
        if apply_filter
            _dz_data = _z_data[2] - _z_data[1]
            _n_z_data = length(_z_data)
            _dz_heights = heights[2] - heights[1]
            _n_heights = length(heights)
            _t_data_c = [0.5 * (_t_data[i] + _t_data[i + 1]) for i in 1:(length(_t_data) - 1)]
            _data_c = [_f(t, z) for z in _z_data, t in _t_data_c]
            _heights_f = collect(range(heights[1] - _dz_heights / 2, heights[end] + _dz_heights / 2, _n_heights + 1))
            _z_data_f = collect(range(_z_data[1] - _dz_data / 2, _z_data[end] + _dz_data / 2, _n_z_data + 1))
            _output[var] = filter_field(_data_c, _t_data, times, _z_data_f, _heights_f)
        else
            _output[var] = [_f(t, z) for z in heights, t in times]
        end
    end

    return _output
end

function filter_field(data::Array{FT}, t_data::Array, times::Array{FT}, z_data::Array, heights::Array{FT})
    @assert issubset(Set(times), Set(t_data))
    @assert issubset(Set(heights), Set(z_data))
    @assert size(data) == (length(z_data) - 1, length(t_data) - 1)

    _nz_calib = length(heights) - 1
    _nt_calib = length(times) - 1
    _output = zeros(_nz_calib, _nt_calib)
    for i in 1:_nz_calib
        ind_ini_z = findall(x -> x == heights[i], z_data)[1]
        ind_end_z = findall(x -> x == heights[i + 1], z_data)[1] - 1
        for j in 1:_nt_calib
            ind_ini_t = findall(x -> x == times[j], t_data)[1]
            ind_end_t = findall(x -> x == times[j + 1], t_data)[1] - 1
            _output[i, j] = mean(data[ind_ini_z:ind_end_z, ind_ini_t:ind_end_t])
        end
    end
    return _output
end

function get_validation_samples(config)
    _ϕ_true = params_validation(config).ϕ

    # Generate reference sample
    _ϕ_names = collect(keys(config["prior"]["parameters"]))
    _G_t = run_dyn_model(_ϕ_true, _ϕ_names, config)

    # Generate (artificial) truth samples
    _n_samples = config["observations"]["number_of_samples"]
    _samples = zeros(length(_G_t), _n_samples)

    # Since KiD is deterministic we add an artificial noise to G_t to generate a set of data with noise
    _Γy = convert(Array, Diagonal((_G_t .* _G_t .+ eps(Float64))[:]))
    _μ = zeros(length(_G_t))
    scov_G_ratio = config["observations"]["scov_G_ratio"]
    for i in 1:_n_samples
        _samples[:, i] = _G_t .+ rand(MvNormal(_μ, _Γy)) .* scov_G_ratio
    end

    return _samples
end

function params_validation(config, priors = nothing)
    ϕ_mean = [v.mean for v in values(config["prior"]["parameters"])]
    ϕ = ϕ_mean .* (1 - config["observations"]["true_values_offset"])
    θ = priors === nothing ? 0 : transform_constrained_to_unconstrained(priors, ϕ)
    return (ϕ = ϕ, θ = θ)
end
