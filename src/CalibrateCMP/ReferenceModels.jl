
function get_obs!(config::Dict)

    FT = config["observations"]["data_type"]
    _variables::Vector{String} = config["observations"]["data_names"]
    _n_heights = get_numbers_from_config(config).n_heights
    if config["model"]["model"] == "Box"
        _z_min = FT(0)
        _z_max = FT(1)
    else
        _z_min = config["model"]["z_min"]
        _z_max = config["model"]["z_max"]
    end
    _dz = (_z_max - _z_min) ./ _n_heights
    _heights::Vector{Vector{FT}} = collect.(range.(_z_min .+ _dz ./ 2, _z_max .- _dz ./ 2, _n_heights))
    _times::Vector{FT} = collect(config["model"]["t_calib"])

    if config["observations"]["data_source"] == "file"
        _cases::Vector{<:Any} = config["observations"]["cases"]
        _obs::Matrix{FT} =
            get_obs_matrix(_cases, _variables, _heights, _times; apply_filter = config["model"]["filter"]["apply"])
    elseif config["observations"]["data_source"] == "perfect_model"
        _obs = get_validation_samples(config)
    else
        @error("Invalid data source!")
    end

    return _obs

end

function get_obs_matrix(
    cases::Vector{<:Any},
    variables::Vector{String},
    heights::Vector{Vector{FT}},
    times::Vector{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _data_matrix = Matrix{FT}
    for (i, case) in enumerate(cases)
        _dir::String = case.dir

        if "t_cal" in collect(keys(case))
            times = case.t_cal
        end
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
    variables::Vector{String},
    heights::Vector{Vector{FT}},
    times::Vector{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _n_times = apply_filter ? length(times) - 1 : length(times)
    _n_single_data_vector_size::Int = sum(length.(heights)) * _n_times

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
    variables::Vector{String},
    heights::Vector{Vector{FT}},
    times::Vector{FT};
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
    variables::Vector{String},
    heights::Vector{Vector{FT}},
    times::Vector{FT};
    apply_filter::Bool = false,
) where {FT <: Real}

    _output = Dict()

    _data_pysdm = read_pysdm_data(filename).variables
    _t_data = _data_pysdm["time"]

    if "height" in keys(_data_pysdm)
        _z_data = _data_pysdm["height"]
    else
        # if height is not provided in pysdm data then the simulation is 0D
        _z_data = [FT(0.5)]
    end
    if "qv" in keys(_data_pysdm)
        _rv = _data_pysdm["qv"]
    elseif "water_vapour_mixing_ratio" in keys(_data_pysdm)
        _rv = _data_pysdm["water_vapour_mixing_ratio"]
    else
        _rv = zeros(FT, size(_data_pysdm["qc"]))
    end
    if ~("rhod" in keys(_data_pysdm))
        attr = read_pysdm_data(filename).attributes
        if "rhod" in keys(attr)
            _data_pysdm["rhod"] = attr["rhod"]
        end
    end

    for (i, var) in enumerate(variables)
        if var == "qt"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = _r_tot ./ (1 .+ _r_tot)
        elseif var == "qv"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = _rv ./ (1 .+ _r_tot)
        elseif var == "ql"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = _data_pysdm["qc"] .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qr"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = _data_pysdm["qr"] .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qlr"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = (_data_pysdm["qc"] .+ _data_pysdm["qr"]) .* 1e-3 ./ (1 .+ _r_tot)
        elseif var == "qtr"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = (_r_tot .+ _data_pysdm["qr"] .* 1e-3) ./ (1 .+ _r_tot)
        elseif var == "rho"
            _r_tot = _rv .+ _data_pysdm["qc"] .* 1e-3
            _data = _data_pysdm["rhod"] .* (1 .+ _r_tot)
        elseif var == "rt"
            _data = _rv .+ _data_pysdm["qc"] .* 1e-3
        elseif var == "rv"
            _data = _rv
        elseif var == "rl"
            _data = _data_pysdm["qc"] .* 1e-3
        elseif var == "rr"
            _data = _data_pysdm["qr"] .* 1e-3
        elseif var == "rlr"
            _data = (_data_pysdm["qc"] .+ _data_pysdm["qr"]) .* 1e-3
        elseif var == "rtr"
            _data = _rv .+ (_data_pysdm["qc"] .+ _data_pysdm["qr"]) .* 1e-3
        elseif var == "Nl"
            _data = _data_pysdm["nc"]
        elseif var == "Nr"
            _data = _data_pysdm["nr"]
        elseif var == "Na"
            _data = _data_pysdm["na"]
        elseif var == "rainrate"
            _data =
                _data_pysdm["qr"] .* 1e-3 .* _data_pysdm["rhod"] .* _data_pysdm["rain averaged terminal velocity"] .*
                FT(3600)
        elseif var == "Z"
            _data = _data_pysdm["radar_refl"]
        elseif var == "reff"
            _data = _data_pysdm["effective_radius"]
        elseif var == "rainrate_surface"
            _rainrate = 
                _data_pysdm["qr"] .* 1e-3 .* _data_pysdm["rhod"] .* _data_pysdm["rain averaged terminal velocity"] .*
                FT(3600)
            rainrate = 1.5 .* _rainrate[:, 1] - 0.5 .* _rainrate[:, 2]
            rainrate[findall(x -> x < 0, rainrate)] .= FT(0)
            _data = rainrate
        elseif var == "reff_top"
            _data = _data_pysdm["reff"]
        elseif var == "Z_top"
            _data = _data_pysdm["Z_top"]
        elseif var == "Z_mid"
            _data = _data_pysdm["Z_mid"]
        elseif var == "Z_bottom"
            _data = _data_pysdm["Z_bottom"]
        elseif var == "Z_1"
            _data = _data_pysdm["Z_1"]
        elseif var == "Z_2"
            _data = _data_pysdm["Z_2"]
        else
            _data = _data_pysdm[var]
        end

        _heights = heights[i]
        if length(size(_data)) == 1
            _f_1d = linear_interpolation(_t_data, _data, extrapolation_bc = Line())
            _f(t, z) = _f_1d(t)
        else
            _f = linear_interpolation((_t_data, _z_data), _data, extrapolation_bc = Line())
        end
        if apply_filter
            _n_z_data = length(_z_data)
            _dz_data = _n_z_data == 1 ? FT(0.5) : _z_data[2] - _z_data[1]
            _z_data_f = collect(range(_z_data[1] - _dz_data / 2, _z_data[end] + _dz_data / 2, _n_z_data + 1))

            _n_heights = length(_heights)
            _dz_heights = _n_heights == 1 ? FT(0.5) : _heights[2] - _heights[1]
            _heights_f = collect(range(_heights[1] - _dz_heights / 2, _heights[end] + _dz_heights / 2, _n_heights + 1))

            _t_data_c = [0.5 * (_t_data[i] + _t_data[i + 1]) for i in 1:(length(_t_data) - 1)]
            _data_c = [_f(t, z) for z in _z_data, t in _t_data_c]
            _output[var] = filter_field(_data_c, _t_data, times, _z_data_f, _heights_f)
        else
            _output[var] = [_f(t, z) for z in _heights, t in times]
        end
    end

    return _output
end

function filter_field(data::Array{FT}, t_data::Array, times::Array{FT}, z_data::Array, heights::Array{FT})
    @assert issorted(heights)
    @assert length(Set(heights)) == length(heights)
    @assert issorted(z_data)
    @assert length(Set(z_data)) == length(z_data)
    @assert z_data[1] <= heights[1]
    @assert heights[end] <= z_data[end]
    @assert issorted(times)
    @assert length(Set(times)) == length(times)
    @assert issorted(t_data)
    @assert length(Set(t_data)) == length(t_data)
    @assert t_data[1] <= times[1]
    @assert times[end] <= t_data[end]
    @assert size(data) == (length(z_data) - 1, length(t_data) - 1)

    _nz_calib = length(heights) - 1
    _nt_calib = length(times) - 1
    _output = zeros(_nz_calib, _nt_calib)
    for i in 1:_nz_calib
        ind_ini_z = findall(x -> x > heights[i], z_data)[1] - 1
        ind_end_z = findall(x -> x >= heights[i + 1], z_data)[1] - 1
        for j in 1:_nt_calib
            ind_ini_t = findall(x -> x > times[j], t_data)[1] - 1
            ind_end_t = findall(x -> x >= times[j + 1], t_data)[1] - 1

            _s_ij = FT(0)
            _s_times_output_ij = FT(0)
            for k in ind_ini_z:ind_end_z
                _z1 = (k == ind_ini_z) ? heights[i] : z_data[k]
                _z2 = (k == ind_end_z) ? heights[i + 1] : z_data[k + 1]
                for l in ind_ini_t:ind_end_t
                    _t1 = (l == ind_ini_t) ? times[j] : t_data[l]
                    _t2 = (l == ind_end_t) ? times[j + 1] : t_data[l + 1]
                    _s = (_z2 - _z1) * (_t2 - _t1)
                    _s_ij = _s_ij + _s
                    _s_times_output_ij = _s_times_output_ij + _s * data[k, l]
                end
            end
            _output[i, j] = _s_times_output_ij / _s_ij

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
    Random.seed!(config["observations"]["random_seed"])
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
