function get_limits(u::Vector{Matrix{Float64}})
    _n_matrices = length(u)
    _n_rows = size(u[1])[1]

    _minima_vec = [minimum(u_i, dims = 2) for u_i in u]
    _maxima_vec = [maximum(u_i, dims = 2) for u_i in u]
    _minima_mat = Matrix{Float64}(undef, _n_rows, _n_matrices)
    _maxima_mat = Matrix{Float64}(undef, _n_rows, _n_matrices)
    for i in 1:_n_matrices
        _minima_mat[:, i] = _minima_vec[i]
        _maxima_mat[:, i] = _maxima_vec[i]
    end
    _minima = minimum(_minima_mat, dims = 2)
    _maxima = maximum(_maxima_mat, dims = 2)
    _ext = [(_minima[i], _maxima[i]) for i in 1:_n_rows]
    return _ext
end

function make_filter_props(
    n_z,
    t_calib;
    apply = false,
    nz_per_filtered_cell = ones(Int, length(n_z)),
    nt_per_filtered_cell = 1,
)

    @assert all(n_z .% nz_per_filtered_cell .== 0)

    filter = Dict()
    filter["apply"] = apply
    filter["nz_per_filtered_cell"] = nz_per_filtered_cell
    filter["nt_per_filtered_cell"] = nt_per_filtered_cell
    filter["nz_unfiltered"] = n_z
    filter["nz_filtered"] = Int.(n_z ./ nz_per_filtered_cell)
    filter["nt_filtered"] = length(t_calib) - 1

    saveat = Float64[]
    for i in 1:(length(t_calib) - 1)
        dt = (t_calib[i + 1] - t_calib[i]) / filter["nt_per_filtered_cell"]
        saveat = [saveat; collect(range(t_calib[i] + dt / 2, t_calib[i + 1] - dt / 2, filter["nt_per_filtered_cell"]))]
    end
    filter["saveat_t"] = saveat

    return filter
end

function get_numbers_from_config(config::Dict)

    n_cases = length(config["observations"]["cases"])
    n_heights =
        config["model"]["filter"]["apply"] ? config["model"]["filter"]["nz_filtered"] :
        config["model"]["filter"]["nz_unfiltered"]
    _n_times_default =
        config["model"]["filter"]["apply"] ? length(config["model"]["t_calib"]) - 1 : length(config["model"]["t_calib"])
    n_times = Int[]
    for case in config["observations"]["cases"]
        if :t_cal in collect(keys(case))
            _n_times_case_i = config["model"]["filter"]["apply"] ? length(case.t_cal) - 1 : length(case.t_cal)
        else
            _n_times_case_i = _n_times_default
        end
        n_times = [n_times; _n_times_case_i]
    end

    @assert length(n_heights) == length(config["observations"]["data_names"])
    return (; n_cases, n_heights, n_times)
end

function get_case_i_vec(vec::Vector{FT}, i::Int, n_single_case::Vector{Int}) where {FT <: Real}
    return vec[(sum(n_single_case[1:(i - 1)]) + 1):sum(n_single_case[1:i])]
end

function get_single_case_fields(vec_single_case::Vector{FT}, n_heights::Vector{Int}, n_times::Int) where {FT <: Real}

    n_variables = length(n_heights)
    n_single_time = sum(n_heights)

    if length(vec_single_case) != n_times * n_single_time
        error("Inconsistent vector size and provided numbers!")
    end

    m_ = reshape(vec_single_case, n_single_time, n_times)
    fields = ntuple(n_variables) do j
        _last_ind = j == 1 ? 0 : sum(n_heights[1:(j - 1)])
        m_[(_last_ind + 1):(_last_ind + n_heights[j]), :]
    end
    return fields
end

function make_block_diagonal_matrix(a::AbstractMatrix, b::AbstractMatrix)

    if (a isa Diagonal) && (b isa Diagonal)
        return Diagonal([diag(a); diag(b)])
    end

    (nr_a, nc_a) = size(a)
    (nr_b, nc_b) = size(b)

    nr_c = nr_a + nr_b
    nc_c = nc_a + nc_b
    c = zeros(nr_c, nc_c)
    c[1:nr_a, 1:nc_a] = a
    c[(nr_a + 1):nr_c, (nc_a + 1):nc_c] = b
    return c
end

"""
    compute_error_metrics(ϕ_values, ϕ_names, config, RS)

Returns loss ||G - y||_Γ, error with std normalization ||G - y||, and error based on max normalization ||G - y||.

# Inputs:
- `ϕ_values` :: vector containing parameter values
- `ϕ_names` :: vector containing parameter names
- `config` :: dictionary containing settings
- `RS` :: ReferenceStatistics object containing statistics of the observations
"""
function compute_error_metrics(
    ϕ_values::Array{FT, 1},
    ϕ_names::Array{String, 1},
    config::Dict,
    RS::ReferenceStatistics,
) where {FT <: Real}

    n_cases = length(RS.case_numbers)

    G = run_dyn_model(ϕ_values, ϕ_names, config; case_numbers = RS.case_numbers)
    yg_diff = RS.centered ? G : G - RS.obs_mean
    yg_diff_mean_normalized = normalize_sim(yg_diff, RS, norm = RS.var_mean_max)
    yg_diff_std_normalized = normalize_sim(yg_diff, RS, norm = RS.var_std_max)
    mse_m = dot(yg_diff_mean_normalized, yg_diff_mean_normalized) / n_cases
    mse_s = dot(yg_diff_std_normalized, yg_diff_std_normalized) / n_cases

    G_n = pca_transform(normalize_sim(G, RS), RS)
    yg_diff_n = G_n - RS.y
    loss = dot(yg_diff_n, RS.Γ \ yg_diff_n) / n_cases

    return (loss = loss, mse_m = mse_m, mse_s = mse_s)
end

"""
    compute_error_metrics(ϕ_values, ϕ_names, config, obs)

Returns loss ||G - y||_Γ, error with std normalization ||G - y||, and error based on max normalization ||G - y||.
This function computes error metrics by computing errors over all cases seperately and then averaging.

# Inputs:
- `ϕ_values` :: vector containing parameter values
- `ϕ_names` :: vector containing parameter names
- `config` :: dictionary containing settings
- `obs` :: matrix of observations
"""
function compute_error_metrics(
    ϕ_values::Array{FT, 1},
    ϕ_names::Array{String, 1},
    config::Dict,
    obs::Matrix{FT},
) where {FT <: Real}
    n_cases = length(config["observations"]["cases"])

    ref_stats_list = make_ref_stats_list(obs, config["statistics"], get_numbers_from_config(config)...)
    model_error = zeros(3)
    for i in 1:n_cases
        model_error_ = compute_error_metrics(ϕ_values, ϕ_names, config, ref_stats_list[i])
        model_error[1] = model_error[1] + model_error_.loss
        model_error[2] = model_error[2] + model_error_.mse_m
        model_error[3] = model_error[3] + model_error_.mse_s
    end
    model_error[1] = model_error[1] / n_cases
    model_error[2] = model_error[2] / n_cases
    model_error[3] = model_error[3] / n_cases
    return (loss = model_error[1], mse_m = model_error[2], mse_s = model_error[3])
end
