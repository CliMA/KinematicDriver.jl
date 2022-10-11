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

function make_filter_props(n_elem, t_calib; apply = false, nz_per_filtered_cell = 1, nt_per_filtered_cell = 1)

    @assert n_elem % nz_per_filtered_cell == 0

    filter = Dict()
    filter["apply"] = apply
    filter["nz_per_filtered_cell"] = nz_per_filtered_cell
    filter["nt_per_filtered_cell"] = nt_per_filtered_cell
    filter["nz_filtered"] = Int(n_elem / nz_per_filtered_cell)
    filter["nt_filtered"] = length(t_calib) - 1

    saveat = Float64[]
    for i in 1:(length(t_calib) - 1)
        dt = (t_calib[i + 1] - t_calib[i]) / filter["nt_per_filtered_cell"]
        saveat = [saveat; collect(range(t_calib[i] + dt / 2, t_calib[i + 1] - dt / 2, filter["nt_per_filtered_cell"]))]
    end
    filter["saveat_t"] = saveat

    return filter
end
