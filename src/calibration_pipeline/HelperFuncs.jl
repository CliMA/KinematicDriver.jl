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
