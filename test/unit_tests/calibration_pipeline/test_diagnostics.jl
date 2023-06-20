@testset "lwp rwp and rainrate" begin
    #setup
    config = get_config()
    config["observations"]["data_names"] = ["rho", "ql", "qr", "rain averaged terminal velocity"]
    (n_c, n_v, n_z, n_t) = KID.get_numbers_from_config(config)
    n_vht = n_v * n_z * n_t
    vec = rand(n_c * n_vht)
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z
    _lwp = zeros(n_t, n_c)
    _rwp = zeros(n_t, n_c)
    _rainrate = zeros(n_t, n_c)

    for i in 1:n_c
        v_ = vec[((i - 1) * n_vht + 1):(i * n_vht)]
        m_ = reshape(v_, n_v * n_z, n_t)
        ρ = m_[1:n_z, :]
        ql = m_[(n_z + 1):(2 * n_z), :]
        qr = m_[(2 * n_z + 1):(3 * n_z), :]
        _lwp[:, i] = sum(ql .* ρ, dims = 1) .* dz
        _rwp[:, i] = sum(qr .* ρ, dims = 1) .* dz
        ρ0 = 1.5 * ρ[1, :] - 0.5 * ρ[2, :]
        qr0 = 1.5 * qr[1, :] - 0.5 * qr[2, :]
        vt0 = 1.5 * m_[3 * n_z + 1, :] - 0.5 * m_[3 * n_z + 2, :]
        ρ0[findall(x -> x < 0, ρ0)] .= Float64(0)
        qr0[findall(x -> x < 0, qr0)] .= Float64(0)
        vt0[findall(x -> x < 0, vt0)] .= Float64(0)
        _rainrate[:, i] = qr0 .* ρ0 .* vt0 .* 3600
    end

    #action
    lwp = KID.lwp(vec, config)
    rwp = KID.rwp(vec, config)
    rainrate = KID.rainrate(vec, config, height = 0.0, isref = true)

    #test
    @test lwp == _lwp
    @test rwp == _rwp
    @test rainrate == _rainrate

    #setup
    mat = reshape(vec, n_z, n_v, n_t, n_c)
    mat[:, 3, :, :] .= 0.0
    vec = reshape(mat, n_z * n_v * n_t * n_c, 1)[:]

    #action
    lwp = KID.lwp(vec, config)
    rwp = KID.rwp(vec, config)
    rainrate = KID.rainrate(vec, config, height = 0.0, isref = false)

    #test
    @test rwp == zeros(n_t, n_c)
    @test rainrate == zeros(n_t, n_c)
end
