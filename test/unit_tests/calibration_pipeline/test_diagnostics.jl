@testset "lwp rwp and rainrate" begin
    #setup
    config = get_config()
    config["observations"]["data_names"] = ["rho", "ql", "qr", "rainrate"]
    (n_c, n_z, n_t) = KCP.get_numbers_from_config(config)
    n_single_case = sum(n_z) * n_t
    vec = rand(n_c * n_single_case)
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z[1]
    _lwp = zeros(n_t, n_c)
    _rwp = zeros(n_t, n_c)
    _rainrate = zeros(n_t, n_c)

    for i in 1:n_c
        v_ = KCP.get_case_i_vec(vec, i, n_single_case)
        fields_ = KCP.get_single_case_fields(v_, n_z, n_t)
        ρ = fields_[1]
        ql = fields_[2]
        qr = fields_[3]
        _lwp[:, i] = sum(ql .* ρ, dims = 1) .* dz
        _rwp[:, i] = sum(qr .* ρ, dims = 1) .* dz
        _rr = fields_[4]
        rainrate_0 = 1.5 * _rr[1, :] - 0.5 * _rr[2, :]
        rainrate_0[findall(x -> x < 0, rainrate_0)] .= Float64(0)
        _rainrate[:, i] = rainrate_0
    end

    #action
    lwp = KCP.lwp(vec, config)
    rwp = KCP.rwp(vec, config)
    rainrate = KCP.rainrate(vec, config, height = 0.0)

    #test
    @test lwp == _lwp
    @test rwp == _rwp
    @test rainrate == _rainrate

    #action
    config["observations"]["data_names"] = ["rho", "ql", "qr"]

    #test
    @test_throws Exception KCP.rainrate(vec, config, height = 0.0)

    #action
    config["observations"]["data_names"] = ["ql", "qr", "rainrate"]

    #test
    @test_throws Exception KCP.lwp(vec, config)
    @test_throws Exception KCP.rwp(vec, config)
end
