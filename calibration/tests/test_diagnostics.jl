TT.@testset "liquid water path, rainwater path and rainrate" begin
    #setup
    config = get_config()
    config["observations"]["data_names"] = ["rho", "ql", "qr", "rainrate"]
    config["model"]["filter"] =
        make_filter_props(config["model"]["n_elem"] .* ones(Int, 4), config["model"]["t_calib"])
    (n_c, n_z, n_t) = get_numbers_from_config(config)
    n_single_case = sum(n_z) * n_t
    vec = rand(sum(n_single_case))
    dz = (config["model"]["z_max"] - config["model"]["z_min"]) / n_z[1]
    _lwp = []
    _rwp = []
    _rr = []

    for i in 1:n_c
        v_ = get_case_i_vec(vec, i, n_single_case)
        fields_ = get_single_case_fields(v_, n_z, n_t[i])
        ρ = fields_[1]
        ql = fields_[2]
        qr = fields_[3]
        push!(_lwp, sum(ql .* ρ, dims = 1) .* dz)
        push!(_rwp, sum(qr .* ρ, dims = 1) .* dz)
        rr_0 = 1.5 * fields_[4][1, :] - 0.5 * fields_[4][2, :]
        rr_0[findall(x -> x < 0, rr_0)] .= Float64(0)
        push!(_rr, rr_0)
    end

    #action
    lwp = liquid_water_path(vec, config)
    rwp = rainwater_path(vec, config)
    rr = rainrate(vec, config, height = 0.0)

    #test
    TT.@test lwp == _lwp
    TT.@test rwp == _rwp
    TT.@test rr == _rr

    #action
    config["observations"]["data_names"] = ["rho", "ql", "qr"]

    #test
    TT.@test_throws Exception rainrate(vec, config, height = 0.0)

    #action
    config["observations"]["data_names"] = ["ql", "qr", "rainrate"]

    #test
    TT.@test_throws Exception liquid_water_path(vec, config)
    TT.@test_throws Exception rainwater_path(vec, config)
end
