"""
Test that initial profiles match between CliMA and PySDM
"""

include("data_utils.jl")
const FT = Float64

function compare_profiles(; is_dry_flag::Bool)
    # Computational domain ...
    z_min = FT(0)
    z_max = FT(2220)
    n_elem = 222
    # ... and the created coordinates
    space, face_space = KiD.make_function_space(FT, z_min, z_max, n_elem)
    coord = CC.Fields.coordinate_field(space)
    face_coord = CC.Fields.coordinate_field(face_space)

    # Solve the initial value problem for density profile
    ρ_profile = KiD.ρ_ivp(FT, params, dry = is_dry_flag)
    # Create the initial condition profiles
    init = map(coord -> KiD.init_1d_column(FT, params, ρ_profile, coord.z, dry = is_dry_flag), coord)

    # Store the CliMA KiD initial profiles
    z_centers = parent(CC.Fields.coordinate_field(space))
    T = parent(init.T)
    p = parent(init.p)
    ρ = parent(init.ρ)
    θ_dry = parent(init.θ_dry)
    q_vap = parent(init.q_tot) - parent(init.q_liq) - parent(init.q_ice)
    q_liq = parent(init.q_liq)
    KM_data = (; z_centers, q_vap, ρ, θ_dry, T, p, q_liq)

    # Read in the PySDM KiD initial profiles
    sdm_case = is_dry_flag ? "dry" : "wet"
    sdm_data = load_sdm_data(sdm_case)

    # Interpolate both (we are not doing computations on the same grid)
    KM_T = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.T[:], extrapolation_bc = IP.Line())
    KM_p = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.p[:], extrapolation_bc = IP.Line())
    KM_ρ = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.ρ[:], extrapolation_bc = IP.Line())
    KM_q_liq = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.q_liq[:], extrapolation_bc = IP.Line())
    KM_q_vap = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.q_vap[:], extrapolation_bc = IP.Line())
    KM_θ_dry = IP.LinearInterpolation(KM_data.z_centers[:], KM_data.θ_dry[:], extrapolation_bc = IP.Line())

    SD_T = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.T_sdm, extrapolation_bc = IP.Line())
    SD_p = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.P_sdm, extrapolation_bc = IP.Line())
    SD_ρ = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.rho_sdm, extrapolation_bc = IP.Line())
    SD_q_liq = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.ql_sdm, extrapolation_bc = IP.Line())
    SD_q_vap = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.qv_sdm, extrapolation_bc = IP.Line())
    SD_θ_dry = IP.LinearInterpolation(sdm_data.z_sdm, sdm_data.thetad_sdm, extrapolation_bc = IP.Line())

    # Test that the initial profiles match between CliMA and PySDM
    test_title = "Case: " * sdm_case * ". Initial profiles of water (vapour and liquid):"
    @testset "$test_title" begin
        z_test = z_min:100:z_max
        @test all(isapprox(KM_q_vap(z_test), SD_q_vap(z_test), rtol = 1e-6))
        @test all(isapprox(KM_q_liq(z_test), SD_q_liq(z_test), rtol = 1e-6))
    end

    test_title = "Case: " * sdm_case * ". Initial profiles of (T, p, ρ):"
    if is_dry_flag
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test all(isapprox(KM_T(z_test), SD_T(z_test), rtol = 1e-3))
            @test all(isapprox(KM_p(z_test), SD_p(z_test), rtol = 1e-3))
            @test all(isapprox(KM_ρ(z_test), SD_ρ(z_test), rtol = 1e-3))
        end
    else
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test_broken all(isapprox(KM_T(z_test), SD_T(z_test), rtol = 1e-3))
            @test_broken all(isapprox(KM_p(z_test), SD_p(z_test), rtol = 1e-3))
            @test_broken all(isapprox(KM_ρ(z_test), SD_ρ(z_test), rtol = 1e-3))
        end
    end

    test_title = "Case: " * sdm_case * ". Initial profiles of θ_dry:"
    if is_dry_flag
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test all(isapprox(KM_θ_dry(z_test), SD_θ_dry(z_test), rtol = 1e-4))
        end
    else
        @testset "$test_title" begin
            z_test = z_min:100:z_max
            @test_broken all(isapprox(KM_θ_dry(z_test), SD_θ_dry(z_test), rtol = 1e-4))
        end
    end

    # Plot the profiles - TODO connect with buildkite artifacts
    plot_comparison(KM_data, sdm_case = sdm_case, dir = sdm_case)
end

compare_profiles(is_dry_flag = true)
compare_profiles(is_dry_flag = false)
