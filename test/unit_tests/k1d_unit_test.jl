params = (common_params, kid_params, thermo_params, air_params, activation_params)

@testset "Parameter overwrites" begin

    @test TD.Parameters.R_d(thermo_params) == 8.314462618 / 0.02896998
    @test TD.Parameters.R_v(thermo_params) == 8.314462618 / 0.018015

    @test TD.Parameters.MSLP(thermo_params) == 100000.0

end

@testset "Make space" begin

    @test K1D.make_function_space(FT, 0, 100, 10) isa
          Tuple{CC.Spaces.CenterFiniteDifferenceSpace, CC.Spaces.FaceFiniteDifferenceSpace}
end

@testset "Make rhs function" begin

    rhs = K1D.make_rhs_function(equil_moist, precip_1m)
    @test typeof(rhs) <: Function
end

@testset "Initialise aux" begin

    _ip = (;
        ρ = 1.2,
        ρ_dry = 1.185,
        θ_liq_ice = 350.0,
        q_tot = 15e-3 / 1.2,
        q_liq = 1e-3 / 1.2,
        q_ice = 2e-3 / 1.2,
        q_rai = 1e-3 / 1.2,
        q_sno = 2e-3 / 1.2,
        q_rim = 0.5e-3 / 1.2,
        q_liqonice = 0.5e-3 / 1.2,
        q_vap = 1e-2 / 1.2,
        ρq_tot = 15e-3,
        ρq_liq = 1e-3,
        ρq_ice = 2e-3,
        ρq_rai = 1e-3,
        ρq_sno = 2e-3,
        ρq_rim = 0.5e-3,
        ρq_liqonice = 0.5e-3,
        ρq_vap = 1e-2,
        p = 101300.0,
        T = 280.0,
        θ_dry = 360.0,
        N_liq = 1e8,
        N_ice = 1e8,
        N_rai = 1e4,
        N_sno = 1e4,
        N_aer = 1e10,
        B_rim = 0.5e-3 / 900,
        zero = 0.0,
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(0),
        CC.Geometry.ZPoint{FT}(1),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = CC.Spaces.FaceFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)
    ip = map(coord -> _ip, coord)

    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, no_precip)
    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, CO.EquilibriumMoisture())
    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, CO.NonEquilibriumMoisture())

    ms_styles = [equil_moist, nequil_moist]
    ps_styles = [no_precip, precip_1m, precip_2m, precip_2m]
    for (ms, ps) in zip(ms_styles, ps_styles)
        aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, ms, ps)
        @test aux isa NamedTuple
        @test LA.norm(aux.precip_sources) == 0
        if ms == nequil_moist
            @test LA.norm(aux.cloud_sources) == 0
        else
            @test aux.cloud_sources === nothing
        end
    end
    ms = p3_moist
    ps = precip_p3
    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, ms, ps)
    @test aux isa NamedTuple
    @test LA.norm(aux.precip_sources) == 0
end

@testset "advection_tendency" begin

    space, face_space = K1D.make_function_space(FT, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    init =
        map(coord -> CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, equil_moist, no_precip)
    Y = CO.initialise_state(equil_moist, no_precip, init)
    dY = Y / 10
    @test_throws Exception K1D.advection_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K1D.advection_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, t)

    ms_styles = [equil_moist, nequil_moist]
    ps_styles = [no_precip, precip_1m, precip_2m, precip_2m]
    for (ms, ps) in zip(ms_styles, ps_styles)
        aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, ms, ps)
        Y = CO.initialise_state(ms, ps, init)
        dY = Y / 10
        ρw = 0.0
        @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
        aux.prescribed_velocity.ρw0 = ρw
        K1D.advection_tendency!(ms, dY, Y, aux, t)
        K1D.advection_tendency!(ps, dY, Y, aux, t)
        @test dY ≈ Y / 10 atol = eps(FT) * 10
    end
end

@testset "advection_tendency p3" begin

    space, face_space = K1D.make_function_space(FT, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ice_start = false
    _magnitude = Float64(0.5)
    _q_flux = Float64(0.65e-4)
    _N_flux = Float64(40000)
    _F_rim = Float64(0.2)
    _F_liq = Float64(0.2)
    _ρ_r_init = Float64(900)
    init = map(
        coord -> CO.p3_initial_condition(
            FT,
            kid_params,
            thermo_params,
            coord.z;
            _q_init = _q_flux,
            _N_init = _N_flux,
            _F_rim = _F_rim,
            _F_liq = _F_liq,
            _ρ_r = _ρ_r_init,
            z_top = FT(3000),
            ice_start = ice_start,
        ),
        coord,
    )
    t = 13.0

    # eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, p3_moist, precip_p3)
    Y = CO.initialise_state(p3_moist, precip_p3, init)
    dY = Y / 10
    @test_throws Exception K1D.advection_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K1D.advection_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, t)

    ms = p3_moist
    ps = precip_p3
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, ms, ps)
    Y = CO.initialise_state(ms, ps, init)
    dY = Y / 10
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw
    K1D.advection_tendency!(ms, dY, Y, aux, t)
    K1D.advection_tendency!(ps, dY, Y, aux, t)
    @test dY >= Y / 10
end

@testset "aerosol activation for 2M schemes" begin
    #setup
    _ip = (;
        ρ = 1.2,
        ρ_dry = 1.185,
        θ_liq_ice = 350.0,
        q_tot = 15e-3 / 1.2,
        q_liq = 1e-3 / 1.2,
        q_ice = 2e-3 / 1.2,
        q_rai = 1e-3 / 1.2,
        q_sno = 2e-3 / 1.2,
        ρq_tot = 15e-3,
        ρq_liq = 1e-3,
        ρq_ice = 2e-3,
        ρq_rai = 1e-3,
        ρq_sno = 2e-3,
        p = 101300.0,
        T = 280.0,
        θ_dry = 360.0,
        N_liq = 1e8,
        N_ice = 1e8,
        N_rai = 1e4,
        N_sno = 1e4,
        N_aer = 5e7,
        zero = 0.0,
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(0),
        CC.Geometry.ZPoint{FT}(1),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = CC.Spaces.FaceFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)
    ip = map(coord -> _ip, coord)

    get_value(x) = parent(CC.Fields.field_values(CC.Fields.level(x, 1)))

    TS = CO.TimeStepping(FT(10), FT(10), FT(20))

    Y = CO.initialise_state(equil_moist, precip_2m, ip)
    dY = similar(Y)
    aux = K1D.initialise_aux(
        FT,
        ip,
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        TS,
        0.0,
        face_space,
        equil_moist,
        precip_2m,
    )

    K1D.precompute_aux_activation!(precip_2m, dY, Y, aux, 13)

    #test
    @test all(isfinite, get_value(aux.activation_sources.N_aer))
    @test all(isfinite, get_value(aux.activation_sources.N_liq))
    @test get_value(aux.activation_sources.N_aer) == -get_value(aux.activation_sources.N_liq)
end
