@testset "Make space" begin

    @test K2D.make_function_space(FT, xlim = (0, 100.0), zlim = (0, 200.0), helem = 16, velem = 32) isa
          Tuple{CC.Spaces.ExtrudedFiniteDifferenceSpace, CC.Spaces.FaceExtrudedFiniteDifferenceSpace}
end

@testset "Make rhs function" begin

    rhs = K2D.make_rhs_function(equil_moist_ρθq, precip_1m)
    @test typeof(rhs) <: Function
end

@testset "Initialise aux" begin

    ip = (;
        ρ = [1.0, 1.0],
        ρ_dry = [0.999, 1.0],
        θ_liq_ice = [440.0, 450.0],
        q_tot = [1e-3, 0.0],
        q_liq = [0.0, 0.0],
        q_ice = [0.0, 0.0],
        p = [101300.0, 90000.0],
        T = [300.0, 290.0],
        θ_dry = [440.0, 450],
        q_rai = [0.0, 0.0],
        q_sno = [0.0, 0.0],
        N_liq = [0.0, 0.0],
        N_rai = [0.0, 0.0],
        N_aer_0 = [0.0, 0.0],
        N_aer = [0.0, 0.0],
        S_ql_moisture = [0.0, 0.0],
        S_qi_moisture = [0.0, 0.0],
        S_qt_precip = [0.0, 0.0],
        S_ql_precip = [0.0, 0.0],
        S_qi_precip = [0.0, 0.0],
        S_qr_precip = [0.0, 0.0],
        S_qs_precip = [0.0, 0.0],
        S_Nl_precip = [0.0, 0.0],
        S_Nr_precip = [0.0, 0.0],
        S_Na_activation = [0.0, 0.0],
        S_Nl_activation = [0.0, 0.0],
        term_vel_rai = [0.0, 0.0],
        term_vel_sno = [0.0, 0.0],
        term_vel_N_rai = [0.0, 0.0],
    )
    space, face_space = K2D.make_function_space(FT)

    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, no_precip)
    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, K1D.EquilibriumMoisture())
    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, K1D.NonEquilibriumMoisture())

    aux = K2D.initialise_aux(FT, ip, params..., 100.0, 200.0, 0.0, 0.0, space, face_space, equil_moist_ρθq)
    @test aux isa NamedTuple
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.aerosol_variables isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test aux.precip_velocities isa CC.Fields.FieldVector
    @test aux.prescribed_velocity isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_variables) == 0
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

end

@testset "advection tendency" begin

    space, face_space = K2D.make_function_space(FT)
    coords = CC.Fields.coordinate_field(space)
    face_coords = CC.Fields.coordinate_field(face_space)
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params)
    init = map(coord -> CO.initial_condition_1d(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z), coords)
    t = 1.1 * kid_params.t1

    @test_throws Exception K2D.advection_tendency!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K2D.advection_tendency!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    # eq
    aux = K2D.initialise_aux(FT, init, params..., 3000.0, 3000.0, 0.0, 0.0, space, face_space, equil_moist_ρθq)
    K2D.precompute_aux_prescribed_velocity!(aux, t)

    Y = CO.initialise_state(equil_moist_ρθq, no_precip, init)
    dY = Y / 10
    K2D.advection_tendency!(equil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = CO.initialise_state(equil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K2D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    # Non-eq
    aux = K2D.initialise_aux(FT, init, params..., 3000.0, 3000.0, 0.0, 0.0, space, face_space, nequil_moist_ρθq)
    K2D.precompute_aux_prescribed_velocity!(aux, t)

    Y = CO.initialise_state(nequil_moist_ρθq, no_precip, init)
    dY = Y / 10
    K2D.advection_tendency!(nequil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = CO.initialise_state(nequil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K2D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

end
