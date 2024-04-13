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

    rhs = K1D.make_rhs_function(equil_moist_ρθq, precip_1m)
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
    space, face_space = K1D.make_function_space(FT, 0, 100, 10)

    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, no_precip)
    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, CO.EquilibriumMoisture())
    @test_throws Exception K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, CO.NonEquilibriumMoisture())

    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, equil_moist_ρθq)
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

    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, equil_moist_ρdTq)
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

    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, nequil_moist_ρθq)
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

    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, nequil_moist_ρdTq)
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

@testset "advection_tendency" begin

    space, face_space = K1D.make_function_space(FT, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = K1D.ρ_ivp(FT, kid_params, thermo_params)
    init = map(coord -> K1D.init_1d_column(FT, common_params, kid_params, thermo_params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, equil_moist_ρθq)
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = CO.initialise_state(equil_moist_ρθq, no_precip, init)
    dY = Y / 10

    @test_throws Exception K1D.advection_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K1D.advection_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, t)

    K1D.advection_tendency!(equil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = CO.initialise_state(equil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K1D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    # Non-eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, nequil_moist_ρθq)
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = CO.initialise_state(nequil_moist_ρθq, no_precip, init)
    dY = Y / 10

    K1D.advection_tendency!(nequil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = CO.initialise_state(nequil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K1D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

end

@testset "aerosol activation for 2M schemes" begin
    #setup
    q_tot = 1e-2
    q_liq = 1e-3
    N_aer = 5e7
    N_aer_0 = 1e8
    T = 280.0
    p = 1e5
    ρ = 1.0
    ρw = 2.0
    dt = 1.0

    #action
    tmp = K1D.aerosol_activation_helper(kid_params, thermo_params, air_params, activation_params, q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, ρw, dt)

    #test
    @test tmp isa NamedTuple
    @test tmp.S_Nl ≈ -tmp.S_Na
    @test 0 < tmp.S_Nl < N_aer_0

    #action
    tmp = K1D.aerosol_activation_helper(kid_params, thermo_params, air_params, activation_params, q_tot, q_liq, N_aer, N_aer_0, 290.0, p, ρ, ρw, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(FT)
    @test tmp.S_Na ≈ 0.0 atol = eps(FT)

    #action
    tmp = K1D.aerosol_activation_helper(kid_params, thermo_params, air_params, activation_params, q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, 0.0, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(FT)
    @test tmp.S_Na ≈ 0.0 atol = eps(FT)

end
