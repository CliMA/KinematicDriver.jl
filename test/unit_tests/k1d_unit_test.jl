const kid_dir = pkgdir(Kinematic1D)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

# override the defaults
default_toml_dict = CP.create_toml_dict(FT, dict_type = "alias")
toml_dict = override_toml_dict(@__DIR__, default_toml_dict)

# create all the parameters structs ...
thermo_params = create_thermodynamics_parameters(toml_dict)
kid_params = create_kid_parameters(toml_dict)
air_params = CMP.AirProperties(FT, toml_dict)
activation_params = CMP.AerosolActivationParameters(FT, toml_dict)
params = (kid_params, thermo_params, air_params, activation_params)
# ... for cloud condensate options ...
equil_moist_ρθq = CO.EquilibriumMoisture_ρθq()
equil_moist_ρdTq = CO.EquilibriumMoisture_ρdTq()
nequil_moist_ρθq = CO.NonEquilibriumMoisture_ρθq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
nequil_moist_ρdTq = CO.NonEquilibriumMoisture_ρdTq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
# ... and precipitation options
no_precip = CO.NoPrecipitation()
precip_0m = CO.Precipitation0M(CMP.Parameters0M(FT, toml_dict))
precip_1m = CO.Precipitation1M(
    CMP.CloudLiquid(FT, toml_dict),
    CMP.CloudIce(FT, toml_dict),
    CMP.Rain(FT, toml_dict),
    CMP.Snow(FT, toml_dict),
    CMP.CollisionEff(FT, toml_dict),
    CMP.KK2000(FT, toml_dict),
    CMP.Blk1MVelType(FT, toml_dict),
)
precip_2m = CO.Precipitation2M(CMP.SB2006(FT, toml_dict), CMP.SB2006VelType(FT, toml_dict))

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

@testset "Initialise state" begin

    @test_throws Exception K1D.initialise_state(K1D.AbstractMoistureStyle(), K1D.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0, N_liq = 0, N_rai = 0, N_aer = 0)

    state = K1D.initialise_state(equil_moist_ρθq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = K1D.initialise_state(equil_moist_ρθq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = K1D.initialise_state(equil_moist_ρθq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = K1D.initialise_state(equil_moist_ρθq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = K1D.initialise_state(equil_moist_ρdTq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = K1D.initialise_state(equil_moist_ρdTq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = K1D.initialise_state(equil_moist_ρdTq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = K1D.initialise_state(equil_moist_ρdTq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = K1D.initialise_state(nequil_moist_ρθq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = K1D.initialise_state(nequil_moist_ρθq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = K1D.initialise_state(nequil_moist_ρθq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = K1D.initialise_state(nequil_moist_ρθq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = K1D.initialise_state(nequil_moist_ρdTq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = K1D.initialise_state(nequil_moist_ρdTq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = K1D.initialise_state(nequil_moist_ρdTq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = K1D.initialise_state(nequil_moist_ρdTq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

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
        S_Na = [0.0, 0.0],
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

@testset "Zero tendencies" begin

    ip = (;
        ρ = [1.0, 1.0],
        ρ_dry = [0.999, 1.0],
        θ_liq_ice = [440.0, 450.0],
        q_tot = [1e-3, 0.0],
        q_liq = [0.0, 0.0],
        q_ice = [0.0, 0.0],
        ρq_tot = [1e-3, 0.0],
        ρq_liq = [0.0, 0.0],
        ρq_ice = [0.0, 0.0],
        ρq_rai = [0.0, 0.0],
        ρq_sno = [0.0, 0.0],
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
        S_Na = [0.0, 0.0],
    )
    space, face_space = K1D.make_function_space(FT, 0, 100, 10)
    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, equil_moist_ρθq)
    Y = K1D.initialise_state(equil_moist_ρθq, no_precip, ip)

    dY = (; ρq_tot = [10.0, 13.0])
    @test_throws Exception K1D.zero_tendencies!(CO.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception K1D.zero_tendencies!(CO.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)

    dY = (; ρq_tot = [10.0, 13.0])
    K1D.zero_tendencies!(equil_moist_ρθq, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    aux = K1D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, nequil_moist_ρθq)
    Y = K1D.initialise_state(nequil_moist_ρθq, precip_1m, ip)
    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
    )
    K1D.zero_tendencies!(nequil_moist_ρθq, dY, Y, aux, 1.0)
    K1D.zero_tendencies!(precip_1m, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        N_liq = [1e8, 1e7],
        N_rai = [1e4, 1e4],
        N_aer = [1e9, 1e8],
    )
    K1D.zero_tendencies!(nequil_moist_ρθq, dY, Y, aux, 1.0)
    K1D.zero_tendencies!(precip_2m, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

end

@testset "Tendency helper functions" begin

    ρ = 1.2
    ρ_dry = 1.185
    θ_liq_ice = 350.0
    ρq_tot = 15e-3
    ρq_liq = 1e-3
    ρq_ice = 2e-3
    ρq_rai = 1e-3
    ρq_sno = 2e-3
    q_tot = ρq_tot / ρ
    q_liq = ρq_liq / ρ
    q_ice = ρq_ice / ρ
    q_rai = ρq_rai / ρ
    q_sno = ρq_sno / ρ
    N_liq = 1e8
    N_rai = 1e4
    T = 280.0

    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    tmp = K1D.moisture_helper_vars_eq_ρθq(thermo_params, ρq_tot, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = K1D.moisture_helper_vars_eq_ρdTq(thermo_params, ρq_tot, ρ_dry, T)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = K1D.moisture_helper_vars_neq_ρθq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T .+ tmp.θ_dry)
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = K1D.moisture_helper_vars_neq_ρdTq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ_dry, T)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice .+ tmp.θ_dry)
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = K1D.moisture_helper_sources(thermo_params, nequil_moist_ρθq, ρ, T, q_tot, q_liq, q_ice)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_ice)

    tmp = K1D.precip_helper_vars(ρq_rai, ρq_sno, ρ)
    @test !isnan(tmp.q_rai .+ tmp.q_sno)
    @test tmp.q_rai == q_rai
    @test tmp.q_sno == q_sno

    dt = 0.5
    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    tmp = K1D.precip_helper_sources!(precip_0m, thermo_params, ts, q_tot, q_liq, q_ice, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot == tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot == -(tmp.S_q_rai + tmp.S_q_sno)

    tmp = K1D.precip_helper_sources!(
        precip_1m,
        thermo_params,
        kid_params,
        air_params,
        ts,
        q_tot,
        q_liq,
        q_ice,
        q_rai,
        q_sno,
        T,
        ρ,
        dt,
    )
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_ice + tmp.S_q_vap
    @test tmp.S_q_tot ≈ -(tmp.S_q_rai + tmp.S_q_sno)

    tmp = K1D.precip_helper_sources!(
        precip_2m,
        thermo_params,
        kid_params,
        air_params,
        q_tot,
        q_liq,
        q_rai,
        N_liq,
        N_rai,
        T,
        ρ,
        dt,
    )
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_rai)
    @test !isnan(tmp.S_N_liq .+ tmp.S_N_rai)
    @test !isnan(tmp.term_vel_rai .+ tmp.term_vel_N_rai)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_vap
    @test tmp.S_q_tot ≈ -tmp.S_q_rai

end

@testset "advection_tendency and sources_tendency" begin

    space, face_space = K1D.make_function_space(FT, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = K1D.ρ_ivp(FT, kid_params, thermo_params)
    init = map(coord -> K1D.init_1d_column(FT, kid_params, thermo_params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, equil_moist_ρθq)
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = K1D.initialise_state(equil_moist_ρθq, no_precip, init)
    dY = Y / 10

    @test_throws Exception advection_tendency!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception advection_tendency!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception K1D.sources_tendency!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K1D.sources_tendency!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    K1D.advection_tendency!(equil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = K1D.initialise_state(equil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K1D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    # Non-eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, nequil_moist_ρθq)
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = K1D.initialise_state(nequil_moist_ρθq, no_precip, init)
    dY = Y / 10

    K1D.advection_tendency!(nequil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = K1D.initialise_state(nequil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K1D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

end

@testset "precompute aux" begin

    space, face_space = K1D.make_function_space(FT, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = K1D.ρ_ivp(FT, kid_params, thermo_params)
    init = map(coord -> K1D.init_1d_column(FT, kid_params, thermo_params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, equil_moist_ρθq)

    Y = K1D.initialise_state(equil_moist_ρθq, no_precip, init)
    dY = Y / 10

    @test_throws Exception K1D.precompute_aux_thermo!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K1D.precompute_aux_precip!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    K1D.precompute_aux_thermo!(equil_moist_ρθq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    K1D.precompute_aux_thermo!(equil_moist_ρdTq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    # Non-eq
    aux = K1D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, nequil_moist_ρθq)
    Y = K1D.initialise_state(nequil_moist_ρθq, no_precip, init)
    dY = Y / 10

    K1D.precompute_aux_thermo!(nequil_moist_ρθq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    K1D.precompute_aux_thermo!(nequil_moist_ρdTq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    aux = K1D.initialise_aux(FT, init, params..., (; dt = 2.0), 0.0, face_space, nequil_moist_ρθq)
    Y = K1D.initialise_state(nequil_moist_ρθq, precip_2m, init)
    dY = Y / 10

    K1D.precompute_aux_precip!(precip_2m, dY, Y, aux, t)
    tmp = @. K1D.precip_helper_sources!(
        precip_2m,
        aux.thermo_params,
        aux.kid_params,
        aux.air_params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.precip_variables.N_liq,
        aux.precip_variables.N_rai,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
    )

    @test aux.precip_sources.S_q_rai ≈ tmp.S_q_rai atol = eps(FT) * 10
    @test aux.precip_sources.S_q_tot ≈ tmp.S_q_tot atol = eps(FT) * 10
    @test aux.precip_sources.S_q_liq ≈ tmp.S_q_liq atol = eps(FT) * 10
    @test aux.precip_sources.S_N_liq ≈ tmp.S_N_liq atol = eps(FT) * 10
    @test aux.precip_sources.S_N_rai ≈ tmp.S_N_rai atol = eps(FT) * 10
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
    tmp = K1D.aerosol_activation_helper(params..., q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, ρw, dt)

    #test
    @test tmp isa NamedTuple
    @test tmp.S_Nl ≈ -tmp.S_Na
    @test 0 < tmp.S_Nl < N_aer_0

    #action
    tmp = K1D.aerosol_activation_helper(params..., q_tot, q_liq, N_aer, N_aer_0, 290.0, p, ρ, ρw, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(FT)
    @test tmp.S_Na ≈ 0.0 atol = eps(FT)

    #action
    tmp = K1D.aerosol_activation_helper(params..., q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, 0.0, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(FT)
    @test tmp.S_Na ≈ 0.0 atol = eps(FT)

end
