"""
    Some elementary unit tests
"""

using Test

import LinearAlgebra

import CLIMAParameters
import ClimaCore
import Thermodynamics
import Kinematic1D
import CloudMicrophysics

const kid_dir = pkgdir(Kinematic1D)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

const LA = LinearAlgebra
const CP = CLIMAParameters
const CC = ClimaCore
const TD = Thermodynamics
const KID = Kinematic1D
const KP = KID.Parameters
const CMT = CloudMicrophysics.CommonTypes

# Create all the params boxes and overwrite the defaults to match PySDM
toml_dict = CP.create_toml_dict(Float64; dict_type = "alias")
params = create_parameter_set(@__DIR__, toml_dict, Float64)
thermo_params = KP.thermodynamics_params(params)

@testset "Moisture and precipitation types" begin

    @test KID.EquilibriumMoisture <: KID.AbstractMoistureStyle
    @test KID.NonEquilibriumMoisture <: KID.AbstractMoistureStyle
    @test KID.EquilibriumMoisture_ρθq <: KID.EquilibriumMoisture
    @test KID.EquilibriumMoisture_ρdTq <: KID.EquilibriumMoisture
    @test KID.NonEquilibriumMoisture_ρθq <: KID.NonEquilibriumMoisture
    @test KID.NonEquilibriumMoisture_ρdTq <: KID.NonEquilibriumMoisture

    @test KID.NoPrecipitation <: KID.AbstractPrecipitationStyle
    @test KID.Precipitation0M <: KID.AbstractPrecipitationStyle
    @test KID.Precipitation1M <: KID.AbstractPrecipitationStyle
    @test KID.Precipitation2M <: KID.AbstractPrecipitationStyle

end

@testset "Parameter overwrites" begin

    @test KP.R_d(params) == 8.314462618 / 0.02896998
    @test KP.R_v(params) == 8.314462618 / 0.018015

    @test KP.MSLP(params) == 100000.0

end

@testset "Make space" begin

    @test KID.make_function_space(Float64, 0, 100, 10) isa
          Tuple{CC.Spaces.CenterFiniteDifferenceSpace, CC.Spaces.FaceFiniteDifferenceSpace}
end

@testset "Make rhs function" begin

    rhs = KID.make_rhs_function(
        KID.EquilibriumMoisture_ρθq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
    )
    @test typeof(rhs) <: Function

end

@testset "Initialise state" begin

    @test_throws Exception KID.initialise_state(KID.AbstractMoistureStyle(), KID.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0, N_liq = 0, N_rai = 0, N_aer = 0)

    state = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρθq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρθq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρdTq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρdTq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρdTq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρdTq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρθq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρθq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = KID.initialise_state(KID.NonEquilibriumMoisture_ρdTq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KID.initialise_state(KID.NonEquilibriumMoisture_ρdTq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρdTq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρdTq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        initial_profiles,
    )
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
    space, face_space = KID.make_function_space(Float64, 0, 100, 10)

    @test_throws Exception KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.NoPrecipitation())
    @test_throws Exception KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture())
    @test_throws Exception KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture())

    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρθq())
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

    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρdTq())
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

    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
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

    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρdTq())
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
    space, face_space = KID.make_function_space(Float64, 0, 100, 10)
    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρθq())
    Y = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), ip)

    dY = (; ρq_tot = [10.0, 13.0])
    @test_throws Exception KID.zero_tendencies!(KID.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception KID.zero_tendencies!(KID.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)

    dY = (; ρq_tot = [10.0, 13.0])
    KID.zero_tendencies!(KID.EquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    aux = KID.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
    Y = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρθq(),
        KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()),
        ip,
    )
    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
    )
    KID.zero_tendencies!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    KID.zero_tendencies!(KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType()), dY, Y, aux, 1.0)
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
    KID.zero_tendencies!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    KID.zero_tendencies!(KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()), dY, Y, aux, 1.0)
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

    tmp = KID.moisture_helper_vars_eq_ρθq(params, ρq_tot, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = KID.moisture_helper_vars_eq_ρdTq(params, ρq_tot, ρ_dry, T)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = KID.moisture_helper_vars_neq_ρθq(params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T .+ tmp.θ_dry)
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = KID.moisture_helper_vars_neq_ρdTq(params, ρq_tot, ρq_liq, ρq_ice, ρ_dry, T)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice .+ tmp.θ_dry)
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = KID.moisture_helper_sources(params, ρ, T, q_tot, q_liq, q_ice)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_ice)

    tmp = KID.precip_helper_vars(ρq_rai, ρq_sno, ρ)
    @test !isnan(tmp.q_rai .+ tmp.q_sno)
    @test tmp.q_rai == q_rai
    @test tmp.q_sno == q_sno

    dt = 0.5
    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    tmp = KID.precip_helper_sources_0M!(params, ts, q_tot, q_liq, q_ice, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot == tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot == -(tmp.S_q_rai + tmp.S_q_sno)

    ps = KID.Precipitation1M(KID.OneMomentRainFormation(), CMT.Blk1MVelType())
    tmp = KID.precip_helper_sources_1M!(params, ts, q_tot, q_liq, q_ice, q_rai, q_sno, T, ρ, dt, ps)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_ice + tmp.S_q_vap
    @test tmp.S_q_tot ≈ -(tmp.S_q_rai + tmp.S_q_sno)

    ps = KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType())
    tmp = KID.precip_helper_sources_2M!(params, q_tot, q_liq, q_rai, N_liq, N_rai, T, ρ, dt, ps)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_rai)
    @test !isnan(tmp.S_N_liq .+ tmp.S_N_rai)
    @test !isnan(tmp.term_vel_rai .+ tmp.term_vel_N_rai)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_vap
    @test tmp.S_q_tot ≈ -tmp.S_q_rai

end

@testset "advection_tendency and sources_tendency" begin

    space, face_space = KID.make_function_space(Float64, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = KID.ρ_ivp(Float64, params)
    init = map(coord -> KID.init_1d_column(Float64, params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρθq())
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10

    @test_throws Exception advection_tendency!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception advection_tendency!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception KID.sources_tendency!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KID.sources_tendency!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    KID.advection_tendency!(KID.EquilibriumMoisture_ρθq(), dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(Float64) * 10

    Y = KID.initialise_state(
        KID.EquilibriumMoisture_ρθq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        init,
    )
    dY = Y / 10
    KID.advection_tendency!(KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()), dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(Float64) * 10

    # Non-eq
    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
    ρw = 0.0
    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

    Y = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10

    KID.advection_tendency!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(Float64) * 10


    Y = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρθq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        init,
    )
    dY = Y / 10
    KID.advection_tendency!(KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()), dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(Float64) * 10

end

@testset "precompute aux" begin

    space, face_space = KID.make_function_space(Float64, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = KID.ρ_ivp(Float64, params)
    init = map(coord -> KID.init_1d_column(Float64, params, ρ_profile, coord.z), coord)
    t = 13.0

    # eq
    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρθq())

    Y = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10

    @test_throws Exception KID.precompute_aux_thermo!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KID.precompute_aux_precip!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    KID.precompute_aux_thermo!(KID.EquilibriumMoisture_ρθq(), dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    KID.precompute_aux_thermo!(KID.EquilibriumMoisture_ρdTq(), dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    # Non-eq
    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
    Y = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10

    KID.precompute_aux_thermo!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    KID.precompute_aux_thermo!(KID.NonEquilibriumMoisture_ρdTq(), dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    aux = KID.initialise_aux(Float64, init, params, (; dt = 2.0), 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
    Y = KID.initialise_state(
        KID.NonEquilibriumMoisture_ρθq(),
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
        init,
    )
    dY = Y / 10

    KID.precompute_aux_precip!(KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()), dY, Y, aux, t)
    tmp = @. KID.precip_helper_sources_2M!(
        aux.params,
        aux.moisture_variables.q_tot,
        aux.moisture_variables.q_liq,
        aux.precip_variables.q_rai,
        aux.precip_variables.N_liq,
        aux.precip_variables.N_rai,
        aux.moisture_variables.T,
        aux.moisture_variables.ρ,
        aux.TS.dt,
        KID.Precipitation2M(CMT.SB2006Type(), CMT.SB2006VelType()),
    )
    @test aux.precip_sources.S_q_rai ≈ tmp.S_q_rai atol = eps(Float64) * 10
    @test aux.precip_sources.S_q_tot ≈ tmp.S_q_tot atol = eps(Float64) * 10
    @test aux.precip_sources.S_q_liq ≈ tmp.S_q_liq atol = eps(Float64) * 10
    @test aux.precip_sources.S_N_liq ≈ tmp.S_N_liq atol = eps(Float64) * 10
    @test aux.precip_sources.S_N_rai ≈ tmp.S_N_rai atol = eps(Float64) * 10

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
    tmp = KID.aerosol_activation_helper(params, q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, ρw, dt)

    #test
    @test tmp isa NamedTuple
    @test tmp.S_Nl ≈ -tmp.S_Na
    @test 0 < tmp.S_Nl < N_aer_0

    #action
    tmp = KID.aerosol_activation_helper(params, q_tot, q_liq, N_aer, N_aer_0, 290.0, p, ρ, ρw, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(Float64)
    @test tmp.S_Na ≈ 0.0 atol = eps(Float64)

    #action
    tmp = KID.aerosol_activation_helper(params, q_tot, q_liq, N_aer, N_aer_0, T, p, ρ, 0.0, dt)

    #test
    @test tmp.S_Nl ≈ 0.0 atol = eps(Float64)
    @test tmp.S_Na ≈ 0.0 atol = eps(Float64)

end

# calibration pipeline unit tests
include("./calibration_pipeline/run_calibration_pipeline_unit_tests.jl")
