"""
    Some elementary unit tests
"""

using Test

import LinearAlgebra

import CLIMAParameters
import ClimaCore
import Thermodynamics
import Kinematic1D

const kid_dir = pkgdir(Kinematic1D)
include(joinpath(kid_dir, "test", "create_parameters.jl"))

const LA = LinearAlgebra
const CP = CLIMAParameters
const CC = ClimaCore
const TD = Thermodynamics
const KID = Kinematic1D
const KP = KID.Parameters

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

    rhs = KID.make_rhs_function(KID.EquilibriumMoisture_ρθq(), KID.Precipitation1M(KID.OneMomentRainFormation()))
    @test typeof(rhs) <: Function

end

@testset "Initialise state" begin

    @test_throws Exception KID.initialise_state(KID.AbstractMoistureStyle(), KID.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0)

    state = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρθq(),
        KID.Precipitation1M(KID.OneMomentRainFormation()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρdTq(), KID.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(KID.EquilibriumMoisture_ρdTq(), KID.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KID.initialise_state(
        KID.EquilibriumMoisture_ρdTq(),
        KID.Precipitation1M(KID.OneMomentRainFormation()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

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
        KID.Precipitation1M(KID.OneMomentRainFormation()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

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
        KID.Precipitation1M(KID.OneMomentRainFormation()),
        initial_profiles,
    )
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

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
        S_ql_moisture = [0.0, 0.0],
        S_qi_moisture = [0.0, 0.0],
        S_qt_precip = [0.0, 0.0],
        S_ql_precip = [0.0, 0.0],
        S_qi_precip = [0.0, 0.0],
        S_qr_precip = [0.0, 0.0],
        S_qs_precip = [0.0, 0.0],
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
        S_ql_moisture = [0.0, 0.0],
        S_qi_moisture = [0.0, 0.0],
        S_qt_precip = [0.0, 0.0],
        S_ql_precip = [0.0, 0.0],
        S_qi_precip = [0.0, 0.0],
        S_qr_precip = [0.0, 0.0],
        S_qs_precip = [0.0, 0.0],
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
    Y = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.Precipitation1M(KID.OneMomentRainFormation()), ip)
    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
    )
    KID.zero_tendencies!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    KID.zero_tendencies!(KID.Precipitation1M(KID.OneMomentRainFormation()), dY, Y, aux, 1.0)
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
    ts = @. TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    tmp = KID.precip_helper_sources_0M!(params, ts, q_tot, q_liq, q_ice, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot == tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot == -(tmp.S_q_rai + tmp.S_q_sno)


    ps = KID.Precipitation1M(KID.OneMomentRainFormation())
    tmp = KID.precip_helper_sources_1M!(params, ts, q_tot, q_liq, q_ice, q_rai, q_sno, T, ρ, dt, ps)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_ice + tmp.S_q_vap
    @test tmp.S_q_tot ≈ -(tmp.S_q_rai + tmp.S_q_sno)

end

@testset "precompute_aux, advection_tendency and sources_tendency" begin

    space, face_space = KID.make_function_space(Float64, 0, 100, 5)
    coord = CC.Fields.coordinate_field(space)
    ρ_profile = KID.ρ_ivp(Float64, params)
    init = map(coord -> KID.init_1d_column(Float64, params, ρ_profile, coord.z), coord)

    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.EquilibriumMoisture_ρθq())
    Y = KID.initialise_state(KID.EquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10
    t = 13.0

    @test_throws Exception KID.precompute_aux_thermo!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KID.precompute_aux_precip!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception advection_tendency!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception advection_tendency!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception KID.sources_tendency!(KID.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KID.sources_tendency!(KID.AbstractPrecipitationStyle(), dY, Y, aux, t)

    KID.advection_tendency!(KID.EquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    @test dY == Y / 10

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

    aux = KID.initialise_aux(Float64, init, params, 0.0, 0.0, face_space, KID.NonEquilibriumMoisture_ρθq())
    Y = KID.initialise_state(KID.NonEquilibriumMoisture_ρθq(), KID.NoPrecipitation(), init)
    dY = Y / 10
    t = 13.0

    KID.advection_tendency!(KID.NonEquilibriumMoisture_ρθq(), dY, Y, aux, 1.0)
    @test dY == Y / 10

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

end

# calibration pipeline unit tests
include("./calibration_pipeline/run_calibration_pipeline_unit_tests.jl")
