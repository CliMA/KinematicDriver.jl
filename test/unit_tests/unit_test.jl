"""
    Some elementary unit tests
"""

using Test

import LinearAlgebra

import CLIMAParameters
import ClimaCore
import Thermodynamics

include("../../src/Kinematic1D.jl")
include("../create_parameters.jl")

const LA = LinearAlgebra
const CP = CLIMAParameters
const CC = ClimaCore
const TD = Thermodynamics

const KiD = Kinematic1D
const KP = KiD.Parameters

# Create all the params boxes and overwrite the defaults to match PySDM
toml_dict = CP.create_toml_dict(Float64; dict_type = "alias")
params = create_parameter_set(@__DIR__, toml_dict, Float64)
thermo_params = KP.thermodynamics_params(params)

@testset "Moisture and precipitation types" begin

    @test KiD.EquilibriumMoisture <: KiD.AbstractMoistureStyle
    @test KiD.NonEquilibriumMoisture <: KiD.AbstractMoistureStyle

    @test KiD.NoPrecipitation <: KiD.AbstractPrecipitationStyle
    @test KiD.Precipitation0M <: KiD.AbstractPrecipitationStyle
    @test KiD.Precipitation1M <: KiD.AbstractPrecipitationStyle

end

@testset "Parameter overwrites" begin

    @test KP.R_d(params) == 8.314462618 / 0.02896998
    @test KP.R_v(params) == 8.314462618 / 0.018015

    @test KP.MSLP(params) == 100000.0

end

@testset "Make space" begin

    @test KiD.make_function_space(Float64, 0, 100, 10) isa
          Tuple{CC.Spaces.CenterFiniteDifferenceSpace, CC.Spaces.FaceFiniteDifferenceSpace}
end

@testset "Make rhs function" begin

    rhs = KiD.make_rhs_function(KiD.EquilibriumMoisture(), KiD.Precipitation1M())
    @test typeof(rhs) <: Function

end

@testset "Initialise state" begin

    @test_throws Exception KiD.initialise_state(KiD.AbstractMoistureStyle(), KiD.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0)

    state = KiD.initialise_state(KiD.EquilibriumMoisture(), KiD.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KiD.initialise_state(KiD.EquilibriumMoisture(), KiD.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = KiD.initialise_state(KiD.EquilibriumMoisture(), KiD.Precipitation1M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = KiD.initialise_state(KiD.NonEquilibriumMoisture(), KiD.NoPrecipitation(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KiD.initialise_state(KiD.NonEquilibriumMoisture(), KiD.Precipitation0M(), initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = KiD.initialise_state(KiD.NonEquilibriumMoisture(), KiD.Precipitation1M(), initial_profiles)
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
    space, face_space = KiD.make_function_space(Float64, 0, 100, 10)

    @test_throws Exception KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.NoPrecipitation())

    aux = KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.EquilibriumMoisture())
    @test aux isa CC.Fields.FieldVector
    @test aux.constants isa CC.Fields.FieldVector
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test aux.precip_velocities isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_variables) == 0
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

    aux = KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.NonEquilibriumMoisture())
    @test aux isa CC.Fields.FieldVector
    @test aux.constants isa CC.Fields.FieldVector
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test aux.precip_velocities isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_variables) == 0
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

end

@testset "Zero tendencies" begin

    ip = (;
        ρ = [1.0, 1.0],
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
    space, face_space = KiD.make_function_space(Float64, 0, 100, 10)
    aux = KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.EquilibriumMoisture())
    Y = KiD.initialise_state(KiD.EquilibriumMoisture(), KiD.NoPrecipitation(), ip)

    dY = (; ρq_tot = [10.0, 13.0])
    @test_throws Exception KiD.zero_tendencies!(KiD.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception KiD.zero_tendencies!(KiD.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)

    dY = (; ρq_tot = [10.0, 13.0])
    KiD.zero_tendencies!(KiD.EquilibriumMoisture(), dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    aux = KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.NonEquilibriumMoisture())
    Y = KiD.initialise_state(KiD.NonEquilibriumMoisture(), KiD.Precipitation1M(), ip)
    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
    )
    KiD.zero_tendencies!(KiD.NonEquilibriumMoisture(), dY, Y, aux, 1.0)
    KiD.zero_tendencies!(KiD.Precipitation1M(), dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

end

@testset "Tendency helper functions" begin

    ρ = 1.2
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

    tmp = KiD.moisture_helper_vars_eq(params, ρq_tot, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice + tmp.T + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.ts)

    tmp = KiD.moisture_helper_vars_neq(params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.T .+ tmp.θ_dry)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.ts)

    tmp = KiD.moisture_helper_sources(params, ρ, T, q_tot, q_liq, q_ice)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_ice)

    tmp = KiD.precip_helper_vars(ρq_rai, ρq_sno, ρ)
    @test !isnan(tmp.q_rai .+ tmp.q_sno)
    @test tmp.q_rai == q_rai
    @test tmp.q_sno == q_sno

    dt = 0.5
    ts = @. TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    tmp = KiD.precip_helper_sources_0M!(params, ts, q_tot, q_liq, q_ice, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot == tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot == -(tmp.S_q_rai + tmp.S_q_sno)

    tmp = KiD.precip_helper_sources_1M!(params, ts, q_tot, q_liq, q_ice, q_rai, q_sno, T, ρ, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot ≈ tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot ≈ -(tmp.S_q_rai + tmp.S_q_sno)

end

@testset "Errors from precompute_aux, advection_tendency and sources_tendency" begin

    ip = (;
        ρ = [1.0, 1.0],
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
    space, face_space = KiD.make_function_space(Float64, 0, 100, 10)
    aux = KiD.initialise_aux(Float64, ip, params, 0.0, 0.0, face_space, KiD.EquilibriumMoisture())
    Y = KiD.initialise_state(KiD.EquilibriumMoisture(), KiD.NoPrecipitation(), ip)
    dY = (; ρq_tot = [10.0, 13.0])
    t = 13.0

    @test_throws Exception KiD.precompute_aux_thermo!(KiD.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KiD.precompute_aux_precip!(KiD.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception advection_tendency!(KiD.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception advection_tendency!(KiD.AbstractPrecipitationStyle(), dY, Y, aux, t)

    @test_throws Exception KiD.sources_tendency!(KiD.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception KiD.sources_tendency!(KiD.AbstractPrecipitationStyle(), dY, Y, aux, t)

end
