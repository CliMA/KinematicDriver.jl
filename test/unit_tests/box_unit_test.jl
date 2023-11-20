# override the defaults
toml_dict = CP.create_toml_dict(FT, dict_type = "alias")

# create box params, time stepping and model settings containing initial conditions
box_params = BX.Parameters.BoxModelParameters(FT(1.22), Int(true), Int(true), FT(1e8))
time_stepping = BX.TimeStepping(FT(1), FT(10), FT(100))
model_settings = Dict("init_q_liq" => 1e-3, "init_q_rai" => FT(0), "init_N_liq" => 1e8, "init_N_rai" => FT(0))
# precipitation options
precip_1m = BX.Precipitation1M(
    CMP.CloudLiquid(FT, toml_dict),
    CMP.Rain(FT, toml_dict),
    CMP.CollisionEff(FT, toml_dict),
    CMP.KK2000(FT, toml_dict),
    CMP.Blk1MVelType(FT, toml_dict),
)
precip_2m = BX.Precipitation2M(CMP.SB2006(FT, toml_dict))

@testset "Precipitation types" begin

    #setup
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    p1m = "Precipitation1M"
    p2m = "Precipitation2M"
    rf_1 = "CliMA_1M"
    rf_2 = "KK2000"
    rf_3 = "B1994"
    rf_4 = "LD2004"
    rf_5 = "TC1980"
    rf_6 = "VarTimeScaleAcnv"
    rf_7 = "SB2006"

    #action
    precip_1m_1 = BX.get_precipitation_type(FT, p1m, rf_1, toml_dict)
    precip_1m_2 = BX.get_precipitation_type(FT, p1m, rf_2, toml_dict)
    precip_1m_3 = BX.get_precipitation_type(FT, p1m, rf_3, toml_dict)
    precip_1m_4 = BX.get_precipitation_type(FT, p1m, rf_4, toml_dict)
    precip_1m_5 = BX.get_precipitation_type(FT, p1m, rf_5, toml_dict)
    precip_1m_6 = BX.get_precipitation_type(FT, p1m, rf_6, toml_dict)
    precip_2m = BX.get_precipitation_type(FT, p2m, rf_7, toml_dict)
    
    #setup
    @test BX.Precipitation1M <: BX.AbstractPrecipitationStyle
    @test BX.Precipitation2M <: BX.AbstractPrecipitationStyle

    @test_throws Exception BX.get_precipitation_type(FT, "_", rf_1, toml_dict)
    @test_throws Exception BX.get_precipitation_type(FT, p1m, "_", toml_dict)
    @test_throws Exception BX.get_precipitation_type(FT, p2m, rf_1, toml_dict)
    @test precip_1m_1 isa BX.Precipitation1M
    @test precip_1m_2 isa BX.Precipitation1M
    @test precip_1m_3 isa BX.Precipitation1M
    @test precip_1m_4 isa BX.Precipitation1M
    @test precip_1m_5 isa BX.Precipitation1M
    @test precip_1m_6 isa BX.Precipitation1M
    @test precip_2m isa BX.Precipitation2M
    
end

@testset "Make rhs function" begin

    rhs = BX.make_rhs_function(precip_1m)
    @test typeof(rhs) <: Function
end

@testset "Initialise state" begin

    @test_throws Exception BX.initialise_state(BX.AbstractPrecipitationStyle(), 0)

    state = BX.initialise_state(precip_1m, model_settings)
    @test state isa Vector
    @test length(state) == 2
    @test state[1] ≈ 1e-3
    @test state[2] ≈ FT(0)

    state = BX.initialise_state(precip_2m, model_settings)
    @test state isa Vector
    @test length(state) == 4
    @test state[1] ≈ 1e-3
    @test state[2] ≈ FT(0)
    @test state[3] ≈ 1e8
    @test state[4] ≈ FT(0)

end

@testset "Initialise aux" begin

    aux = BX.initialise_aux(box_params, time_stepping)
    @test aux isa NamedTuple
    @test aux.box_params isa BX.Parameters.BoxModelParameters
    @test aux.TS isa BX.TimeStepping

end

@testset "Tendency helper functions" begin

    q_liq = 1e-3
    q_rai = 1e-3
    N_liq = 1e8
    N_rai = 1e4
    ρ = 1.22
    dt = 1.0

    tmp = BX.precip_helper_sources!(precip_1m, box_params, q_liq, q_rai, ρ, dt)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_rai)
    @test tmp.S_q_liq ≈ -tmp.S_q_rai

    tmp = BX.precip_helper_sources!(precip_2m, box_params, q_liq, q_rai, N_liq, N_rai, ρ, dt)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_rai)
    @test !isnan(tmp.S_N_liq .+ tmp.S_N_rai)
    @test tmp.S_q_liq ≈ -tmp.S_q_rai

end

@testset "precompute aux" begin

    #setup
    aux = BX.initialise_aux(box_params, time_stepping)
    t = 13.0
    Y = BX.initialise_state(precip_1m, model_settings)
    dY = Y ./ 10
    ρ = BX.Parameters.ρ_air(aux.box_params) / (1 - Y[1])
    
    # action
    BX.precompute_aux_precip!(precip_1m, dY, Y, aux, t)
    tmp = BX.precip_helper_sources!(precip_1m, aux.box_params, Y[1], Y[2], ρ, aux.TS.dt)

    # test
    @test_throws Exception BX.precompute_aux_precip!(BX.AbstractPrecipitationStyle(), dY, Y, aux, t)
    @test dY[1] ≈ tmp.S_q_liq atol = eps(FT) * 10
    @test dY[2] ≈ tmp.S_q_rai atol = eps(FT) * 10

    # setup
    Y = BX.initialise_state(precip_2m, model_settings)
    dY = Y / 10

    # action
    BX.precompute_aux_precip!(precip_2m, dY, Y, aux, t)
    tmp = BX.precip_helper_sources!(precip_2m, aux.box_params, Y[1], Y[2], Y[3], Y[4], ρ, aux.TS.dt)

    # test
    @test dY[1] ≈ tmp.S_q_liq atol = eps(FT) * 10
    @test dY[2] ≈ tmp.S_q_rai atol = eps(FT) * 10
    @test dY[3] ≈ tmp.S_N_liq atol = eps(FT) * 10
    @test dY[4] ≈ tmp.S_N_rai atol = eps(FT) * 10
end