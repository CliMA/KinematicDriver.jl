params = (common_params, thermo_params, air_params, activation_params)

@testset "Precipitation types" begin

    #setup
    toml_dict = CP.create_toml_dict(FT)
    pn = "NoPrecipitation"
    p0m = "Precipitation0M"
    p1m = "Precipitation1M"
    p2m = "Precipitation2M"
    rf_1 = "CliMA_1M"
    rf_2 = "KK2000"
    rf_3 = "B1994"
    rf_4 = "LD2004"
    rf_5 = "TC1980"
    rf_6 = "VarTimeScaleAcnv"
    rf_7 = "SB2006"
    st_1 = "CliMA_1M"
    st_2 = "Chen2022"

    #action
    precip_pn = CO.get_precipitation_type(pn, toml_dict)
    precip_0m = CO.get_precipitation_type(p0m, toml_dict)
    precip_1m_1 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_1)
    precip_1m_2 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_2)
    precip_1m_3 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_3)
    precip_1m_4 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_4)
    precip_1m_5 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_5)
    precip_1m_6 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_6, sedimentation_choice = st_1)
    precip_1m_7 = CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = rf_6, sedimentation_choice = st_2)
    precip_2m = CO.get_precipitation_type(p2m, toml_dict, rain_formation_choice = rf_7)

    #test
    @test_throws Exception CO.get_precipitation_type("_", toml_dict, rain_formation_choice = rf_1)
    @test_throws Exception CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = "_")
    @test_throws Exception CO.get_precipitation_type(p2m, toml_dict, rain_formation_choice = rf_1)
    @test precip_pn isa CO.NoPrecipitation
    @test precip_0m isa CO.Precipitation0M
    @test precip_1m_1 isa CO.Precipitation1M
    @test precip_1m_2 isa CO.Precipitation1M
    @test precip_1m_3 isa CO.Precipitation1M
    @test precip_1m_4 isa CO.Precipitation1M
    @test precip_1m_5 isa CO.Precipitation1M
    @test precip_1m_6 isa CO.Precipitation1M
    @test precip_1m_7 isa CO.Precipitation1M
    @test precip_2m isa CO.Precipitation2M

end

@testset "Moisture type" begin

    #setup
    toml_dict = CP.create_toml_dict(FT)
    eqm = "EquilibriumMoisture"
    neqm = "NonEquilibriumMoisture"

    #action
    moisture_eq = CO.get_moisture_type(eqm, toml_dict)
    moisture_neq = CO.get_moisture_type(neqm, toml_dict)

    #test
    @test CO.EquilibriumMoisture <: CO.AbstractMoistureStyle
    @test CO.NonEquilibriumMoisture <: CO.AbstractMoistureStyle
    @test CO.EquilibriumMoisture_ρθq <: CO.EquilibriumMoisture
    @test CO.EquilibriumMoisture_ρdTq <: CO.EquilibriumMoisture
    @test CO.NonEquilibriumMoisture_ρθq <: CO.NonEquilibriumMoisture
    @test CO.NonEquilibriumMoisture_ρdTq <: CO.NonEquilibriumMoisture

    @test_throws Exception CO.get_moisture_type("_", toml_dict)
    @test moisture_eq isa CO.EquilibriumMoisture_ρdTq
    @test moisture_neq isa CO.NonEquilibriumMoisture_ρdTq

end

@testset "Initialise state" begin

    @test_throws Exception CO.initialise_state(CO.AbstractMoistureStyle(), CO.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0, N_liq = 0, N_rai = 0, N_aer = 0)

    state = CO.initialise_state(equil_moist_ρθq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist_ρθq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist_ρθq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(equil_moist_ρθq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(equil_moist_ρdTq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist_ρdTq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist_ρdTq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(equil_moist_ρdTq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(nequil_moist_ρθq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist_ρθq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist_ρθq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(nequil_moist_ρθq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(nequil_moist_ρdTq, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist_ρdTq, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist_ρdTq, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(nequil_moist_ρdTq, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0
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
        S_Na_activation = [0.0, 0.0],
        S_Nl_activation = [0.0, 0.0],
        term_vel_rai = [0.0, 0.0],
        term_vel_sno = [0.0, 0.0],
        term_vel_N_rai = [0.0, 0.0],
    )

    @test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, no_precip)
    @test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, CO.EquilibriumMoisture())
    @test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, CO.NonEquilibriumMoisture())

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist_ρθq)
    @test aux isa NamedTuple
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.aerosol_variables isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist_ρdTq)
    @test aux isa NamedTuple
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.aerosol_variables isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist_ρθq)
    @test aux isa NamedTuple
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.aerosol_variables isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist_ρdTq)
    @test aux isa NamedTuple
    @test aux.moisture_variables isa CC.Fields.FieldVector
    @test aux.precip_variables isa CC.Fields.FieldVector
    @test aux.moisture_sources isa CC.Fields.FieldVector
    @test aux.aerosol_variables isa CC.Fields.FieldVector
    @test aux.precip_sources isa CC.Fields.FieldVector
    @test LA.norm(aux.precip_sources) == 0
    @test LA.norm(aux.moisture_sources) == 0

end

@testset "Zero tendencies, sources tendencies" begin

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
        S_ql_moisture = [1.0, 2.0],
        S_qi_moisture = [3.0, 4.0],
        S_qt_precip = [5.0, 6.0],
        S_ql_precip = [7.0, 8.0],
        S_qi_precip = [9.0, 10.0],
        S_qr_precip = [11.0, 12.0],
        S_qs_precip = [13.0, 14.0],
        S_Nl_precip = [15.0, 16.0],
        S_Nr_precip = [17.0, 18.0],
        S_Na_activation = [19.0, 20.0],
        S_Nl_activation = [21.0, 22.0],
        term_vel_rai = [0.0, 0.0],
        term_vel_sno = [0.0, 0.0],
        term_vel_N_rai = [0.0, 0.0],
    )
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist_ρθq)
    Y = CO.initialise_state(equil_moist_ρθq, no_precip, ip)

    dY = (; ρq_tot = [10.0, 13.0])
    @test_throws Exception CO.zero_tendencies!(CO.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception CO.zero_tendencies!(CO.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)

    dY = (; ρq_tot = [10.0, 13.0])
    CO.zero_tendencies!(equil_moist_ρθq, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist_ρθq)
    Y = CO.initialise_state(nequil_moist_ρθq, precip_1m, ip)
    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
    )
    CO.zero_tendencies!(nequil_moist_ρθq, dY, Y, aux, 1.0)
    CO.zero_tendencies!(precip_1m, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    dY = (;
        ρq_tot = [10.0, 13.0],
        ρq_liq = [5.0, 10.0],
        ρq_ice = [44.0, -42.0],
        ρq_rai = [1.0, 2.0],
        ρq_sno = [1.0, 1.0],
        N_liq = [1e8, 1e7],
        N_rai = [1e4, 1e4],
        N_aer = [1e9, 1e8],
    )
    CO.zero_tendencies!(nequil_moist_ρθq, dY, Y, aux, 1.0)
    CO.zero_tendencies!(precip_2m, dY, Y, aux, 1.0)
    @test LA.norm(dY) == 0

    @test_throws Exception CO.sources_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception CO.sources_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)
    CO.sources_tendency!(nequil_moist_ρθq, dY, Y, aux, 1.0)
    CO.sources_tendency!(precip_2m, dY, Y, aux, 1.0)
    @test dY.ρq_tot == [5.0, 6.0]
    @test dY.ρq_liq == [8.0, 10.0]
    @test dY.ρq_ice == [12.0, 14.0]
    @test dY.ρq_rai == [11.0, 12.0]
    @test dY.ρq_sno == [13.0, 14.0]
    @test dY.N_liq == [36.0, 38.0]
    @test dY.N_rai == [17.0, 18.0]
    @test dY.N_aer == [19.0, 20.0]

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

    tmp = CO.moisture_helper_vars_eq_ρθq(thermo_params, ρq_tot, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = CO.moisture_helper_vars_eq_ρdTq(thermo_params, ρq_tot, ρ_dry, T)
    @test !isnan(tmp.q_tot + tmp.q_liq + tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice + tmp.θ_dry)
    @test tmp.q_tot >= 0.0
    @test tmp.q_liq >= 0.0
    @test tmp.q_ice >= 0.0
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = CO.moisture_helper_vars_neq_ρθq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ, θ_liq_ice)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ_dry .+ tmp.p .+ tmp.T .+ tmp.θ_dry)
    @test tmp.ρ_dry == ρ - ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.T == TD.air_temperature(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, tmp.T, tmp.ρ_dry)

    tmp = CO.moisture_helper_vars_neq_ρdTq(thermo_params, ρq_tot, ρq_liq, ρq_ice, ρ_dry, T)
    @test !isnan(tmp.q_tot .+ tmp.q_liq .+ tmp.q_ice .+ tmp.ρ .+ tmp.p .+ tmp.θ_liq_ice .+ tmp.θ_dry)
    @test tmp.ρ == ρ_dry + ρq_tot
    @test tmp.p == TD.air_pressure(thermo_params, tmp.ts)
    @test tmp.θ_liq_ice == TD.liquid_ice_pottemp(thermo_params, tmp.ts)
    @test tmp.θ_dry == TD.dry_pottemp(thermo_params, T, ρ_dry)

    tmp = CO.moisture_helper_sources(thermo_params, nequil_moist_ρθq, ρ, T, q_tot, q_liq, q_ice)
    @test !isnan(tmp.S_q_liq .+ tmp.S_q_ice)

    tmp = CO.precip_helper_vars(ρq_rai, ρq_sno, ρ)
    @test !isnan(tmp.q_rai .+ tmp.q_sno)
    @test tmp.q_rai == q_rai
    @test tmp.q_sno == q_sno

    dt = 0.5
    ts = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q_tot)
    tmp = CO.precip_helper_sources!(precip_0m, thermo_params, ts, q_tot, q_liq, q_ice, dt)
    @test !isnan(tmp.S_q_tot .+ tmp.S_q_liq .+ tmp.S_q_ice .+ tmp.S_q_rai .+ tmp.S_q_sno)
    @test tmp.S_q_tot == tmp.S_q_liq + tmp.S_q_ice
    @test tmp.S_q_tot == -(tmp.S_q_rai + tmp.S_q_sno)

    tmp = CO.precip_helper_sources!(
        precip_1m,
        common_params,
        thermo_params,
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

    tmp = CO.precip_helper_sources!(
        precip_2m,
        common_params,
        thermo_params,
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

@testset "precompute aux" begin

    _ip = (;
        ρ = 1.0,
        ρ_dry = 0.999,
        θ_liq_ice = 440.0,
        q_tot = 1e-3,
        q_liq = 0.0,
        q_ice = 0.0,
        ρq_tot = 1e-3,
        ρq_liq = 0.0,
        ρq_ice = 0.0,
        ρq_rai = 0.0,
        ρq_sno = 0.0,
        p = 101300.0,
        T = 300.0,
        θ_dry = 440.0,
        q_rai = 0.0,
        q_sno = 0.0,
        N_liq = 0.0,
        N_rai = 0.0,
        N_aer_0 = 0.0,
        N_aer = 0.0,
        S_ql_moisture = 0.0,
        S_qi_moisture = 0.0,
        S_qt_precip = 0.0,
        S_ql_precip = 0.0,
        S_qi_precip = 0.0,
        S_qr_precip = 0.0,
        S_qs_precip = 0.0,
        S_Nl_precip = 0.0,
        S_Nr_precip = 0.0,
        S_Na_activation = 0.0,
        S_Nl_activation = 0.0,
        term_vel_rai = 0.0,
        term_vel_sno = 0.0,
        term_vel_N_rai = 0.0,
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(0),
        CC.Geometry.ZPoint{FT}(1),
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)
    ip = map(coord -> _ip, coord)

    t = 13.0

    # eq
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist_ρθq)
    Y = CO.initialise_state(equil_moist_ρθq, no_precip, ip)
    dY = Y / 10

    @test_throws Exception CO.precompute_aux_thermo!(CO.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception CO.precompute_aux_precip!(CO.AbstractPrecipitationStyle(), dY, Y, aux, t)

    CO.precompute_aux_thermo!(equil_moist_ρθq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    CO.precompute_aux_thermo!(equil_moist_ρdTq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    # Non-eq
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist_ρθq)
    Y = CO.initialise_state(nequil_moist_ρθq, no_precip, ip)
    dY = Y / 10

    CO.precompute_aux_thermo!(nequil_moist_ρθq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ_dry == aux.moisture_variables.ρ .- Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.T == TD.air_temperature.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    CO.precompute_aux_thermo!(nequil_moist_ρdTq, dY, Y, aux, t)
    @test aux.moisture_variables.ρ == aux.moisture_variables.ρ_dry .+ Y.ρq_tot
    @test aux.moisture_variables.p == TD.air_pressure.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_liq_ice == TD.liquid_ice_pottemp.(thermo_params, aux.moisture_variables.ts)
    @test aux.moisture_variables.θ_dry ==
          TD.dry_pottemp.(thermo_params, aux.moisture_variables.T, aux.moisture_variables.ρ_dry)

    aux = CO.initialise_aux(FT, ip, params..., (; dt = 2.0), 0.0, nequil_moist_ρθq)
    Y = CO.initialise_state(nequil_moist_ρθq, precip_2m, ip)
    dY = Y / 10

    CO.precompute_aux_precip!(precip_2m, dY, Y, aux, t)
    tmp = @. CO.precip_helper_sources!(
        precip_2m,
        aux.common_params,
        aux.thermo_params,
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
