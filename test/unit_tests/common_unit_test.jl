params = (common_params, thermo_params, air_params, activation_params)

TT.@testset "Precipitation types" begin

    #setup
    toml_dict = CP.create_toml_dict(FT)
    pn = "NoPrecipitation"
    p0m = "Precipitation0M"
    p1m = "Precipitation1M"
    p2m = "Precipitation2M"
    pp3 = "PrecipitationP3"
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
    precip_p3 = CO.get_precipitation_type(pp3, toml_dict)

    #test
    TT.@test_throws Exception CO.get_precipitation_type("_", toml_dict, rain_formation_choice = rf_1)
    TT.@test_throws Exception CO.get_precipitation_type(p1m, toml_dict, rain_formation_choice = "_")
    TT.@test_throws Exception CO.get_precipitation_type(p2m, toml_dict, rain_formation_choice = rf_1)
    TT.@test precip_pn isa CO.NoPrecipitation
    TT.@test precip_0m isa CO.Precipitation0M
    TT.@test precip_1m_1 isa CO.Precipitation1M
    TT.@test precip_1m_2 isa CO.Precipitation1M
    TT.@test precip_1m_3 isa CO.Precipitation1M
    TT.@test precip_1m_4 isa CO.Precipitation1M
    TT.@test precip_1m_5 isa CO.Precipitation1M
    TT.@test precip_1m_6 isa CO.Precipitation1M
    TT.@test precip_1m_7 isa CO.Precipitation1M
    TT.@test precip_2m isa CO.Precipitation2M
    TT.@test precip_p3 isa CO.PrecipitationP3

end

TT.@testset "Moisture type" begin

    #setup
    toml_dict = CP.create_toml_dict(FT)
    eqm = "EquilibriumMoisture"
    neqm = "NonEquilibriumMoisture"
    mp3 = "MoistureP3"

    #action
    moisture_eq = CO.get_moisture_type(eqm, toml_dict)
    moisture_neq = CO.get_moisture_type(neqm, toml_dict)
    moisture_p3 = CO.get_moisture_type(mp3, toml_dict)

    #test
    TT.@test CO.EquilibriumMoisture <: CO.AbstractMoistureStyle
    TT.@test CO.NonEquilibriumMoisture <: CO.AbstractMoistureStyle
    TT.@test CO.MoistureP3 <: CO.AbstractMoistureStyle

    TT.@test_throws Exception CO.get_moisture_type("_", toml_dict)
    TT.@test moisture_eq isa CO.EquilibriumMoisture
    TT.@test moisture_neq isa CO.NonEquilibriumMoisture
    TT.@test moisture_p3 isa CO.MoistureP3

end

TT.@testset "Initialise state" begin

    TT.@test_throws Exception CO.initialise_state(CO.AbstractMoistureStyle(), CO.AbstractPrecipitationStyle(), 0)

    mesh = CC.CommonGrids.DefaultZMesh(FT; z_min = 0, z_max = 1, z_elem = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)
    ic_0d = (;
        ρq_tot = FT(0), ρq_liq = FT(0), ρq_ice = FT(0), ρq_rai = FT(0), ρq_sno = FT(0),
        N_liq = FT(0), N_rai = FT(0), N_aer = FT(0), N_ice = FT(0),
        ρq_rim = FT(0), ρq_liqonice = FT(0), B_rim = FT(0),
        q_vap = FT(0), ρq_vap = FT(0), q_rai = FT(0),
    )
    initial_profiles = map(Returns(ic_0d), coord)

    state = CO.initialise_state(equil_moist, no_precip, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist, precip_0m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist, precip_1m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_rai) == 0
    TT.@test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(equil_moist, precip_2m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_rai) == 0
    TT.@test LA.norm(state.N_liq) == 0
    TT.@test LA.norm(state.N_rai) == 0
    TT.@test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(nequil_moist, no_precip, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_liq) == 0
    TT.@test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist, precip_0m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_liq) == 0
    TT.@test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist, precip_1m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_liq) == 0
    TT.@test LA.norm(state.ρq_ice) == 0
    TT.@test LA.norm(state.ρq_rai) == 0
    TT.@test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(nequil_moist, precip_2m, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_liq) == 0
    TT.@test LA.norm(state.ρq_rai) == 0
    TT.@test LA.norm(state.N_liq) == 0
    TT.@test LA.norm(state.N_rai) == 0
    TT.@test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(p3_moist, precip_p3, initial_profiles)
    TT.@test state isa CC.Fields.Field
    TT.@test LA.norm(state.ρq_tot) == 0
    TT.@test LA.norm(state.ρq_liq) == 0
    TT.@test LA.norm(state.ρq_rai) == 0
    TT.@test LA.norm(state.ρq_ice) == 0
    TT.@test LA.norm(state.ρq_rim) == 0
    TT.@test LA.norm(state.ρq_liqonice) == 0
    TT.@test LA.norm(state.N_ice) == 0
    TT.@test LA.norm(state.B_rim) == 0
    TT.@test LA.norm(state.N_liq) == 0
    TT.@test LA.norm(state.N_rai) == 0
    TT.@test LA.norm(state.N_aer) == 0

end


TT.@testset "Tendency helper functions" begin

    _ip = (;
        ρ = 1.2,
        ρ_dry = 1.185,
        θ_liq_ice = 350.0,
        q_tot = 15e-3 / 1.2,
        q_liq = 2e-3 / 1.2,
        q_ice = 2e-3 / 1.2,
        q_rai = 1e-3 / 1.2,
        q_sno = 0.0,
        q_rim = 0.5e-3 / 1.2,
        q_liqonice = 0.5e-3 / 1.2,
        q_vap = 1e-2 / 1.2,
        ρq_tot = 15e-3,
        ρq_liq = 2e-3,
        ρq_ice = 2e-3,
        ρq_rai = 1e-3,
        ρq_sno = 0.0,
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
    mesh = CC.CommonGrids.DefaultZMesh(FT; z_min = 0, z_max = 1, z_elem = 1)
    space = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    coord = CC.Fields.coordinate_field(space)
    ip = map(Returns(_ip), coord)

    #
    # test initialize aux
    #

    TT.@test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, no_precip, no_precip)
    TT.@test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, CO.AbstractMoistureStyle(), no_precip)

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, precip_0m)
    TT.@test aux isa NamedTuple
    TT.@test aux.thermo_variables isa CC.Fields.Field
    TT.@test aux.microph_variables isa CC.Fields.Field
    TT.@test aux.cloud_sources isa Nothing
    TT.@test LA.norm(aux.precip_sources) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_2m)
    TT.@test aux isa NamedTuple
    TT.@test aux.thermo_variables isa CC.Fields.Field
    TT.@test aux.microph_variables isa CC.Fields.Field
    TT.@test LA.norm(aux.cloud_sources) == 0
    TT.@test LA.norm(aux.precip_sources) == 0

    #
    # test zero tendencies
    #

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, no_precip)
    Y = CO.initialise_state(equil_moist, no_precip, ip)
    dY = similar(Y)
    CO.zero_tendencies!(dY)
    TT.@test LA.norm(dY) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_1m)
    Y = CO.initialise_state(nequil_moist, precip_1m, ip)
    dY = similar(Y)
    CO.zero_tendencies!(dY)
    TT.@test LA.norm(dY) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, p3_moist, precip_p3)
    Y = CO.initialise_state(p3_moist, precip_p3, ip)
    dY = similar(Y)
    CO.zero_tendencies!(dY)
    TT.@test LA.norm(dY) == 0

    TT.@test_throws Exception CO.cloud_sources_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    TT.@test_throws Exception CO.precip_sources_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)

    #
    # test precompute functions
    #

    # test if throws error
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, no_precip)
    Y = CO.initialise_state(equil_moist, no_precip, ip)
    TT.@test_throws Exception CO.precompute_aux_thermo!(
        CO.AbstractMoistureStyle(),
        CO.AbstractPrecipitationStyle(),
        Y,
        aux,
    )
    TT.@test_throws Exception CO.precompute_aux_precip!(CO.AbstractPrecipitationStyle(), Y, aux)

    # test precompute_aux_thermo
    for ms in (equil_moist, nequil_moist)
        Y = CO.initialise_state(ms, no_precip, ip)
        aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, ms, no_precip)
        CO.precompute_aux_thermo!(ms, no_precip, Y, aux)
        (; ρ_dry, p, T, θ_dry, θ_liq_ice, ts, ρ, ρ_dry) = aux.thermo_variables
        (; q_tot, q_liq, q_ice) = aux.microph_variables
        for (el,) in CC.Fields.field_iterator(aux.thermo_variables)
            if el != aux.thermo_variables.ts
                TT.@test all(isfinite, parent(el))
            end
        end
        TT.@test all(≥(0), parent(q_tot))
        TT.@test all(≥(0), parent(q_liq))
        TT.@test all(≥(0), parent(q_ice))

        TT.@test ρ_dry == ρ .- Y.ρq_tot
        TT.@test p == TD.air_pressure.(thermo_params, ts)
        TT.@test T == TD.air_temperature.(thermo_params, ts)
        TT.@test θ_dry == TD.dry_pottemp.(thermo_params, T, ρ_dry)
        #TODO - check Thermodynamics?
        #TT.@test get_value(θ_liq_ice)[1] == TD.liquid_ice_pottemp(thermo_params, get_value(ts)[1])
    end

    # test p3 precompute_aux_thermo
    ms = p3_moist
    Y = CO.initialise_state(ms, precip_p3, ip)
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, ms, precip_p3)
    CO.precompute_aux_thermo!(ms, precip_p3, Y, aux)
    (; ρ_dry, p, T, θ_dry, θ_liq_ice, ts, ρ, ρ_dry) = aux.thermo_variables
    (; q_tot, q_liq, q_ice) = aux.microph_variables
    for (el,) in CC.Fields.field_iterator(aux.thermo_variables)
        if el != aux.thermo_variables.ts
            TT.@test all(isfinite, parent(el))
        end
    end
    TT.@test all(≥(0), parent(q_tot))
    TT.@test all(≥(0), parent(q_liq))
    TT.@test all(≥(0), parent(q_ice))

    TT.@test ρ_dry == ρ .- Y.ρq_tot
    TT.@test p == TD.air_pressure.(thermo_params, ts)
    TT.@test T == TD.air_temperature.(thermo_params, ts)
    TT.@test θ_dry == TD.dry_pottemp.(thermo_params, T, ρ_dry)
    #TODO - check Thermodynamics?
    #TT.@test get_value(θ_liq_ice)[1] == TD.liquid_ice_pottemp(thermo_params, get_value(ts)[1])

    # test precompute_aux_precip
    for ps in (precip_1m, precip_2m)
        Y = CO.initialise_state(equil_moist, ps, ip)
        aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, ps)
        CO.precompute_aux_thermo!(equil_moist, ps, Y, aux)
        CO.precompute_aux_precip!(ps, Y, aux)
        for (el,) in CC.Fields.field_iterator(merge.(aux.velocities, aux.microph_variables))
            TT.@test all(isfinite, parent(el))
            TT.@test all(≥(0), parent(el))
        end
    end

    # TODO - tmp remove P3 from tests
    #Y = CO.initialise_state(p3_moist, precip_p3, ip)
    #aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, p3_moist, precip_p3)
    #CO.precompute_aux_thermo!(p3_moist, precip_p3, Y, aux)
    #CO.precompute_aux_precip!(precip_p3, Y, aux)
    #for el in merge(aux.velocities, aux.microph_variables)
    #    TT.@test all(isfinite, get_value(el))
    #    TT.@test all(get_value(el) .>= FT(0))
    #end

    # test precompute_aux_moisture_sources
    Y = CO.initialise_state(nequil_moist, precip_1m, ip)
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_1m)
    CO.precompute_aux_thermo!(nequil_moist, precip_1m, Y, aux)
    CO.precompute_aux_moisture_sources!(nequil_moist, aux)
    TT.@test all(isfinite, parent(aux.cloud_sources.q_ice))
    TT.@test all(isfinite, parent(aux.cloud_sources.q_liq))

    # test precompute_aux_precip_sources
    for ps in (precip_0m, precip_1m, precip_2m)
        Y = CO.initialise_state(equil_moist, precip_1m, ip)
        TS = CO.TimeStepping(FT(10), FT(10), FT(20))
        aux = CO.initialise_aux(FT, ip, params..., TS, 0.0, equil_moist, ps)
        CO.precompute_aux_thermo!(equil_moist, ps, Y, aux)
        CO.precompute_aux_precip_sources!(ps, aux)
        if ps isa CO.Precipitation0M
            (; q_tot, q_liq, q_ice) = aux.precip_sources
            TT.@test all(isfinite, parent(q_tot))
            TT.@test all(isfinite, parent(q_liq))
            TT.@test all(isfinite, parent(q_ice))
            TT.@test parent(q_tot) == parent(q_liq .+ q_ice)
        elseif ps isa CO.Precipitation1M
            (; q_tot, q_liq, q_ice, q_rai, q_sno) = aux.precip_sources
            TT.@test all(isfinite, parent(q_tot))
            TT.@test all(isfinite, parent(q_liq))
            TT.@test all(isfinite, parent(q_ice))
            TT.@test all(isfinite, parent(q_rai))
            TT.@test all(isfinite, parent(q_sno))
            TT.@test all(iszero, parent(q_tot))
        else
            (; q_tot, q_liq, q_rai, N_aer, N_liq, N_rai) = aux.precip_sources
            TT.@test all(isfinite, parent(q_tot))
            TT.@test all(isfinite, parent(q_liq))
            TT.@test all(isfinite, parent(q_rai))
            TT.@test all(isfinite, parent(N_aer))
            TT.@test all(isfinite, parent(N_liq))
            TT.@test all(isfinite, parent(N_rai))
            TT.@test all(iszero, parent(q_tot))
        end
    end
end
