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

    @test_throws Exception CO.get_moisture_type("_", toml_dict)
    @test moisture_eq isa CO.EquilibriumMoisture
    @test moisture_neq isa CO.NonEquilibriumMoisture

end

@testset "Initialise state" begin

    @test_throws Exception CO.initialise_state(CO.AbstractMoistureStyle(), CO.AbstractPrecipitationStyle(), 0)

    initial_profiles = (; ρq_tot = 0, ρq_liq = 0, ρq_ice = 0, ρq_rai = 0, ρq_sno = 0, N_liq = 0, N_rai = 0, N_aer = 0)

    state = CO.initialise_state(equil_moist, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0

    state = CO.initialise_state(equil_moist, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(equil_moist, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.N_liq) == 0
    @test LA.norm(state.N_rai) == 0
    @test LA.norm(state.N_aer) == 0

    state = CO.initialise_state(nequil_moist, no_precip, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist, precip_0m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0

    state = CO.initialise_state(nequil_moist, precip_1m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
    @test LA.norm(state.ρq_ice) == 0
    @test LA.norm(state.ρq_rai) == 0
    @test LA.norm(state.ρq_sno) == 0

    state = CO.initialise_state(nequil_moist, precip_2m, initial_profiles)
    @test state isa CC.Fields.FieldVector
    @test LA.norm(state.ρq_tot) == 0
    @test LA.norm(state.ρq_liq) == 0
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
        zero = [0.0, 0.0],
    )

    @test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, no_precip, no_precip)
    @test_throws Exception CO.initialise_aux(FT, ip, params..., 0.0, 0.0, CO.AbstractMoistureStyle(), no_precip)

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, precip_0m)
    @test aux isa NamedTuple
    @test aux.thermo_variables isa NamedTuple
    @test aux.microph_variables isa NamedTuple
    @test aux.cloud_sources isa Nothing
    @test LA.norm(aux.precip_sources) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_2m)
    @test aux isa NamedTuple
    @test aux.thermo_variables isa NamedTuple
    @test aux.microph_variables isa NamedTuple
    @test LA.norm(aux.cloud_sources) == 0
    @test LA.norm(aux.precip_sources) == 0
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
        zero = [1.3, 4.2],
    )
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, no_precip)
    Y = CO.initialise_state(equil_moist, no_precip, ip)
    dY = similar(Y)
    CO.zero_tendencies!(dY)
    @test LA.norm(dY) == 0

    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_1m)
    Y = CO.initialise_state(nequil_moist, precip_1m, ip)
    dY = similar(Y)
    CO.zero_tendencies!(dY)
    @test LA.norm(dY) == 0

    @test_throws Exception CO.cloud_sources_tendency!(CO.AbstractMoistureStyle(), dY, Y, aux, 1.0)
    @test_throws Exception CO.precip_sources_tendency!(CO.AbstractPrecipitationStyle(), dY, Y, aux, 1.0)
end

@testset "Tendency helper functions" begin

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
        N_aer = 1e10,
        N_aer_0 = 1e10,
        zero = 0.0,
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

    get_value(x) = parent(CC.Fields.field_values(CC.Fields.level(x, 1)))

    # test if throws error
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, no_precip)
    Y = CO.initialise_state(equil_moist, no_precip, ip)
    @test_throws Exception CO.precompute_aux_thermo!(CO.AbstractMoistureStyle(), Y, aux)
    @test_throws Exception CO.precompute_aux_precip!(CO.AbstractPrecipitationStyle(), Y, aux)

    # test precompute_aux_thermo
    for ms in (equil_moist, nequil_moist)
        Y = CO.initialise_state(ms, no_precip, ip)
        aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, ms, no_precip)
        CO.precompute_aux_thermo!(ms, Y, aux)
        (; ρ_dry, p, T, θ_dry, θ_liq_ice, ts, ρ, ρ_dry) = aux.thermo_variables
        (; q_tot, q_liq, q_ice) = aux.microph_variables
        for el in aux.thermo_variables
            if el != aux.thermo_variables.ts
                @test all(isfinite, get_value(el))
            end
        end
        @test get_value(q_tot)[1] >= 0.0
        @test get_value(q_liq)[1] >= 0.0
        @test get_value(q_ice)[1] >= 0.0

        @test ρ_dry == ρ .- Y.ρq_tot
        @test p == TD.air_pressure.(thermo_params, ts)
        @test T == TD.air_temperature.(thermo_params, ts)
        @test θ_dry == TD.dry_pottemp.(thermo_params, T, ρ_dry)
        #TODO - check Thermodynamics?
        #@test get_value(θ_liq_ice)[1] == TD.liquid_ice_pottemp(thermo_params, get_value(ts)[1])
    end

    # test precompute_aux_precip
    for ps in (precip_1m, precip_2m)
        Y = CO.initialise_state(equil_moist, ps, ip)
        aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, equil_moist, ps)
        CO.precompute_aux_thermo!(equil_moist, Y, aux)
        CO.precompute_aux_precip!(ps, Y, aux)
        for el in merge(aux.velocities, aux.microph_variables)
            @test all(isfinite, get_value(el))
            @test all(get_value(el) .>= FT(0))
        end
    end

    # test precompute_aux_moisture_sources
    Y = CO.initialise_state(nequil_moist, precip_1m, ip)
    aux = CO.initialise_aux(FT, ip, params..., 0.0, 0.0, nequil_moist, precip_1m)
    CO.precompute_aux_thermo!(nequil_moist, Y, aux)
    CO.precompute_aux_moisture_sources!(nequil_moist, aux)
    @test all(isfinite, get_value(aux.cloud_sources.q_ice))
    @test all(isfinite, get_value(aux.cloud_sources.q_liq))

    # test precompute_aux_precip_sources
    for ps in (precip_0m, precip_1m, precip_2m)
        Y = CO.initialise_state(equil_moist, precip_1m, ip)
        TS = CO.TimeStepping(FT(10), FT(10), FT(20))
        aux = CO.initialise_aux(FT, ip, params..., TS, 0.0, equil_moist, ps)
        CO.precompute_aux_thermo!(equil_moist, Y, aux)
        CO.precompute_aux_precip_sources!(ps, aux)
        if ps isa CO.Precipitation0M
            @test all(isfinite, get_value(aux.precip_sources.q_tot))
            @test all(isfinite, get_value(aux.precip_sources.q_liq))
            @test all(isfinite, get_value(aux.precip_sources.q_ice))
            @test get_value(aux.precip_sources.q_tot) ==
                  get_value(aux.precip_sources.q_liq) + get_value(aux.precip_sources.q_ice)
        elseif ps isa CO.Precipitation1M
            @test all(isfinite, get_value(aux.precip_sources.q_tot))
            @test all(isfinite, get_value(aux.precip_sources.q_liq))
            @test all(isfinite, get_value(aux.precip_sources.q_ice))
            @test all(isfinite, get_value(aux.precip_sources.q_rai))
            @test all(isfinite, get_value(aux.precip_sources.q_sno))
            @test all(
                isapprox.(
                    get_value(aux.precip_sources.q_tot),
                    -get_value(aux.precip_sources.q_rai) - get_value(aux.precip_sources.q_sno),
                    rtol = sqrt(eps(FT)),
                ),
            )

        else
            @test all(isfinite, get_value(aux.precip_sources.q_tot))
            @test all(isfinite, get_value(aux.precip_sources.q_liq))
            @test all(isfinite, get_value(aux.precip_sources.q_rai))
            @test all(isfinite, get_value(aux.precip_sources.N_aer))
            @test all(isfinite, get_value(aux.precip_sources.N_liq))
            @test all(isfinite, get_value(aux.precip_sources.N_rai))
            @test get_value(aux.precip_sources.q_tot) == -get_value(aux.precip_sources.q_rai)
        end
    end
end
