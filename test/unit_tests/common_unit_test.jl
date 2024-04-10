@testset "Precipitation types" begin

    #setup
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
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
    precip_pn = CO.get_precipitation_type(FT, pn, toml_dict)
    precip_0m = CO.get_precipitation_type(FT, p0m, toml_dict)
    precip_1m_1 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_1)
    precip_1m_2 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_2)
    precip_1m_3 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_3)
    precip_1m_4 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_4)
    precip_1m_5 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_5)
    precip_1m_6 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_6, sedimentation_choice = st_1)
    precip_1m_7 = CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = rf_6, sedimentation_choice = st_2)
    precip_2m = CO.get_precipitation_type(FT, p2m, toml_dict, rain_formation_choice = rf_7)

    #test
    @test_throws Exception CO.get_precipitation_type(FT, "_", toml_dict, rain_formation_choice = rf_1)
    @test_throws Exception CO.get_precipitation_type(FT, p1m, toml_dict, rain_formation_choice = "_")
    @test_throws Exception CO.get_precipitation_type(FT, p2m, toml_dict, rain_formation_choice = rf_1)
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
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    eqm = "EquilibriumMoisture"
    neqm = "NonEquilibriumMoisture"

    #action
    moisture_eq = CO.get_moisture_type(FT, eqm, toml_dict)
    moisture_neq = CO.get_moisture_type(FT, neqm, toml_dict)

    #test
    @test CO.EquilibriumMoisture <: CO.AbstractMoistureStyle
    @test CO.NonEquilibriumMoisture <: CO.AbstractMoistureStyle
    @test CO.EquilibriumMoisture_ρθq <: CO.EquilibriumMoisture
    @test CO.EquilibriumMoisture_ρdTq <: CO.EquilibriumMoisture
    @test CO.NonEquilibriumMoisture_ρθq <: CO.NonEquilibriumMoisture
    @test CO.NonEquilibriumMoisture_ρdTq <: CO.NonEquilibriumMoisture

    @test_throws Exception K1D.get_moisture_and_precipitation_types(FT, "_", toml_dict)
    @test moisture_eq isa CO.EquilibriumMoisture_ρdTq
    @test moisture_neq isa CO.NonEquilibriumMoisture_ρdTq

end