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
equil_moist_ρθq = K1D.EquilibriumMoisture_ρθq()
nequil_moist_ρθq = K1D.NonEquilibriumMoisture_ρθq(CMP.CloudLiquid(FT, toml_dict), CMP.CloudIce(FT, toml_dict))
# # ... and precipitation options
no_precip = K1D.NoPrecipitation()
precip_1m = K1D.Precipitation1M(
    CMP.CloudLiquid(FT, toml_dict),
    CMP.CloudIce(FT, toml_dict),
    CMP.Rain(FT, toml_dict),
    CMP.Snow(FT, toml_dict),
    CMP.CollisionEff(FT, toml_dict),
    CMP.KK2000(FT, toml_dict),
    CMP.Blk1MVelType(FT, toml_dict),
)
precip_2m = K1D.Precipitation2M(CMP.SB2006(FT, toml_dict), CMP.SB2006VelType(FT, toml_dict))

@testset "Make space" begin

    @test K2D.make_function_space(FT, xlim = (0, 100.0), zlim = (0, 200.0), helem = 16, velem = 32) isa
          Tuple{CC.Spaces.ExtrudedFiniteDifferenceSpace, CC.Spaces.FaceExtrudedFiniteDifferenceSpace}
end

@testset "Make rhs function" begin

    rhs = K2D.make_rhs_function(equil_moist_ρθq, precip_1m)
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
        S_Na = [0.0, 0.0],
    )
    space, face_space = K2D.make_function_space(FT)

    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, no_precip)
    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, K1D.EquilibriumMoisture())
    @test_throws Exception K2D.initialise_aux(FT, ip, params..., 0.0, 0.0, face_space, K1D.NonEquilibriumMoisture())

    aux = K2D.initialise_aux(FT, ip, params..., 100.0, 200.0, 0.0, 0.0, space, face_space, equil_moist_ρθq)
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

@testset "advection tendency" begin

    space, face_space = K2D.make_function_space(FT)
    coords = CC.Fields.coordinate_field(space)
    face_coords = CC.Fields.coordinate_field(face_space)
    ρ_profile = K1D.ρ_ivp(FT, kid_params, thermo_params)
    init = map(coord -> K2D.init_2d_domain(FT, kid_params, thermo_params, ρ_profile, coord.x, coord.z), coords)
    t = 1.1 * kid_params.t1

    @test_throws Exception K2D.advection_tendency!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    @test_throws Exception K2D.advection_tendency!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    # eq
    aux = K2D.initialise_aux(FT, init, params..., 3000.0, 3000.0, 0.0, 0.0, space, face_space, equil_moist_ρθq)
    K2D.precompute_aux_prescribed_velocity!(aux, t)

    Y = K1D.initialise_state(equil_moist_ρθq, no_precip, init)
    dY = Y / 10
    K2D.advection_tendency!(equil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = K1D.initialise_state(equil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K2D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    # Non-eq
    aux = K2D.initialise_aux(FT, init, params..., 3000.0, 3000.0, 0.0, 0.0, space, face_space, nequil_moist_ρθq)
    K2D.precompute_aux_prescribed_velocity!(aux, t)

    Y = K1D.initialise_state(nequil_moist_ρθq, no_precip, init)
    dY = Y / 10
    K2D.advection_tendency!(nequil_moist_ρθq, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

    Y = K1D.initialise_state(nequil_moist_ρθq, precip_2m, init)
    dY = Y / 10
    K2D.advection_tendency!(precip_2m, dY, Y, aux, t)
    @test dY ≈ Y / 10 atol = eps(FT) * 10

end
