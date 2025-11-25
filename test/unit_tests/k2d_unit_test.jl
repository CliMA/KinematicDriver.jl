TT.@testset "Make space" begin

    TT.@test K2D.make_function_space(FT, xlim = (0, 100.0), zlim = (0, 200.0), helem = 16, velem = 32) isa
             Tuple{CC.Spaces.ExtrudedFiniteDifferenceSpace, CC.Spaces.FaceExtrudedFiniteDifferenceSpace}
end

TT.@testset "Make rhs function" begin

    rhs = K2D.make_rhs_function(equil_moist, precip_1m)
    TT.@test typeof(rhs) <: Function
end

TT.@testset "Initialise aux" begin

    space, face_space = K2D.make_function_space(FT)
    coords = CC.Fields.coordinate_field(space)
    face_coords = CC.Fields.coordinate_field(face_space)
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params, "ShipwayHill")
    init = CO.initial_condition_1d.(FT, common_params, kid_params, thermo_params, (ρ_profile,), coords.z, "ShipwayHill")
    t = 1.1 * kid_params.t1

    TT.@test_throws Exception K2D.initialise_aux(FT, init, params..., 0.0, 0.0, face_space, no_precip, "ShipwayHill")
    TT.@test_throws Exception K2D.initialise_aux(
        FT,
        init,
        params...,
        0.0,
        0.0,
        face_space,
        K1D.EquilibriumMoisture(),
        "ShipwayHill",
    )
    TT.@test_throws Exception K2D.initialise_aux(
        FT, init, params..., 0.0, 0.0, face_space, K1D.NonEquilibriumMoisture(), "ShipwayHill",
    )


    aux = K2D.initialise_aux(
        FT,
        init,
        params...,
        100.0,
        200.0,
        0.0,
        0.0,
        space,
        face_space,
        equil_moist,
        precip_1m,
        "ShipwayHill",
    )

    TT.@test aux isa NamedTuple
    TT.@test isnothing(aux.cloud_sources)
    TT.@test LA.norm(aux.precip_sources) == 0
end

TT.@testset "advection tendency" begin

    space, face_space = K2D.make_function_space(FT)
    coords = CC.Fields.coordinate_field(space)
    face_coords = CC.Fields.coordinate_field(face_space)
    ρ_profile = CO.ρ_ivp(FT, kid_params, thermo_params, "ShipwayHill")
    init = CO.initial_condition_1d.(FT, common_params, kid_params, thermo_params, (ρ_profile,), coords.z, "ShipwayHill")
    t = 1.1 * kid_params.t1

    TT.@test_throws Exception K2D.advection_tendency!(K1D.AbstractMoistureStyle(), dY, Y, aux, t)
    TT.@test_throws Exception K2D.advection_tendency!(K1D.AbstractPrecipitationStyle(), dY, Y, aux, t)

    ms_styles = [equil_moist, nequil_moist]
    ps_styles = [no_precip, precip_2m]
    for (ms, ps) in zip(ms_styles, ps_styles)
        aux =
            K2D.initialise_aux(FT, init, params..., 3000.0, 3000.0, 0.0, 0.0, space, face_space, ms, ps, "ShipwayHill")
        K2D.precompute_aux_prescribed_velocity!(aux, t)
        Y = CO.initialise_state(ms, ps, init)
        dY = Y ./ 10
        K2D.advection_tendency!(ms, dY, Y, aux, t)
        K2D.advection_tendency!(ps, dY, Y, aux, t)
        TT.@test dY ≈ Y ./ 10 atol = eps(FT) * 10
    end
end
