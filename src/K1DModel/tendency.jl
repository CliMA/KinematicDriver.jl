"""
Aerosol activation tendencies
"""

@inline function find_cloud_base(S_Nl, z, cloud_base_S_Nl_and_z)
    # Find S_Nl and z at cloud base:
    z_min = minimum(parent(z)) # height at first level (only works without topography)
    CC.Operators.column_mapreduce!(
        tuple,  # At every point map the input S_Nl and z into a tuple that will be used by reduce
        ((target_S_Nl_value, target_z_value), (S_Nl_value, z_value)) -> ifelse(
            target_z_value == z_min && S_Nl_value >= 1e7,
            (S_Nl_value, z_value),
            (target_S_Nl_value, target_z_value),
        ), # reduction function
        cloud_base_S_Nl_and_z, # destination for output (a field of tuples)
        S_Nl, # input
        z, # input
    )
end

@inline function precompute_aux_activation!(sp::CO.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("activation_tendency not implemented for a given $sp")
end
@inline function precompute_aux_activation!(
    ::Union{CO.NoPrecipitation, CO.Precipitation0M, CO.Precipitation1M},
    dY,
    Y,
    aux,
    t,
) end
@inline function precompute_aux_activation!(::CO.Precipitation2M, dY, Y, aux, t)

    (; thermo_params, activation_params, air_params, kid_params) = aux
    (; T, p, ρ) = aux.thermo_variables
    (; q_tot, q_liq, q_ice, N_aer, N_liq) = aux.microph_variables
    (; r_dry, std_dry, κ) = kid_params
    (; ρw) = aux.prescribed_velocity
    (; dt) = aux.TS

    f_interp = CC.Operators.InterpolateF2C()

    FT = eltype(thermo_params)

    S_liq = aux.scratch.tmp
    N_act = aux.scratch.tmp2
    S_Nl = aux.scratch.tmp3

    cloud_base_S_Nl_and_z = aux.scratch.tmp_surface

    vol_mix_ratio = (FT(1),)
    mass_mix_ratio = (FT(1),)
    molar_mass = (FT(0.42),) # molar mass not needed for 1 aerosol type
    kappa = (κ,)

    S_eltype = eltype(aux.activation_sources)
    to_sources(args...) = S_eltype(tuple(args...))

    # Update N_aer for netcdf output
    @. N_aer = Y.N_aer

    # Compute supersaturation
    @. S_liq = TD.supersaturation(thermo_params, TD.PhasePartition(q_tot, q_liq, q_ice), ρ, T, TD.Liquid())

    # Total number of activated aerosol particles
    @. N_act = CMAA.total_N_activated(
        activation_params,
        CMAM.AerosolDistribution(
            CMAM.Mode_κ(r_dry, std_dry, N_aer, (vol_mix_ratio,), (mass_mix_ratio,), (molar_mass,), (kappa,)),
        ),
        air_params,
        thermo_params,
        T,
        p,
        f_interp(ρw.components.data.:1) / ρ,
        TD.PhasePartition(q_tot, q_liq, q_ice),
    )
    # Convert the total activated number to tendency
    @. S_Nl = ifelse(S_liq < 0 || isnan(N_act), FT(0), max(FT(0), N_act - N_liq) / dt)

    # Use the S_Nl tendnecy at cloud base
    z = CC.Fields.coordinate_field(S_Nl).z # height
    find_cloud_base(S_Nl, z, cloud_base_S_Nl_and_z)

    @. S_Nl = ifelse(z == last(cloud_base_S_Nl_and_z), S_Nl, FT(0))
    @. aux.activation_sources = to_sources(-S_Nl, S_Nl)
end

# helper functions for precomputing aux activation for cloudy
@inline function get_aerosol_activation_rate(
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    q_tot,
    q_liq,
    N_liq,
    N_aer,
    T,
    p,
    ρ,
    ρw,
    dt,
)
    FT = eltype(q_tot)
    S_Nl::FT = FT(0)

    q = TD.PhasePartition(q_tot, q_liq, FT(0))
    S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    (; r_dry, std_dry, κ) = kid_params
    w = ρw / ρ

    aerosol_distribution =
        CMAM.AerosolDistribution((CMAM.Mode_κ(r_dry, std_dry, N_aer, (FT(1),), (FT(1),), (FT(0),), (κ,)),))
    N_act = CMAA.N_activated_per_mode(activation_params, aerosol_distribution, air_params, thermo_params, T, p, w, q)[1]

    S_Nl = ifelse(S < 0 || isnan(N_act), FT(0), max(FT(0), N_act - N_liq) / dt)

    return S_Nl
end
@inline function get_activation_sources(S_Nl, cloudy_params)
    FT = eltype(S_Nl)
    r = FT(2.0e-6) # 2.0 μm
    v = 4 / 3 * π * r^3
    m = v * FT(1000)
    shape = 2
    S_act = ntuple(length(cloudy_params.mom_norms)) do k
        if k == 1
            S_Nl
        elseif k == 2
            S_Nl * m * shape
        elseif k == 3 && cloudy_params.NProgMoms[1] == 3
            S_Nl * m^(k - 1) * shape * (shape + 1)
        else
            FT(0)
        end
    end

    return S_act
end
@inline function precompute_aux_activation!(::CO.CloudyPrecip, dY, Y, aux, t)

    (; kid_params, thermo_params, air_params, activation_params, cloudy_params) = aux
    (; q_tot, q_liq, N_liq, N_aer) = aux.microph_variables
    (; T, p, ρ) = aux.thermo_variables
    (; ρw) = aux.prescribed_velocity
    (; dt) = aux.TS
    S_Nl = aux.scratch.tmp
    cloud_base_S_Nl_and_z = aux.scratch.tmp_surface
    FT = eltype(thermo_params)
    f_interp = CC.Operators.InterpolateF2C()

    @. N_aer = Y.N_aer
    @. S_Nl = get_aerosol_activation_rate(
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        q_tot,
        q_liq,
        N_liq,
        N_aer,
        T,
        p,
        ρ,
        f_interp.(ρw.components.data.:1),
        dt,
    )
    # Use the S_Nl tendnecy at cloud base 
    z = CC.Fields.coordinate_field(S_Nl).z # height
    find_cloud_base(S_Nl, z, cloud_base_S_Nl_and_z)
    @. S_Nl = ifelse(z == last(cloud_base_S_Nl_and_z), S_Nl, FT(0))

    @. aux.activation_sources.N_aer = -1 * S_Nl
    @. aux.activation_sources.activation = get_activation_sources(S_Nl, cloudy_params)
    @. aux.activation_sources.ρq_vap = -aux.activation_sources.activation.:2
end

"""
    Prescribed momentum flux as a function of time
"""
@inline function ρw_helper(t, w1, t1)
    return t < t1 ? w1 * sin(pi * t / t1) : 0.0
end

@inline function precompute_aux_prescribed_velocity!(aux, t)

    FT = eltype(aux.microph_variables.q_tot)
    ρw = FT(ρw_helper(t, aux.kid_params.w1, aux.kid_params.t1))

    @. aux.prescribed_velocity.ρw = CC.Geometry.WVector.(ρw)
    aux.prescribed_velocity.ρw0 = ρw

end

"""
   Advection Equation: ∂ϕ/dt = -∂(vΦ)
"""
@inline function advection_tendency!(sm::CO.AbstractMoistureStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sm")
end
@inline function advection_tendency!(sp::CO.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("advection_tendency not implemented for a given $sp")
end
@inline function advection_tendency!(::CO.EquilibriumMoisture, dY, Y, aux, t)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    @. dY.ρq_tot += -∂(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_tot))

    if Bool(aux.kid_params.qtot_flux_correction)
        fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_tot)
    end

    return dY
end
@inline function advection_tendency!(::CO.CloudyMoisture, dY, Y, aux, t)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    @. dY.ρq_vap += -∂(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_vap))

    if Bool(aux.kid_params.qtot_flux_correction)
        fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
        @. dY.ρq_vap += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_vap)
    end

    return dY
end
@inline function advection_tendency!(::CO.NonEquilibriumMoisture, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()

    ∂_qt = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    ∂_ql = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * FT(0.0))),
        top = CC.Operators.Extrapolate(),
    )
    ∂_qi = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * FT(0.0))),
        top = CC.Operators.Extrapolate(),
    )

    @. dY.ρq_tot += -∂_qt(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_tot))
    @. dY.ρq_liq += -∂_ql(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_liq))
    @. dY.ρq_ice += -∂_qi(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_ice))

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    if Bool(aux.kid_params.qtot_flux_correction)
        @. dY.ρq_tot += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_tot)
    end
    @. dY.ρq_liq += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_liq)
    @. dY.ρq_ice += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_ice)

    return dY
end
@inline function advection_tendency!(::Union{CO.NoPrecipitation, CO.Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::CO.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )

    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )
    @. dY.ρq_sno +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
            ) * If(Y.ρq_sno),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )
    @. dY.ρq_sno += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
        ),
        Y.ρq_sno,
    )


    return dY
end
@inline function advection_tendency!(::CO.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )
    ∂ₐ = CC.Operators.DivergenceF2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    @. dY.N_aer += -∂ₐ((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) * If(Y.N_aer))
    @. dY.N_liq += -∂((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) * If(Y.N_liq))

    @. dY.N_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_N_rai) * FT(-1))
            ) * If(Y.N_rai),
        )
    @. dY.ρq_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_aer += fcc((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)), Y.N_aer)
    @. dY.N_liq += fcc((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)), Y.N_liq)

    @. dY.N_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_N_rai) * FT(-1))
        ),
        Y.N_rai,
    )
    @. dY.ρq_rai += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )

    return dY
end

@inline function advection_tendency!(::CO.CloudyPrecip, dY, Y, aux, t)
    FT = eltype(Y.ρq_vap)
    Nmom = Int(sum(aux.cloudy_params.NProgMoms))

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(FT(0))),
    )
    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    for i in 1:Nmom
        @. dY.moments.:($$i) +=
            -∂(
                (
                    aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                    CC.Geometry.WVector(If(aux.velocities.weighted_vt.:($$i)) * FT(-1))
                ) * If(Y.moments.:($$i)),
            )
        @. dY.moments.:($$i) += fcc(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.weighted_vt.:($$i)) * FT(-1))
            ),
            Y.moments.:($$i),
        )
    end

    return dY
end
