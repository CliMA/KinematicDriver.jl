"""
Aerosol activation helper functions
"""

"""
    find_cloud_base!(S_Nl, z, cloud_base_S_Nl_and_z)
        S_Nl,
        z,
        cloud_base_S_Nl_and_z,
    )

Performs a vertical column-wise reduction to identify the first level (from bottom to top) where `S_Nl > 0`, 
indicating the onset of cloud droplet activation (i.e., cloud base). The result is stored in `cloud_base_S_Nl_and_z`
as a tuple `(S_Nl, z)` for each column.

 # Arguments
 - `S_Nl`: Field of cloud droplet number concentration tendencies [m⁻³/s], from an activation parameterization.
 - `z`: Field of geometric height [m] corresponding to each grid point.
 - `cloud_base_S_Nl_and_z`: Output field to store the `(S_Nl, z)` values at cloud base for each column.
"""
@inline function find_cloud_base!(S_Nl, z, cloud_base_S_Nl_and_z)
    # Find S_Nl and z at cloud base:
    CC.Operators.column_reduce!(
        ((target_S_Nl_value, target_z_value), (S_Nl_value, z_value)) -> ifelse(
            target_S_Nl_value == 0 && S_Nl_value > 0,
            (S_Nl_value, z_value),
            (target_S_Nl_value, target_z_value),
        ), # reduction function
        cloud_base_S_Nl_and_z, # destination for output (a field of tuples)
        Base.broadcasted(tuple, S_Nl, z),
    )
end

# Compute activation fraction
@inline function get_aerosol_activation_rate(
    common_params,
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    q_tot,
    q_liq,  # cloud liquid + rain
    q_ice,  # cloud ice + snow
    N_liq,  # cloud liquid + rain
    N_aer,
    T,
    p,
    ρ,
    ρw,
    dt,
)
    FT = eltype(q_tot)
    S_Nl::FT = FT(0)

    S::FT = CM.ThermodynamicsInterface.supersaturation_over_liquid(
        thermo_params, q_tot, q_liq, q_ice, ρ, T,
    )

    if S < 0 || N_aer < eps(FT)
        return S_Nl
    end

    (; r_dry, std_dry, κ) = kid_params
    w = ρw / ρ

    _aerosol_budget = N_aer + (common_params.open_system_activation ? FT(0) : N_liq)
    _preexisting_liquid_particles = common_params.local_activation ? N_liq : FT(0)

    aerosol_distribution =
        CMAM.AerosolDistribution((CMAM.Mode_κ(r_dry, std_dry, _aerosol_budget, (FT(1),), (FT(1),), (FT(0),), (κ,)),))

    args = (
        activation_params,
        aerosol_distribution,
        air_params,
        thermo_params,
        T,
        p,
        w,
        q_tot,
        q_liq,
        q_ice,
        _preexisting_liquid_particles,
        FT(0),
    ) # Assuming no ice particles
    S_max = CMAA.max_supersaturation(args...)
    N_act = CMAA.total_N_activated(args...)

    # Convert the total activated number to tendency
    S_Nl = ifelse(
        isnan(N_act) || (common_params.local_activation && (S_max < S || N_act < N_liq)),
        FT(0),
        (N_act - N_liq) / dt,
    )
    return S_Nl
end

# compute moment tendencies from S_Nl for cloudy
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

"""
Aerosol activation tendencies
"""
@inline function precompute_aux_activation!(sp::CO.AbstractPrecipitationStyle, dY, Y, aux, t)
    error("activation_tendency not implemented for a given $sp")
end
@inline function precompute_aux_activation!(
    ::Union{CO.NoPrecipitation, CO.Precipitation0M, CO.Precipitation1M, CO.PrecipitationP3},
    dY,
    Y,
    aux,
    t,
) end
@inline function precompute_aux_activation!(::CO.Precipitation2M, dY, Y, aux, t)

    (; thermo_params, activation_params, air_params, kid_params, common_params) = aux
    (; q_tot, q_liq, q_ice, q_rai, q_sno, N_aer, N_liq, N_rai) = aux.microph_variables
    (; T, p, ρ) = aux.thermo_variables
    (; ρw) = aux.prescribed_velocity
    (; dt) = aux.TS
    S_Nl = aux.scratch.tmp
    cloud_base_S_Nl_and_z = aux.scratch.tmp_surface
    FT = eltype(thermo_params)
    f_interp = CC.Operators.InterpolateF2C()

    @. N_aer = Y.N_aer
    @. S_Nl = get_aerosol_activation_rate(
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        q_tot,
        q_liq + q_rai,
        q_ice + q_sno,
        N_liq + N_rai,
        N_aer,
        T,
        p,
        ρ,
        f_interp.(ρw.components.data.:1),
        dt,
    )
    if !common_params.local_activation
        # Use the S_Nl tendnecy at cloud base
        z = CC.Fields.coordinate_field(S_Nl).z # height
        find_cloud_base!(S_Nl, z, cloud_base_S_Nl_and_z)
        @. S_Nl = ifelse(z == last(cloud_base_S_Nl_and_z), S_Nl, FT(0))
    end

    @. aux.activation_sources.N_aer = -1 * !common_params.open_system_activation * S_Nl
    @. aux.activation_sources.N_liq = S_Nl
end

@inline function precompute_aux_activation!(::CO.CloudyPrecip, dY, Y, aux, t)

    (; common_params, kid_params, thermo_params, air_params, activation_params, cloudy_params) = aux
    (; q_tot, q_liq, q_ice, q_rai, N_liq, N_rai, N_aer) = aux.microph_variables
    (; T, p, ρ) = aux.thermo_variables
    (; ρw) = aux.prescribed_velocity
    (; dt) = aux.TS
    S_Nl = aux.scratch.tmp
    cloud_base_S_Nl_and_z = aux.scratch.tmp_surface
    FT = eltype(thermo_params)
    f_interp = CC.Operators.InterpolateF2C()

    @. N_aer = Y.N_aer
    @. S_Nl = get_aerosol_activation_rate(
        common_params,
        kid_params,
        thermo_params,
        air_params,
        activation_params,
        q_tot,
        q_liq + q_rai,
        q_ice,
        N_liq + N_rai,
        N_aer,
        T,
        p,
        ρ,
        f_interp.(ρw.components.data.:1),
        dt,
    )
    # Use the S_Nl tendnecy at cloud base
    z = CC.Fields.coordinate_field(S_Nl).z # height
    find_cloud_base!(S_Nl, z, cloud_base_S_Nl_and_z)
    @. S_Nl = ifelse(z == last(cloud_base_S_Nl_and_z), S_Nl, FT(0))

    @. aux.activation_sources.N_aer = -1 * !common_params.open_system_activation * S_Nl
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
@inline function advection_tendency!(::CO.MoistureP3, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)

    # advect aerosols

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(aux.prescribed_velocity.ρw0 * aux.q_surf)),
        top = CC.Operators.Extrapolate(),
    )
    @. dY.N_aer += -∂(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.N_aer))
    @. dY.ρq_vap += -∂(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) * If(Y.ρq_vap))

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. dY.N_aer += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.N_aer)
    @. dY.ρq_vap += fcc(aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ), Y.ρq_vap)

    return dY
end
@inline function advection_tendency!(::Union{CO.NoPrecipitation, CO.Precipitation0M}, dY, Y, aux, t) end
@inline function advection_tendency!(::CO.Precipitation1M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    S₁ = aux.scratch.tmp
    S₂ = aux.scratch.tmp2

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )

    @. S₁ =
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )
    @. S₂ =
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
            ) * If(Y.ρq_sno),
        )

    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())
    @. S₁ += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )
    @. S₂ += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_sno) * FT(-1))
        ),
        Y.ρq_sno,
    )

    @. dY.ρq_tot += S₁ + S₂
    @. dY.ρq_rai += S₁
    @. dY.ρq_sno += S₂

    return dY
end
@inline function advection_tendency!(::CO.Precipitation2M, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    S = aux.scratch.tmp

    If = CC.Operators.InterpolateC2F()
    ∂ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(0.0)),
    )
    # Aerosol flux at the bottom is ρw0 * N_d * (1-q_surf) / ρ_SDP 
    surf_aero_flux::FT = aux.prescribed_velocity.ρw0 * aux.common_params.prescribed_Nd * (1 - aux.q_surf) / FT(1.225)
    ∂ₐ = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.SetValue(CC.Geometry.WVector(surf_aero_flux)),
        top = CC.Operators.Extrapolate(),
    )

    @. dY.N_aer += -∂ₐ((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) * If(Y.N_aer))
    @. dY.N_liq += -∂((aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) * If(Y.N_liq))

    @. dY.N_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_N_rai) * FT(-1))
            ) * If(Y.N_rai),
        )
    @. S =
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
    @. S += fcc(
        (
            aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
            CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
        ),
        Y.ρq_rai,
    )

    @. dY.ρq_tot += S
    @. dY.ρq_rai += S

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

@inline function advection_tendency!(::CO.PrecipitationP3, dY, Y, aux, t)
    FT = eltype(Y.ρq_tot)
    # TODO - flux magnitude at top should agree
    # with ρ * q instead of being only based on q

    # P3 advection (introducing particles through boundary):
    # TODO - change run_KiD_simulation so that temp. variables
    # (_q_flux, _N_flux, etc) are not hard coded in the advection tendency

    # update aux with state variables
    (; ρ) = aux.thermo_variables
    (; ρq_tot, ρq_liq, ρq_rai, ρq_ice, ρq_rim, ρq_liqonice, N_ice, N_rai, N_aer, N_liq, B_rim, ρq_vap, q_rai) =
        aux.microph_variables
    @. ρq_ice = Y.ρq_ice
    @. ρq_rim = Y.ρq_rim
    @. ρq_liqonice = Y.ρq_liqonice
    @. N_ice = Y.N_ice
    @. B_rim = Y.B_rim

    # choose flux characteristics
    (; ice_start, _magnitude, _q_flux, _N_flux, _F_rim, _F_liq, _ρ_r_init) = aux.p3_boundary_condition
    # if we have initial ice signal then do not introduce boundary flux
    if ice_start
        _magnitude = FT(0)
    end

    If = CC.Operators.InterpolateC2F()

    # define divergence operators with different boundary conditions
    # for each state variable
    ∂q_ice = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _q_flux)),
    )
    ∂q_rim = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _F_rim * (1 - _F_liq) * _q_flux)),
    )
    ∂q_liqonice = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _F_liq * _q_flux)),
    )
    ∂B_rim = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _F_rim * (1 - _F_liq) * _q_flux / _ρ_r_init)),
    )
    ∂N_ice = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _N_flux)),
    )

    # apply divergence
    @. dY.ρq_ice += ∂q_ice(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice))) * If(Y.ρq_ice),
    )
    @. dY.ρq_rim += ∂q_rim(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice))) * If(Y.ρq_rim),
    )
    @. dY.B_rim += ∂B_rim(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice))) * If(Y.B_rim),
    )
    @. dY.ρq_liqonice += ∂q_liqonice(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice))) * If(Y.ρq_liqonice),
    )
    @. dY.N_ice += ∂N_ice(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_N_ice))) * If(Y.N_ice),
    )
    fcc = CC.Operators.FluxCorrectionC2C(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

    # apply flux correction
    @. dY.ρq_ice += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice) * FT(-1))),
        Y.ρq_ice,
    )
    @. dY.ρq_rim += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice) * FT(-1))),
        Y.ρq_rim,
    )
    @. dY.ρq_liqonice += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice) * FT(-1))),
        Y.ρq_liqonice,
    )
    @. dY.B_rim += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice) * FT(-1))),
        Y.B_rim,
    )
    @. dY.N_ice += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_N_ice) * FT(-1))),
        Y.N_ice,
    )

    # 2M rain advection (zero boundary flux):

    @. ρq_rai = Y.ρq_rai
    @. ρq_liq = Y.ρq_liq
    @. N_rai = Y.N_rai
    @. N_liq = Y.N_liq
    @. N_aer = Y.N_aer

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
    @. dY.q_rai +=
        -∂(
            (
                aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ) +
                CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1))
            ) * If(Y.ρq_rai),
        )

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

    # advecting q_tot:

    @. ρq_vap = Y.ρq_vap
    @. ρq_tot = ρq_ice + ρq_rai + ρq_vap + ρq_liq

    ∂q_tot = CC.Operators.DivergenceF2C(
        bottom = CC.Operators.Extrapolate(),
        top = CC.Operators.SetValue(CC.Geometry.WVector(_magnitude * _q_flux)),
    )

    @. dY.ρq_tot += ∂q_tot(
        FT(-1) * (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice))) * If(Y.ρq_ice) +
        CC.Geometry.WVector(If(aux.velocities.term_vel_rai)) * If(ρq_rai),
    )

    @. dY.ρq_tot += fcc(
        (aux.prescribed_velocity.ρw / If(aux.thermo_variables.ρ)) +
        (CC.Geometry.WVector(If(aux.velocities.term_vel_ice) * FT(-1))) +
        CC.Geometry.WVector(If(aux.velocities.term_vel_rai) * FT(-1)),
        Y.ρq_tot,
    )

    return dY
end
