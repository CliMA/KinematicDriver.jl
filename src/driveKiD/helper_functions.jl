"""
    Prescribed momentum flux as a function of time
"""
@inline function ρw_helper(t, w1, t1)
    return t < t1 ? w1 * sin(pi * t / t1) : 0.0
end

"""
    Returns the number of new activated aerosol particles and updates aerosol number density
"""
@inline function aerosol_activation_helper(
    kid_params,
    thermo_params,
    air_params,
    activation_params,
    q_tot,
    q_liq,
    N_aer,
    N_aer_0,
    T,
    p,
    ρ,
    ρw,
    dt,
)

    FT = eltype(q_tot)
    S_Nl::FT = FT(0)
    S_Na::FT = FT(0)

    q = TD.PhasePartition(q_tot, q_liq, FT(0))
    S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    if (S < FT(0))
        return (; S_Nl, S_Na)
    end

    (; r_dry, std_dry, κ) = kid_params
    w = ρw / ρ

    aerosol_distribution = CMAM.AerosolDistribution((CMAM.Mode_κ(r_dry, std_dry, N_aer_0, FT(1), FT(1), FT(0), κ, 1),))
    N_act = CMAA.N_activated_per_mode(activation_params, aerosol_distribution, air_params, thermo_params, T, p, w, q)[1]

    if isnan(N_act)
        return (; S_Nl, S_Na)
    end

    S_Nl = max(0, N_act - (N_aer_0 - N_aer)) / dt
    S_Na = -S_Nl

    return (; S_Nl, S_Na)
end
