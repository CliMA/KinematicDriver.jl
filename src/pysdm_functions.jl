"""
   Helper thermodynamics functions from PySDM.

   We use them in the initialization step only, to make sure the
   initial condition for (ρ, T, and p) is the same between the two models.
"""

function SDM_ρ_dry_of_ρ(ρ, q_vap)
    FT = eltype(ρ)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ / (1 + r_vap)
end

function SDM_ρ_of_ρ_dry(ρ_dry, q_vap)
    FT = eltype(q_vap)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ_dry * (1 + r_vap)
end

function SDM_ρ_dry(params, p, q_vap, θ_std)
    # θ_std - is the potential temperature (not the liquid ice potential temperature)

    FT = eltype(q_vap)
    r_vap::FT = q_vap / (1 - q_vap)

    molmass_ratio::FT = CP.Planet.molmass_ratio(params)
    R_d::FT = CP.Planet.R_d(params)
    cp_d::FT = CP.Planet.cp_d(params)
    p_0::FT = CP.Planet.MSLP(params)

    return p * (1 - 1 / (1 + 1 / molmass_ratio / r_vap)) / ((p / p_0)^(R_d / cp_d) * R_d * θ_std)
end

function SDM_θ_dry(params, θ, q_vap)
    FT = eltype(q_vap)
    r_vap::FT = q_vap / (1 - q_vap)

    R_d::FT = CP.Planet.R_d(params)
    cp_d::FT = CP.Planet.cp_d(params)
    molmass_ratio::FT = CP.Planet.molmass_ratio(params)

    return θ * (1 + r_vap * molmass_ratio)^(R_d / cp_d)
end

function SDM_T(params, θ_dry, ρ_dry)
    FT = eltype(θ_dry)

    R_d::FT = CP.Planet.R_d(params)
    cp_d::FT = CP.Planet.cp_d(params)
    p_0::FT = CP.Planet.MSLP(params)

    return θ_dry * (ρ_dry * θ_dry / p_0 * R_d)^(R_d / cp_d / (1 - R_d / cp_d))
end

function SDM_p(params, ρ_dry, T, q_vap)
    FT = eltype(T)
    r_vap::FT = q_vap / (1 - q_vap)

    R_d::FT = CP.Planet.R_d(params)
    R_v::FT = CP.Planet.R_v(params)

    return ρ_dry * (1 + r_vap) * (R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)) * T
end
