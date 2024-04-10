"""
  Main buliding blocks for the ODE solver
"""

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE. The rhs is assembled via dispatch
   based on the moisture and precipitation types.
"""
function make_rhs_function(ps::CO.AbstractPrecipitationStyle)
    function rhs!(dY, Y, aux, t)

        zero_tendencies!(ps, dY, Y, aux, t)
        update_density!(Y, aux)
        precompute_aux_precip!(ps, dY, Y, aux, t)
        advection_tendency!(ps, dY, Y, aux, t)
        sources_tendency!(ps, dY, Y, aux, t)

    end
    return rhs!
end

"""
    Interface to ODE solver. Returns the ClimaCore.jl FieldVector with
    ODE solver state variables. The state is created via dispatching
    on different moisture and precipitation types
"""
function initialise_state(sp::CO.AbstractPrecipitationStyle, initial_profiles)
    error("initailisation not implemented for a given $sp")
end
function initialise_state(::Union{CO.NoPrecipitation, CO.Precipitation0M}, initial_profiles)
    return CC.Fields.FieldVector(; ρq_tot = initial_profiles.ρq_tot, ρq_liq = initial_profiles.ρq_liq)
end
function initialise_state(::CO.Precipitation1M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_rai = initial_profiles.ρq_rai,
    )
end
function initialise_state(::CO.Precipitation2M, initial_profiles)
    return CC.Fields.FieldVector(;
        ρq_tot = initial_profiles.ρq_tot,
        ρq_liq = initial_profiles.ρq_liq,
        ρq_rai = initial_profiles.ρq_rai,
        N_liq = initial_profiles.N_liq,
        N_rai = initial_profiles.N_rai,
    )
end

"""
   Interface to ODE solver. It initializes the auxiliary state.
   The auxiliary state is created as a ClimaCore FieldVector
   and passed to ODE solver via the `p` parameter of the ODEProblem.
"""
function initialise_aux(FT, ic, kid_params, TS, Stats, face_space)
    term_vel_rai = CC.Geometry.WVector.(zeros(FT, face_space))

    return (;
        moisture_variables = CC.Fields.FieldVector(; ρ_dry = ic.ρ_dry, ρ = ic.ρ, q_tot = ic.q_tot, q_liq = ic.q_liq),
        precip_variables = CC.Fields.FieldVector(; q_rai = ic.q_rai, N_liq = ic.N_liq, N_rai = ic.N_rai),
        precip_sources = CC.Fields.FieldVector(;
            S_q_tot = ic.S_qt,
            S_q_liq = ic.S_ql,
            S_q_rai = ic.S_qr,
            S_N_liq = ic.S_Nl,
            S_N_rai = ic.S_Nr,
        ),
        precip_velocities = CC.Fields.FieldVector(; term_vel_rai = term_vel_rai, term_vel_N_rai = term_vel_rai),
        kid_params = kid_params,
        Stats = Stats,
        TS = TS,
    )
end
