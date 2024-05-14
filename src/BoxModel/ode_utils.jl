"""
  Main buliding blocks for the ODE solver
"""

"""
   Interface to ODE solver. Returns the function needed to compute the
   right hand side of the solved ODE. The rhs is assembled via dispatch
   based on the precipitation types.
"""
function make_rhs_function(ms::CO.AbstractMoistureStyle, ps::CO.AbstractPrecipitationStyle)

    function rhs!(dY, Y, aux, t)

        CO.zero_tendencies!(dY)

        CO.precompute_aux_thermo!(ms, Y, aux)
        CO.precompute_aux_precip!(ps, Y, aux)

        CO.precip_sources_tendency!(ms, ps, dY, Y, aux, t)
    end
    return rhs!
end
