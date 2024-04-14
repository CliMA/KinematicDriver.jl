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

        for eq_style in [ms, ps]
            CO.zero_tendencies!(eq_style, dY, Y, aux, t)
        end

        CO.precompute_aux_thermo!(ms, dY, Y, aux, t)
        CO.precompute_aux_precip!(ps, dY, Y, aux, t)

        for eq_style in [ms, ps]
            CO.sources_tendency!(eq_style, dY, Y, aux, t)
        end

    end
    return rhs!
end