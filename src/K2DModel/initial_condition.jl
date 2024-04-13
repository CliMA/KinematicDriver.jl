"""
   Initial condition computation
"""

"""
    Populate the remaining profiles based on the KiD initial condition
    and the density profile
"""
function init_2d_domain(::Type{FT}, common_params, kid_params, thermo_params, ρ_profile, x, z) where {FT}
  
    return K1D.init_1d_column(FT, common_params, kid_params, thermo_params, ρ_profile, z)

end
