# Define initial condition
function init_1d_column(::Type{FT}, params, z) where {FT}
    # physics parameters
    # R_d::FT = CLIMAParameters.Planet.R_d(params)

    z_0::FT = 0.0
    z_1::FT = 740.0
    z_2::FT = 3260.0
    qv_0::FT = 0.016
    qv_1::FT = 0.0138
    qv_2::FT = 0.0024
    θ_0::FT = 279.9
    θ_1::FT = 279.9
    θ_2::FT = 312.66

    # density, pressure, etc...

    # potential temperature
    θ = z < z_1 ? θ_0 : θ_1 + (θ_2 - θ_1)/(z_2 - z_1) * z

    # water vapour specific humidity (TODO - or is it mixing ratio?)
    qv = z < z_1 ? qv_0 + (qv_1 - qv_0)/(z_1 - z_0) * z : qv_1 + (qv_2 - qv_1)/(z_2 - z_1) * z

    return(θ = θ, qv = qv)
end
