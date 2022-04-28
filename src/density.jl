# Thermodynamic states (TODO)
# - solve for density and pressure as a function of z 
function rhod(T,p,qv)
    return Thermodynamics.air_density(params, T, p, qv)
end

function rhod(θ,qv,rhod0)
    R_d::FT = CLIMAParameters.Planet.R_d(params)
    R_v::FT = CLIMAParameters.Planet.R_v(params)
    c_pv::FT = CLIMAParameters.Planet.cp_v(params)
    c_pd::FT = CLIMAParameters.Planet.cp_d(params)
    g::FT = CLIMAParameters.Planet.grav(params)
    
    function drho_dz(rhod, _, z)
        dql_dz=0
        θ_dry = θ * (1 + qv)^(R_d / c_pd) # TODO: fix this with thermodynamics.jl
        T = θ_dry * ((rhod * θ_dry) / 1000 * R_d)^(R_d/c_pd / (1 - R_d/c_pd))
        p = Thermodynamics.air_pressure_given_θ(params, T, rhod, qv) # TODO: should be moist rhod
        lv = Thermodynamics.latent_heat_vapor(params, T)
        Rq = R_v / (1/qv + 1) + R_d / (1+qv)
        cp = c_pv / (1/qv + 1) + c_pd / (1 + qv)
        rho = Thermodynamics.air_density(params, T, p)

        drhodz = (g / T * rho * (Rq / cp - 1) - p * lv / cp / T^2 * dql_dz) / Rq
        return drhodz
    end

    zspan = (z_min, z_max)
    prob = ODEProblem(drhodz, rhod0, (z_min, z_max))
    sol = solve(problem);

    return z -> sol.u(z)
end


