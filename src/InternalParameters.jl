# Instantiate CliMA Parameters and overwrite the defaults
struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end

CP.Planet.grav(ps::EarthParameterSet) = ps.nt.grav

CP.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CP.gas_constant(ps::EarthParameterSet) = ps.nt.gas_constant

CP.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
CP.Planet.cp_v(ps::EarthParameterSet) = ps.nt.cp_v

CP.Planet.molmass_dryair(ps::EarthParameterSet) = ps.nt.molmass_dryair
CP.Planet.molmass_water(ps::EarthParameterSet) = ps.nt.molmass_water

CP.Planet.R_d(ps::EarthParameterSet) = ps.nt.gas_constant / ps.nt.molmass_dryair
CP.Planet.R_v(ps::EarthParameterSet) = ps.nt.gas_constant / ps.nt.molmass_water

CP.Microphysics.q_liq_threshold(ps::EarthParameterSet) = ps.nt.q_liq_threshold

nt = (;
    grav = 9.80665,
    MSLP = 100000.0,
    gas_constant = 8.314462618,
    molmass_dryair = 0.02896998,
    molmass_water = 0.018015,
    cp_d = 1005.0,
    cp_v = 1850.0,
    q_liq_threshold = 1e-4,
)

params_overwrite = EarthParameterSet(nt)
