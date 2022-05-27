# Instantiate CliMA Parameters and overwrite the defaults
struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end

CP.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CP.gas_constant(ps::EarthParameterSet) = ps.nt.gas_constant
CP.Planet.molmass_dryair(ps::EarthParameterSet) = ps.nt.molmass_dryair
CP.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
nt = (; MSLP = 100000.0, gas_constant = 8.314462618, molmass_dryair = 0.02896998, cp_d = 1005.0)

params_overwrite = EarthParameterSet(nt)
