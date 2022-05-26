
"""
This should be turned into a test file
"""
dry=false

include("../src/KiD.jl")
include("data_utils.jl")

const FT = Float64

# Instantiate CliMA Parameters
struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
CP.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CP.gas_constant(ps::EarthParameterSet) = ps.nt.gas_constant
CP.Planet.molmass_dryair(ps::EarthParameterSet) = ps.nt.molmass_dryair
CP.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
nt = (;
    MSLP = 100000.0,
    gas_constant = 8.314462618,
    molmass_dryair = 0.02896998,
    cp_d = 1005.0,
)
params = EarthParameterSet(nt)

# Set up the computational domain and time step
z_min = FT(0)
z_max = FT(2220)
n_elem = 222 # TODO: run PySDM with 20m resolution
Δt = 1.0
Δt_output = 10 * Δt
t_ini = 0.0
t_end = 10.0 * 60

# Updraft momentum flux terms, surface terms and initial conditions
w1 = 2 # m/s * kg/m3
t1 = 600 # s
w_params = (w1 = w1, t1 = t1)
q_surf = 0.016
ρw0 = 0.0

# initialize the timestepping struct
TS = TimeStepping(FT(Δt), FT(Δt_output), FT(t_end))

# create the coordinates,
space, face_space = make_function_space(z_min, z_max, n_elem)

coord = CC.Fields.coordinate_field(space)
face_coord = CC.Fields.coordinate_field(face_space)

# initialize the netcdf output Stats struct
Stats = NetCDFIO_Stats("Output.nc", 1.0, vec(face_coord), vec(coord))

# solve the initial value problem for density profile
ρ_profile = ρ_ivp(FT, params, dry=dry)
# create the initial condition profiles
init = map(coord -> init_1d_column(FT, params, ρ_profile, coord.z, dry=dry), coord)

# create state vector and apply initial condition
Y = initialise_state(EquilibriumMoisture(), NoPrecipitation(), init)

# create aux vector and apply initial condition
aux = initialise_aux(init, params, w_params, q_surf, ρw0, TS, Stats, face_space)

# output the data for plotting
z_centers = parent(CC.Fields.coordinate_field(space))
q_vap = parent(aux.q_tot) - parent(aux.q_liq) - parent(aux.q_ice) 
        - parent(aux.q_rai) - parent(aux.q_sno)
ρ = parent(aux.ρ)
#θ_dry = parent(aux.θ_dry)
#T = parent(aux.T)
θ_dry = parent(init.θ_dry)
T = parent(init.T)
p = parent(aux.p)
q_liq = parent(aux.q_liq)

KM_data = (; z_centers, q_vap, ρ, θ_dry, T, p, q_liq)
if dry
    sdm_case = "dry"
else
    sdm_case = "wet"
end
plot_comparison(KM_data, sdm_case=sdm_case, dir=sdm_case)

