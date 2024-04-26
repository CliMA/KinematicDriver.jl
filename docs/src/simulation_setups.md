# Simulation setups

## Box model
The Box Model simulation represents a zero-dimensional system that focuses on key processes within thermodynamics and collisional microphysics, designed to capture the dynamics of atmospheric particles under controlled conditions without spatial variation. This model is ideal for studying isolated processes in a simplified environment. It includes microphysical processes such as auto-conversion, which represents the transformation of cloud droplets into raindrops through collision and coalescence, and accretion, which models the capture of smaller cloud droplets by larger raindrops as they fall, crucial for understanding precipitation dynamics.

An example of running a box simulation is provided in the following file:
`test/experiments/box_driver/run_box_simulation.jl`.

## KinematicDriver
The KinematicDriver setup is based on the kinematic framework introduced by [KiD2012](@cite):
Vertical momentum flux is constant with height and varying in time.
Density and temerature profiles are constant and defined by the initial
  condition, which is unsaturated.
As the simulation progresses in time moisture is transported upwards,
  supersaturation grows in the upper part of the domain
  and precipitation is formed.
In the second part of the simulation the vertical momentum flux is switched off,
  leaving only cloud microphysics tendencies acting to change
  the model state.
Below figure shows an example prescribed vertical momentum as a function of time.

```@example example_figure
using Plots
import KinematicDriver.K1DModel as K1D

t_range = range(0, 15 * 60, length=100)
w1 = 2.0
t1 = 600.0
plot(t_range / 60.0, [K1D.ρw_helper(t, w1, t1) for t in t_range], linewidth=3, xlabel="t [min]", ylabel="updraft momentum flux [m/s kg/m3]")
savefig("prescribed_momentum_flux.svg") #hide
```
![](prescribed_momentum_flux.svg)

An example of running a one-dimensional KiD simulation is provided in the following file:
`test/experiments/KiD_driver/KiD_driver.jl`.

## Cloud layer simulations
Cloud Layer Simulations employ a one-dimensional setup to explore atmospheric processes, facilitating a detailed examination of the evolution and interactions within a cloud layer. Initially, the simulation initializes the atmospheric column with a uniform gamma distribution, classifying particles into liquid and rain categories. As the simulation progresses, it incorporates thermodynamics, collisional microphysics, and sedimentation. These elements work together to model the complex behaviors of cloud formation and precipitation. By accounting for the vertical distribution and movement of particles, this simulation offers insights into how processes such as auto-conversion, accretion, and sedimentation influence cloud dynamics and precipitation patterns.

An example of running a cloud layer simulation is provided in the following file:
`test/experiments/KiD_col_sed_driver/run_KiD_col_sed_simulation.jl`.

## Two-dimensional system
The two-dimensional setup models the dynamics of an updraft and a downdraft in moist air, incorporating a range of atmospheric processes including thermodynamics, condensation, aerosol activation, collisional processes, and sedimentation of raindrops. These elements collectively simulate complex cloud and precipitation systems. The velocity profile, which varies across spatial dimensions ``x``, ``z`` and over time ``t``, is determined using the stream function ``\psi``, defined as follows:

```math
\begin{align}
\psi(x, z) = - (\rho w)_0 \frac{W}{\pi} \sin\left(\frac{\pi z}{H}\right)\cos\left(\frac{2\pi x}{W}\right)
\end{align}
```

In this equation, ``W`` denotes the width of the domain, ``H`` represents the height of the domain, and ``(\rho w)_0`` is the updraft amplitude, varying over time according to the following equation:

```math
\begin{align}
(\rho w)_0 = \begin{cases} 
(\rho w)_{max} \sin\left(\frac{\pi t}{t_1}\right), & \text{for } t < t_1, \\
0, & \text{for } t > t_1, 
\end{cases}
\end{align}
```
where ``(\rho w)_{max}`` represents the maximum updraft amplitude. An example of running a two-dimensional simulation is provided in the following file:
`test/experiments/Ki2D_driver/run_kinematic2d_simulation.jl`.

## Changes to CliMA defaults

One of the goals of KinematicDriver.jl
  is to test against [PySDM](https://github.com/atmos-cloud-sim-uj/PySDM)
  and the particle-based implementation of the kinematic model available in
  [PySDM-examples](https://github.com/atmos-cloud-sim-uj/PySDM-examples).
As a result, some CliMA constants were changed from their default values to better match PySDM:


| symbol           |         definition                          | units         | default value      | new value             |
|------------------|---------------------------------------------|---------------|--------------------|-----------------------|
| MSLP             | Mean sea level pressure                     | ``Pa``        | ``101325 ``        | ``10^{5}``            |
| grav             | Gravitational acceleration                  | ``m/s^2``     | ``9.81``           | ``9.80665``           |
| gas_constant     | Universal gas constant                      | ``J/mol/K``   | ``8.3144598``      | ``8.314462618``       |
| kappa_d          | Adiabatic exponent for dry air (2/7)        | ``-``         | ``0.28571428571``  | ``0.2855747338575384``|
| cp_v             | Isobaric specific heat of water vapor       | ``J/kg/K``    | ``1859``           | ``1850``              |
| molmass_dryair   | Molecular mass of dry air                   | ``kg/mol``    | ``0.02897``        | ``0.02896998``        |
| molmass_water    | Molecular mass of water                     | ``kg/mol``    | ``0.01801528``     | ``0.018015``          |
| q_liq_threshold  | Cloud liquid water autoconversion threshold | ``kg/kg``     | ``0.0005``         | ``0.0001``            |
