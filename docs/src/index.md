# Kinematic1D.jl

The Kinematic1D.jl is a single column, prescribed flow driver for testing microphysics schemes.
It uses
  [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl) operators and
  [CloudMicrophysics.jl](https://github.com/CliMA/CloudMicrophysics.jl) tendencies
  to create the numerical problem that is solved using
  [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).
This is a sandbox for testing cloud microphysics schemes.

## Simulation setup

The modeling setup is based on the kinematic framework introduced by [KiD2012](@cite).
The vertical momentum flux is constant with height and varying in time.
The density and temerature profiles are constant and defined by the initial
  condition, which is unsaturated.
As the simulation progresses the moisture is transported upwards leading to
  supersaturation growth in the upper part of the domain
  that leads to the formation of precipitation.

Below figure shows an example vertical momentum evolution.

```@example example_figure
using Plots

include("../../src/helper_functions.jl")

t_range = range(0, 15 * 50, length=100)
w1 = 2.0
t1 = 600.0

plot(t_range / 60.0, [ρw_helper(t, w1, t1) for t in t_range], linewidth=3, xlabel="t [min]", ylabel="updraft momentum flux [m/s kg/m3]")
savefig("prescribed_momentum_flux.svg") #hide
```
![](prescribed_momentum_flux.svg)

## Changes to CliMA defaults

One of the goals is to test against [PySDM](https://github.com/atmos-cloud-sim-uj/PySDM)
  and their implementation of the kinematic model available in
  [PySDM-examples](https://github.com/atmos-cloud-sim-uj/PySDM-examples).
To facilitate the comparison,
  some CliMA constants were changed from their default values to better match PySDM.


| symbol                |         definition                          | units         | default value      | new value             |
|-----------------------|---------------------------------------------|---------------|--------------------|-----------------------|
|``MSLP``               | Mean sea level pressure                     | ``Pa``        | ``101325 ``        | ``10^{5}``            |
|``grav``               | Gravitational acceleration                  | ``m/s^2``     | ``9.81``           | ``9.80665``           |
|``gas\_constant``      | Universal gas constant                      | ``J/mol/K``   | ``8.3144598``      | ``8.314462618``       |
|``kappa\_d``           | Adiabatic exponent for dry air (2/7)        | ``-``         | ``0.28571428571``  | ``0.2855747338575384``|
|``cp\_v``              | Isobaric specific heat of water vapor       | ``J/kg/K``    | ``1859``           | ``1850``              |
|``molmass\_dryair``    | Molecular mass of dry air                   | ``kg/mol``    | ``0.02897``        | ``0.02896998``        |
|``molmass\_water``     | Molecular mass of water                     | ``kg/mol``    | ``0.01801528``     | ``0.018015``          |
| ``q\_liq\_threshold`` | Cloud liquid water autoconversion threshold | ``kg/kg``     | ``0.0005``         | ``0.0001``            |

## Comparison with PySDM

TODO
