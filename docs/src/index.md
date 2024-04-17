# KinematicDriver.jl

`KinematicDriver.jl` provides tools for prescribed flow simulations for testing and calibrating microphysics schemes.
It uses
  [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl) operators and
  [CloudMicrophysics.jl](https://github.com/CliMA/CloudMicrophysics.jl) tendencies
  to create the numerical problem that is then solved using
  [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Documentation outline

  - **Simulation Setups**: Introduces the physical systems that can be simulated using this package.
  - **Calibration Features**: Describes the tools available for calibrating microphysical parameters, based on comparing simulations provided in this package with reference detailed simulations.
  - **References**: Provides a list of publications that support the tools implemented in this package.


## Contributing and license

`KinematicDriver.jl` is developed by the [Climate Modeling Alliance](https://clima.caltech.edu/) and
is released under the Apache License Version 2.0. Please open an issue or reach out to us if you have any questions or comments.

![Clima logo](assets/logo.svg)