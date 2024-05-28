# KinematicDriver.jl

`KinematicDriver.jl` (KiD) is a Julia package that provides tools for prescribed flow simulations, and it is used for testing and calibrating microphysics schemes. It serves as a simple substitute for a dynamical core, acting as a driver for bulk cloud microphysics parameterization schemes. KiD is designed to facilitate the assessment of microphysical parameterizations. By prescribing both the momentum and temperature fields, it prevents any feedback between the dynamics and the microphysics, ensuring that variations in results can only be attributed to microphysics parameterizations. An important characteristic of the KiD model is its computational efficiency, which is crucial for efficient parameter calibrations.
The prescribed flow is described by the equation:
```math
\begin{align}
    \rho w(z,t) = (\rho w)_0 \sin\left(\frac{\pi t}{t_1}\right),\quad 0<t<t_1,
\end{align}
```
where ``t_1`` defines the upper time limit of the updraft, ``\rho`` is the dry air density, ``w`` represents the vertical velocity component and ``(\rho w)_0`` denotes the maximum updraft momentum.

Inside `KinematicDriver.jl` we use [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl) operators and [CloudMicrophysics.jl](https://github.com/CliMA/CloudMicrophysics.jl) tendencies to create the numerical problem that is then solved using [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Documentation outline

  - **Simulation Setups**: Introduces the physical systems that can be simulated using this package.
  - **Calibration Features**: Describes the tools available for calibrating microphysical parameters, based on comparing simulations provided in this package with reference detailed simulations.
  - **References**: Provides a list of publications that support the tools implemented in this package.


## Contributing and license

`KinematicDriver.jl` is developed by the [Climate Modeling Alliance](https://clima.caltech.edu/) and
is released under the Apache License Version 2.0. Please open an issue or reach out to us if you have any questions or comments.

![Clima logo](assets/logo.svg)