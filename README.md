# KinematicDriver.jl

Kinematic 1D driver

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |

[docs-bld-img]: https://github.com/CliMA/KinematicDriver.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/KinematicDriver.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/KinematicDriver.jl/dev/

[gha-ci-img]: https://github.com/CliMA/KinematicDriver.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/KinematicDriver.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/KinematicDriver.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/KinematicDriver.jl

## Installation and running instructions

KinematicDriver.jl is a Julia registered package.
See the [Project.toml](https://github.com/CliMA/CloudMicrophysics.jl/blob/main/Project.toml)
  for a full list of dependencies.
When using the KinematicDriver.jl,
  the easiest way to obtain the package dependencies
  is to use the Julia built-in package manager
  (accessed by pressing `]` in the Julia REPL):
```bash
julia --project

julia> ]
pkg> instantiate
```
Additional dependencies related to plotting and NetCDF output are
  included in the `test` environment.
See the [Pkg docs](https://docs.julialang.org/en/v1/stdlib/Pkg/)
  for an overview of basic package manager features.

The `KiD_driver.jl` inside `test/experiments/KiD_driver` folder is
  the main file for running the simulations.
It accepts some command line arguments,
  see the `--help` for details.
```bash
julia --project=test/ test/experiments/KiD_driver/KiD_driver.jl --help
```

An example command to run the `Kid_driver.jl` from terminal:
```bash
julia --color=yes --project=test test/experiments/KiD_driver/KiD_driver.jl --moisture_choice=NonEquilibriumMoisture --precipitation_choice=Precipitation1M
```

To run `KiD_driver.jl` interactively, open `julia --project=test` and run:
```julia
include("test/experiments/KiD_driver/parse_commandline.jl")

# Set arguments
config = parse_commandline()
config["FLOAT_TYPE"] = "FLOAT32"
config["moisture_choice"] = "EquilibriumMoisture"

include("test/experiments/KiD_driver/KiD_driver.jl")
```

In addition to simulating the 1D rainshaft, KinematicDriver.jl provides tools for calibrating microphysics parameters against available data. The main program for running calibrations is given inside `test/experiments/calibrations` folder. This program is accompanied by the `config.jl` file that defines all the settings for the dynamics, observations, optimization process, and parameters to be calibrated. To run calibrations of microphysics schemes by using KinematicDriver the `config.jl` file needs to be adjusted. Then the program can be called from terminal:
```bash
julia --color=yes --project=test test/experiments/calibrations/run_calibrations.jl
```
