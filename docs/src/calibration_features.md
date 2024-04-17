# Calibration features

Our package provides tools for calibrating microphysical parameters, using one-dimensional simulations and comparing results with reference particle-based simulations. This process aids in fine-tuning the simulation models to better match detailed simulations of particle-based methods, enhancing the accuracy and reliability of simulation outcomes.

## Calibration Methods

We support several calibration methods, including:

- **Ensemble Kalman Inversion (EKI)**: This method uses an ensemble of model states to update estimates of the state and parameters. It is particularly robust against model errors. [EKI](@cite)
- **Unscented Kalman Inversion (UKI)**: UKI is valuable for exploring parameter uncertainties and correlations. [UKI](@cite)
- **Optimization Tools**: We also integrate optimization tools from `Optim.jl` for those who prefer straightforward optimization approaches.

## Configuration File

To use the calibration tools, users must provide a configuration file detailing various elements of the calibration process:

- **Prior Settings**: Defines the prior distributions of parameters.
- **Observations**: Specifies the reference data against which simulations are compared.
- **Calibration Process**: Outlines the methods and approaches used for calibration.
- **Model Parameters**: Details parameters that affect the model's behavior.
- **Statistics**: Specifies the statistical methods to be applied on the data and simulation results.

## Example and Documentation

An example of the configuration file for a simple microphysics calibration can be found at
`test/experiments/calibrations/config.jl`,
and the program to run the calibration at 
`test/experiments/calibrations/run_calibration.jl`.