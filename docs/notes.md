## Jouan 2020 in KiD

### Interpolating data files
1. Using Interpolations.jl and DelimitedFiles.jl
    * There is not a better way of doing this with built-in Julia, however this dependency could be neglected and rewritten 
    * This does allow for the CSV file needed in later work to be seamlessly read in as well, but that file does not have to be a CSV   
### Inital Conditions are a command line option
2. Decided to keep initial conditions in model initation separate from prescribed flow
    * Able to test and debug the two separately, and allows for future reuse if a different data set is desired or a different prescribed flow
    * Unfortunately this has affected all of the model, with the init_param function being taking a new flag in all instances for compilation

![Current Model](contour_plot.png)

### Best Guess on Velocity Field
3. I cannot determine their units, maybe this will resolve when the $\rho$ is implemented in our model 

### Using this notes.md file for updates

4. I can take suggestions. 

    

## open Pr for initial conditions option
### make another pr for documentation about it -> include figures
### make figures about time scaling profile
### m/s in background (velocity max = 4cm/s)
### third PR for velocity
### pr for plotting