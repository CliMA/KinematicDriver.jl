
1. Advection tendencies
- How to advect tuples rather than scalars?

2. Precipitation tendencies
- Avoid re-computing kernel tensors
- Avoid redeclaring other properties such as velocity, kernel function

3. Activation tendencies
- Add the activation tendency to the moments somehow...

Clean up
Remove unnecessary variables and other tracers from the initialisation and state vector


Note: onset of activation is around 180 sec