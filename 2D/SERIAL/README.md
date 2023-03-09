## qasMHD/2D/SERIAL
Serial codes for the following two-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability,
- `MRX`... Magnetic reconnection,
- `OTvortex` ... Orszag-Tang vortex problem,
- `RMI` ... Richtmyer-Meshkov instability,
- `RTI` ... Rayleigh-Taylor instability,
- `blast` ... blast wave propagation problem,
- `loop` ... field loop advection problem.

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition,
- `global.hpp` defines the simulation parameters (number of grid points, time step, plasma parameters, etc.),
- `mhd_fd2d.h` defines macros `RMN`, `ODR`, `R_K`, and `CTW` to select the Riemann solvers, spatial and temporal order of accuracy, and multidimensional upwinding for the Constrained Transport method.

### How to run the simulation
```
>cd OTvortex/
>make
>./a.out
```

The result is stored in `dat/`.

### How to check the result
Execute the python script `batch.py`.
```
>python
>>>exec(open("batch.py").read())
```

