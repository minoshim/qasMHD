## qasMHD/1D/SERIAL
Sereal codes for the following one-dimensional problems are available:
- `shock`... standard shock tube problems
- `wave`... MHD wave propagation problems

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition
- `global.hpp` defines the simulation parameters (number of grid points, plasma parameters, etc.)
- `mhd_fd1d.h` to choose the Riemann solvers, spatial and temporal order of accuracy

## How to run the simulation?
```
>cd shock/
>make
>./a.out
```

The simulation data is stored in the `dat` directory.
