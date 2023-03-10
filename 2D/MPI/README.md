## qasMHD/2D/MPI
Serial codes for the following two-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability[^1],
- `MRX`... Magnetic reconnection[^2],
- `OTvortex` ... Orszag-Tang vortex problem[^1][^2],
- `RMI` ... Richtmyer-Meshkov instability[^1],
- `RTI` ... Rayleigh-Taylor instability,

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition,
- `global.hpp` defines the simulation parameters (number of grid points, time step, plasma parameters, etc.),
- `mhd_fd2d.h` defines macros `RMN`, `ODR`, `R_K`, and `CTW` to select the Riemann solvers, spatial and temporal order of accuracy, and the multidimensional upwinding for the Constrained Transport method.
``
### How to run the simulation
```
>cd OTvortex/
>make
>mpiexec -np 4 -genv OMP_NUM_THREADS 2 ./a.out #for MPICH
>mpiexec -np 4 -x OMP_NUM_THREADS=2 ./a.out #for OpenMPI
```
Here, `4` is the number of MPI processes and `2` is the number of OpenMP threads.

The number of MPI processes should be equal to the value of `mnp=MNP_X*MNP_y` defined in `global.hpp` (otherwise, the simulation does not run).

The result is stored in `dat/`.

### How to check the result
Execute the python script `batch.py`.
```
>python
>>>exec(open("batch.py").read())
```

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
