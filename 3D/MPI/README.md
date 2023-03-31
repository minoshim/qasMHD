## qasMHD/3D/MPI
MPI parallel codes for the following three-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability;
- `MRX`... Magnetic reconnection;
- `OTvortex` ... Orszag-Tang vortex problem;
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem.

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition;
- `global.hpp` defines the simulation parameters (number of grid points, time step, plasma parameters, etc.);
- `mhd_fd3d.h` defines macros `RMN`, `ODR`, `R_K`, and `CTW` to select the Riemann solvers, spatial and temporal order of accuracy, and the multidimensional upwinding for the Constrained Transport method[^1].

### How to run the simulation
```
>cd OTvortex/
>make
>mpiexec -np 4 -genv OMP_NUM_THREADS 2 ./a.out #for MPICH users
>mpiexec -np 4 -x OMP_NUM_THREADS=2 ./a.out    #for OpenMPI users
```
Here `4` is the number of MPI processes and `2` is the number of OpenMP threads, thus 8 CPU cores are used for this example run.

The number of MPI processes should be equal to the value of `mnp` defined in `global.hpp` (otherwise, the simulation does not run).

Users can abort the run by Ctrl+C, and restart it by the same `a.out` and command.

The result is stored in `dat/`.

### How to check the result
Since the raw simulation data stored in `dat/` are MPI-decomposed, users firstly merge them by
```
>./merge.out dat/ dat/
```

Subsequently, execute the python script `batch.py`.
```
>python
>>>exec(open("batch.py").read())
```

[^1]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)