## qasMHD/2D/MPI
MPI parallel codes for the following two-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability[^1][^2];
- `MRX`... Magnetic reconnection[^3];
- `OTvortex` ... Orszag-Tang vortex problem[^1][^3];
- `RMI` ... Richtmyer-Meshkov instability[^1][^2];
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem[^1][^3];
- `loop` ... field loop advection problem[^3].

Users may edit the following files contained in each directory:
- `mhd2d_init_.cpp` defines the initial condition;
- `mhd2d_paras.cpp` defines the simulation parameters (spatial domain and boundary condition);
- `mymacros.hpp` defines macros about simulation space, MPI number of processes, time, and the solver design (Riemann solver, spatial and temporal accuracies, and the multidimensional upwinding for the Constrained Transport method[^3]).

### How to run the simulation
```
>cd OTvortex/
>make
>mpiexec -np 4 -genv OMP_NUM_THREADS 2 ./a.out #for MPICH users
>mpiexec -np 4 -x OMP_NUM_THREADS=2 ./a.out    #for OpenMPI users
```
Here `4` is the number of MPI processes and `2` is the number of OpenMP threads, thus 8 CPU cores are used for this example run.<br>
The number of MPI processes should be equal to the value of the product of `MNP_X` and `MNP_Y` defined in `mymacros.hpp` (otherwise, the simulation does not run).

Users can abort the run by Ctrl+C, and restart it by the same `a.out` and command.

The result is stored in `dat/`.

Users may run `>make clean` to delete object files, and `>make cdata` to delete the result stored in `dat/`.

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

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
