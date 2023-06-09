## qasMHD/3D/MPI
MPI parallel codes for the following three-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability;
- `MRX`... Magnetic reconnection;
- `OTvortex` ... Orszag-Tang vortex problem;
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem;
- `loop` ... field loop advection problem[^1].

Users may edit the following files contained in each directory:
- `mhd3d_init_.cpp` defines the initial condition;
- `mhd3d_paras.cpp` defines the simulation parameters (spatial domain and boundary condition);
- `mymacros.hpp` defines macros about simulation space, MPI number of processes, time, and the solver design (Riemann solver, spatial and temporal accuracies, and the multidimensional upwinding for the Constrained Transport method[^1]).

### How to run the simulation
```
>cd OTvortex/
>make
>mpiexec -np 16 -genv OMP_NUM_THREADS 2 ./a.out #for MPICH users
>mpiexec -np 16 -x OMP_NUM_THREADS=2 ./a.out    #for OpenMPI users
```
Here `16` is the number of MPI processes and `2` is the number of OpenMP threads, thus 32 CPU cores are used for this example run.<br>
The number of MPI processes should be equal to the value of the product of `MNP_X`, `MNP_Y`, and `MPI_Z` defined in `mymacros.hpp` (otherwise, the simulation does not run).

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

### How to visualize in 3D
`python/plt3d.py` module is available, which requires [Mayavi](https://mayavi.readthedocs.io/ja/latest/index.html) (installed by `>pip install mayavi`)
```
>>>from python import plt3d
>>>
>>>plt3d.volume(x,y,z,ro) # Plots iso-surfaces for a 3D volume of data
>>>
>>>fig,obj=plt3d.slice(x,y,z,ro,plane_orientation='x_axes') #Plots 2D image in y-z plane sliced through a 3D volume of data
>>>plt3d.slice(x,y,z,ro,plane_orientation='y_axes',figure=figure) #Overplot in x-z plane
>>>plt3d.slice(x,y,z,ro,plane_orientation='z_axes',figure=figure) #Overplot in x-y plane
```

[^1]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
