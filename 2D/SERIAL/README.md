## qasMHD/2D/SERIAL
Serial codes for the following two-dimensional problems are available:
- `KHI`... Kelvin-Helmholtz instability[^1];
- `MRX`... Magnetic reconnection[^2];
- `OTvortex` ... Orszag-Tang vortex problem[^1][^2];
- `RMI` ... Richtmyer-Meshkov instability[^1];
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem[^1][^2];
- `loop` ... field loop advection problem[^2].

Users may edit the following files contained in each directory:
- `mhd2d_init_.cpp` defines the initial condition;
- `mhd2d_paras.cpp` defines the simulation parameters (spatial domain and boundary condition);
- `mymacros.hpp` defines macros about simulation space, time, and the solver design (Riemann solver, spatial and temporal orders, and the multidimensional upwinding for the Constrained Transport method[^2]).

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

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
