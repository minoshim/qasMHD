## qasMHD/1D/SERIAL
Sereal codes for the following one-dimensional problems are available:
- `shock`... standard shock tube problems[^1]
- `wave`... MHD wave propagation problems[^2]

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition
- `global.hpp` defines the simulation parameters (number of grid points, plasma parameters, etc.)
- `mhd_fd1d.h` defines the Riemann solvers, spatial and temporal order of accuracy

### How to run the simulation
```
>cd shock/
>make
>./a.out
```

The simulation data is stored in the `dat` directory.

### How to check the result
Call the python script `batch.py` or `batch_a.py`
```
>python
>>>exec(open("batch.py").read())
```

[^1]:hoge
[^2]:piyo
