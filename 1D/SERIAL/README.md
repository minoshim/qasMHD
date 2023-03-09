## qasMHD/1D/SERIAL
Serial codes for the following one-dimensional problems are available:
- `shock`... standard shock tube problems[^1],
- `wave`... MHD wave propagation problems[^2].

Users may edit the following files contained in each directory:
- `init.hpp` defines the initial condition,
- `global.hpp` defines the simulation parameters (number of grid points, time step, plasma parameters, etc.),
- `mhd_fd1d.h` defines macros `RMN`, `ODR`, and `R_K` to select the Riemann solvers, spatial and temporal order of accuracy.

### How to run the simulation
```
>cd shock/
>make
>./a.out
```

The result is stored in `dat/`.

### How to check the result
Execute the python script `batch.py` or `batch_a.py`.
```
>python
>>>exec(open("batch.py").read())
```

[^1]: [Miyoshi T. and Kusano K. 2005, JCP](https://www.sciencedirect.com/science/article/pii/S0021999105001142?via%3Dihub)
[^2]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
