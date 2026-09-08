## qasMHD/1D/SERIAL
Serial codes for the following one-dimensional problems are available:

- `shock`... standard shock tube problems[^1];
- `wave`... MHD wave propagation problems[^2];
- `h-shock` ... Hall-MHD shock tube problems;
- `h-wave` ... Hall-MHD wave propagation problems.

### Project structure

- `common/` contains the MHD1D and HMHD1D class and solver implementations shared by all problems;
- `python/` contains shared visualization and spectrum-analysis scripts linked from each problem directory;
- `shock/`, `wave/`, `h-shock/`, and `h-wave/` contain the problem-specific initial conditions, parameters, macros, and build files.

The sources under `common/` are compiled separately for each problem so that its local macro settings are applied.

Users may edit the following files contained in each directory:

- `mhd1d_init_.cpp` defines the initial condition;
- `mhd1d_paras.cpp` (`hmhd1d_paras.cpp` in the Hall-MHD directories) defines the simulation parameters (spatial domain and boundary condition);
- `mymacros.hpp` defines macros about simulation space, time, and the solver design (Riemann solver, spatial and temporal accuracies).

### How to run the simulation
```
>cd shock/
>make
>./a.out
```

The executable is created as `a.out`, while object and dependency files are stored under `build/`.
The result is stored in `dat/`.

Users may run `>make clean` to delete `a.out` and the files generated under `build/`, and `>make cdata` to delete the result stored in `dat/`.

### How to check the result
Execute the python script `batch.py` or `batch_a.py`.
```
>python
>>>exec(open("batch.py").read())
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../license/COPYING) file for details.

[^1]: [Miyoshi T. and Kusano K. 2005, JCP](https://www.sciencedirect.com/science/article/pii/S0021999105001142?via%3Dihub)
[^2]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
