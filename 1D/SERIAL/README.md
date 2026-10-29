## qasMHD/1D/SERIAL

Serial codes for the following one-dimensional problems are available:

- `shock`... standard shock tube problems[^1];
- `wave`... MHD wave propagation problems[^2];
- `h-shock` ... Hall-MHD shock tube problems;
- `h-wave` ... Hall-MHD wave propagation problems.

### Project structure

- `common/` contains the shared `MHD1D` and Hall-MHD `HMHD1D` classes and solvers, plus the build rules in `common/case.mk`;
- each problem directory retains its driver (`main.cpp`), initial conditions, parameters, macros, and a small `Makefile`;
- `python/` contains visualization and spectrum-analysis scripts linked from each problem directory.

The sources in this directory's `common/` are compiled separately for each problem, using that problem's `mymacros.hpp`. They are not a precompiled class library. The repository-level `common/` supplies the numerical kernels in `libqasmhd.a`.

| Problem | Class | Initial-condition source |
| --- | --- | --- |
| shock, wave | `MHD1D` | `mhd1d_init_.cpp` |
| h-shock, h-wave | `HMHD1D` | `mhd1d_init_.cpp` |

The Hall-MHD problems select their additional shared sources with `HALL := 1` in their Makefiles.

### Configuring a problem

- Edit the initial-condition source listed above for the physical setup.
- Edit `mhd1d_paras.cpp` (`hmhd1d_paras.cpp` for Hall-MHD) for the domain and boundary conditions. `setup_grid(xmin, xmax)` takes physical domain bounds excluding ghost cells and sets the cell-centered coordinates, mesh spacing, and initial timestep.
- The output directory defaults to `./dat/`. To change it, assign `fildir` in `paras()` and create the directory before running.
- Edit `mymacros.hpp` for mesh size, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), and `R_K` (1–3) are checked by `static_assert` in the shared class header.

The wave and h-wave problems use `rand_noise_mt()` with a local `std::mt19937` during initialization. Their seed is fixed at `10` in `mhd1d_init_.cpp`; edit that assignment and rebuild to change it. Sequences differ from the legacy `rand_noise()` generator even for the same seed.

### How to run the simulation

From `1D/SERIAL/`:

```sh
cd shock/
make -C ../../../common
make
mkdir -p dat
OMP_NUM_THREADS=2 ./a.out
```

Use the C++ compiler and flags configured in the repository's `Makefile.inc`. This is a single-process program with OpenMP support; no MPI launcher is needed. `OMP_NUM_THREADS=1` runs with one thread.

The executable is `a.out`; object and dependency files are stored in each problem's `build/`. Results are written to `dat/`. Local macro and shared-header edits are tracked by dependency files: rerun `make` after editing. Rebuild the repository-level library after changing its sources; after changing compiler flags, clean and rebuild both the library objects and the case objects.

`make clean` removes `a.out` and the current build's object/dependency files, but preserves results. `make cdata` deletes `dat/*.dat`; use it deliberately. A new simulation can overwrite existing output files, so preserve previous results first.

### How to check the result

From the problem directory, start Python 3 with NumPy and Matplotlib installed:

```sh
python
```

Then execute the linked script in the interactive session:

```pycon
>>> exec(open("batch.py").read())
```

Enter `dat` when prompted for the data directory, then select an output index. The `batch.py` and `python/` links let these commands run directly from each problem directory. Keeping the Python session open allows further inspection of the loaded data and figures.

To load all output times, execute `exec(open("python/batch_a.py").read())` instead.

For wave spectra, see the scripts and examples in the wave and h-wave problem directories.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../license/COPYING) file for details.

[^1]: [Miyoshi T. and Kusano K. 2005, JCP](https://www.sciencedirect.com/science/article/pii/S0021999105001142?via%3Dihub)
[^2]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
