## qasMHD/2D/SERIAL

Serial codes for the following two-dimensional problems are available:

- `KHI`... Kelvin-Helmholtz instability[^1][^2];
- `MRX`... Magnetic reconnection[^3];
- `OTvortex` ... Orszag-Tang vortex problem[^1][^3];
- `RMI` ... Richtmyer-Meshkov instability[^1][^2];
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem[^1][^3];
- `loop` ... field loop advection problem[^3].

### Project structure

- `common/` contains the shared `MHD2D`, dissipative `DMHD2D`, and gravity-aware `GMHD2D` classes and solvers, plus the build rules in `common/case.mk`;
- each problem directory retains its driver (`main.cpp`), initial conditions, parameters, macros, and a small `Makefile`;
- `python/` contains visualization scripts linked from each problem directory.

The sources in this directory's `common/` are compiled separately for each problem, using that problem's `mymacros.hpp`. They are not a precompiled class library. The repository-level `common/` supplies the numerical kernels in `libqasmhd.a`.

| Problem | Class | Initial-condition source |
| --- | --- | --- |
| KHI, OTvortex, blast, loop | `MHD2D` | `mhd2d_init_.cpp` |
| MRX | `DMHD2D` | `dmhd2d_init_.cpp` (also defines `setdc()`) |
| RMI | Local derived class `RMI2D` | `rmi2d_init_.cpp` |
| RTI | `GMHD2D` | `gmhd2d_init_.cpp` (also defines the RTI-specific `bound()`) |

RMI keeps `rmi2d_class.hpp/.cpp` locally for its inflow boundary. MRX and RTI select their additional shared sources with `DISSIPATION := 1` and `GRAVITY := 1`, respectively, in their Makefiles.

### Configuring a problem

- Edit the initial-condition source listed above for the physical setup.
- Edit `mhd2d_paras.cpp` for the domain and boundary conditions. `setup_grid(xmin, xmax, ymin, ymax)` takes physical domain bounds excluding ghost cells and sets the coordinates, mesh spacings, and initial timestep. Optional `xshift, yshift` arguments specify offsets in cell widths: both default to `0.5`; OTvortex explicitly uses `0.0, 0.0` to preserve its original coordinates.
- The output directory defaults to `./dat/`. To change it, assign `fildir` in `paras()` and create the directory before running.
- Edit `mymacros.hpp` for mesh size, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), and `R_K` (1–3) are checked by `static_assert` in the shared class header. `CTW` controls multidimensional CT upwinding[^3].

KHI, MRX, RMI, and RTI use `rand_noise_mt()` with `std::mt19937` when random perturbations are enabled by `RANDOM`. Their seed is currently time-based. For reproducible runs, replace the seed assignment in the case's initial-condition source with a fixed value such as `unsigned seed=10;`, then rebuild. RMI retains its generator between inflow updates; other cases use it only during initialization. Sequences differ from the legacy `rand_noise()` generator even for the same seed.

### How to run the simulation

From `2D/SERIAL/`:

```sh
cd OTvortex/
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

Both `batch.py` and `batch_a.py` automatically subtract `rho*phi_g` when `g_potential.dat` is present in the selected data directory; without it, they use the ordinary MHD pressure formula. Keep each output set with its matching potential file, and do not leave a stale potential file in a non-gravitating run's directory. See [RTI/README.md](RTI/README.md) for the file format.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
