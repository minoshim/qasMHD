## qasMHD/2D/MPI

MPI parallel codes for the following two-dimensional problems are available:

- `KHI`... Kelvin-Helmholtz instability[^1][^2];
- `MRX`... Magnetic reconnection[^3];
- `OTvortex` ... Orszag-Tang vortex problem[^1][^3];
- `RMI` ... Richtmyer-Meshkov instability[^1][^2];
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem[^1][^3];
- `loop` ... field loop advection problem[^3].

### Project structure

- `common/` contains the shared MPI `MHD2D` and dissipative `DMHD2D` classes and solvers, plus the build rules in `common/case.mk`.
- Each problem directory retains its driver (`main.cpp`), initial conditions, parameters, macros, and a small `Makefile`.
- RMI temporarily retains its local `mhd2d_class.cpp` for the inflow boundary; RTI retains `mhd2d_class.cpp` and `mhd2d_solve.cpp` for gravity. Both use the shared `mhd2d_class.hpp`.
- The existing `merge.out`, `batch.py`, and `python/` links remain in each problem directory.

The sources in this directory's `common/` are compiled separately for each problem, using that problem's `mymacros.hpp`. They are not a precompiled class library. The repository-level `common/` and `mpi/` supply the numerical kernels and MPI routines in `libqasmhd.a` and `libmympi.a`.

MRX selects the dissipative sources with `DISSIPATION = 1`. RMI and RTI select their local implementations through `CASE_SRCS` and override `COMMON_NAMES` to avoid duplicate definitions.

### Configuring a problem

- Edit `mhd2d_init_.cpp` for initial conditions, or `dmhd2d_init_.cpp` for MRX (also defines `setdc()`).
- Edit `mhd2d_paras.cpp` for the domain and boundary conditions. `setup_grid(xmin, xmax, ymin, ymax)` takes **global** physical domain bounds excluding ghost cells; it sets local coordinates, mesh spacings, and the initial timestep using the existing MPI decomposition. Optional `xshift, yshift` arguments specify offsets in cell widths: both default to `0.5`; OTvortex uses `0.0, 0.0`.
- The output directory defaults to `./dat/`. To change it, assign `fildir` in `paras()` and create the directory before running.
- Edit `mymacros.hpp` for global mesh size, MPI process counts, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), and `R_K` (1–3) are checked by `static_assert` in the shared class header. `CTW` controls multidimensional CT upwinding[^3].

The case-specific boundary treatments are unchanged. The SERIAL version's RMI/RTI class separation and gravity update fix have not yet been ported. MRX now uses the same diffusion substep safety coefficient of 6 as SERIAL; its time-splitting order is unchanged.

KHI, MRX, RMI, and RTI use `rand_noise_mt()` with `std::mt19937`. Rank 0 broadcasts an unsigned time-based base seed. For reproducible runs, replace `std::time(nullptr)` in the case's initial-condition source with a fixed value such as `10u`, then rebuild. KHI/MRX/RTI seed each generator with `std::seed_seq{seed, mpi_ranx}`, so ranks sharing an X subdomain use identical sequences regardless of their Y coordinate. RMI uses `std::seed_seq{seed, mpi_rank}` for rank-local density noise, and broadcasts interface phases generated on rank 0. RMI retains its generator across the existing inflow calls; the function-local initialization flag and RK-stage inflow timing remain pending the class split. Its generator state is not checkpointed.

Random sequences differ from the legacy generator. Reproducibility requires the same base seed, decomposition, and execution/call history; changing the MPI layout does not preserve the same global random field. Initial perturbation arrays use `std::vector`, with allocation failures stopping all ranks.

### How to run the simulation

From `2D/MPI/`:

```sh
cd OTvortex/
make -C ../../../common
make -C ../../../mpi
make
mkdir -p dat
mpiexec -np 4 env OMP_NUM_THREADS=2 ./a.out
```

Use the MPI C++ compiler `MPICXX` and flags configured in the repository's `Makefile.inc`. The MPI process count must equal `MNP_X*MNP_Y` in `mymacros.hpp`. The example uses four MPI processes, each with two OpenMP threads; `OMP_NUM_THREADS=1` uses one thread per process.

The executable is `a.out`; object and dependency files are stored in each problem's `build/`. Results are written to `dat/`. Local macro and shared-header edits are tracked by dependency files: rerun `make` after editing. Rebuild the repository-level libraries after changing their sources; after changing compiler flags, clean and rebuild both the library objects and the case objects.

Rerunning the same executable and MPI layout loads a saved backup from the output directory when available. Preserve previous results and use a separate directory for a fresh run; see the restart requirements below.

`make clean` removes `a.out` and the current build's object/dependency files, but preserves results. It does not remove objects left directly in a case directory by the previous layout; those files are no longer linked. `make cdata` deletes `dat/*.dat`, including backups; use it deliberately.

### Time stepping and restart

Time integration follows the SERIAL version's bookkeeping. On a new run, the initial CFL timestep is aligned to the output interval with `nrec=max(1,ceil(dtrec/dt))` and `dt=dtrec/nrec`.

- `exec_(0)` keeps this fixed timestep, runs `nrec*N_OUT` steps, and records every `nrec` steps.
- `exec_(1)` retains the existing CFL recalculation every two steps, runs until `tmax`, and shortens the final step to reach `tmax` exactly. Records are written when the time crosses the next output threshold; intermediate record times need not equal that threshold exactly.
- Nonfinite/nonpositive timesteps and unsupported step counts stop all ranks. Common initialization/restart checks are in `common/mhd2d_control.cpp`.

Backup load/save results are collected with `mpi_sync_status()`: all ranks must either find no checkpoint or load one successfully. Missing rank files, invalid metadata, inconsistent rank statuses/times, and I/O errors stop the entire MPI job. A restored timestep is preserved, not realigned; fixed-step checkpoints must already be consistent with the current output interval and counters. Use the same mesh, MPI decomposition, physical parameters, and time-stepping mode/settings for restart. The old format does not encode all of these settings, so validation cannot detect every incompatible configuration.

Output checks use the same path/open/write/close conventions as SERIAL, but invoke `MPI_Abort()` on failure. The MPI file formats are unchanged: normal fields and gravitational potential use single precision, while backup fields use double precision. A new checkpoint is saved only after every rank has completed its normal output.

Backups still overwrite the existing files without retaining generations. This is **not** an atomic multi-file checkpoint: forced termination during a save can leave mixed or incomplete data, and size checks cannot detect every such case. Preserve a known-good checkpoint externally when required. Random-generator state is not stored, so exact restart reproduction of random inflow in RMI is not guaranteed.

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
