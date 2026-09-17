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

- `common/` contains the shared MPI `MHD2D`, dissipative `DMHD2D`, and gravity-aware `GMHD2D` classes and solvers, plus the build rules in `common/case.mk`.
- Each problem directory retains its driver (`main.cpp`), initial conditions, parameters, macros, and a small `Makefile`.
- RMI uses the local derived `RMI2D` class for its inflow boundary, as in SERIAL. RTI uses the shared `GMHD2D`; its initial state and special energy boundary conditions are in `RTI/gmhd2d_init_.cpp`.
- The existing `merge.out`, `batch.py`, and `python/` links remain in each problem directory.

The sources in this directory's `common/` are compiled separately for each problem, using that problem's `mymacros.hpp`. They are not a precompiled class library. The repository-level `common/` and `mpi/` supply the numerical kernels and MPI routines in `libqasmhd.a` and `libmympi.a`.

MRX selects the dissipative sources with `DISSIPATION = 1`; RTI selects the gravity sources with `GRAVITY := 1`, as in SERIAL. RMI lists its derived class and initialization sources in `CASE_SRCS`, alongside the shared MHD implementation.

### Configuring a problem

- Edit `mhd2d_init_.cpp` for initial conditions, `dmhd2d_init_.cpp` for MRX (also defines `setdc()`), `rmi2d_init_.cpp` for RMI (also sets inflow parameters), or `gmhd2d_init_.cpp` for RTI (also defines its special `bound()`).
- Edit `mhd2d_paras.cpp` for the domain and boundary conditions. `setup_grid(xmin, xmax, ymin, ymax)` takes **global** physical domain bounds excluding ghost cells; it sets local coordinates, mesh spacings, and the initial timestep using the existing MPI decomposition. Optional `xshift, yshift` arguments specify offsets in cell widths: both default to `0.5`; OTvortex uses `0.0, 0.0`.
- The output directory defaults to `./dat/`. To change it, assign `fildir` in `paras()` and create the directory before running.
- Edit `mymacros.hpp` for global mesh size, MPI process counts, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), and `R_K` (1–3) are checked by `static_assert` in the shared class header. `CTW` controls multidimensional CT upwinding[^3].

RMI refreshes inflow density once per time step and keeps it through all RK stages, only at the global upper Y boundary. RTI now follows SERIAL's gravity treatment: primitive conversion does not modify stored energy, and gravity sources use the density at the start of each RK stage. Its special boundary correction only applies to the actual conserved state at global Y boundaries. MRX uses the same diffusion substep safety coefficient of 6 as SERIAL; its time-splitting order is unchanged.

KHI, MRX, RMI, and RTI use `rand_noise_mt()` with `std::mt19937`. Rank 0 broadcasts an unsigned time-based base seed. For reproducible runs, replace `std::time(nullptr)` in the case's initial-condition source with a fixed value such as `10u`, then rebuild. KHI/MRX/RTI seed each generator with `std::seed_seq{seed, mpi_ranx}`, so ranks sharing an X subdomain use identical sequences regardless of their Y coordinate. RMI uses `std::seed_seq{seed, mpi_rank}` for rank-local density noise, and broadcasts interface phases generated on rank 0. RMI retains its per-instance generator across time steps; boundary calls within RK stages reuse the cached inflow density. Its generator state is not checkpointed.

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

Build the standalone C++11 merge utility once (from a problem directory):
```sh
make -C ../common
```

The existing `merge.out` links point to this shared executable. It uses `CXX` and the flags in `Makefile.inc`, not an MPI launcher. Its temporary-file handling uses POSIX functions, as supported on Linux/macOS.

Merge rank-local data using either form:
```sh
./merge.out dat             # Read and write in dat (same as: ./merge.out dat dat)
mkdir -p merged
./merge.out dat merged      # Read in dat, write in an existing separate directory
```

Trailing slashes are optional; quote directory names containing spaces. The utility retains the coordinate formats and the eight-field float32 layout of `merge_outdat_XXXXX.dat`, but fields 4 and 5 (`Bx`, `By`, zero-based) now contain **cell-center** values. Other fields, including total energy, are unchanged. It copies `t.dat` and `params.dat` when necessary. Raw rank files are not modified. Stop the simulation before merging; input files must not change during the merge.

When all rank-local `g_potential_XXXXX.dat` files exist, they are also merged into `merge_g_potential.dat` (float32, global `(ny, nx)`, without ghost cells). Both `batch.py` and `python/batch_a.py` read this optional file and subtract `rho*phi_g` when calculating pressure, without changing stored total energy. Missing potential on only some ranks, invalid/nonfinite potential, and incorrect binary sizes are errors. With no rank potentials, the utility reports that gravity correction is disabled; it cannot distinguish a nongravitating run from a gravitational run whose potential files were all lost. It refuses to proceed if an old merged potential exists in the destination but no input potentials exist; use a clean output directory instead of reusing unrelated results.

Malformed metadata, missing rank data, allocation failures and I/O failures return a nonzero exit status. Each output is written to a temporary file and renamed only after successful close. This is not a transaction across all output files or a power-loss guarantee: a failed/interrupted merge can leave a mixture of old and new complete files. Rerun successfully before plotting, preferably into a separate directory.

For each output index, the merger also writes the original CT face fields:

| File | float32 shape | Location |
| --- | --- | --- |
| `merge_outdat_XXXXX.dat` | `(8, ny, nx)` | All fields at cell centers |
| `merge_bx_face_XXXXX.dat` | `(ny, nx+1)` | X faces, including the global right face |
| `merge_by_face_XXXXX.dat` | `(ny+1, nx)` | Y faces, including the global top face |

The center values are `Bx[j,i]=(Bx_face[j,i]+Bx_face[j,i+1])/2` and `By[j,i]=(By_face[j,i]+By_face[j+1,i])/2`. Ghost cells in the rank-local data supply the final right/top faces; at least one ghost cell per direction is required. Shared interfaces are taken from the rank on the positive side. Face values are copied without interpolation. The merger does not impose periodic wrapping at physical boundaries. The two CT companion files add approximately 25% to the eight-field merged storage. Two-point interpolation is an output convention and need not reproduce the solver's higher-order interpolation.

In Python, **`bx, by, bz` are cell-center fields** loaded directly from `merge_outdat`, suitable for ordinary plots and pressure evaluation. `bxct, byct` are the separately loaded staggered fields. `divb` is the second-order CT divergence at all cell centers, with shape `(ny,nx)`. `jz` is the CT curl at interior corners, shape `(ny-1,nx-1)`, with coordinates `xjz, yjz`; outer corners are omitted because tangential ghost data are not included. `batch_a.py` adds a leading time dimension to these arrays. The divergence diagnostic is based on float32 output and second-order differences, not necessarily the solver's full-precision/high-order discrete divergence.

Rebuild `merge.out` and rerun it on raw data before using the updated Python scripts. Old merged fields stored face values under the same eight-field layout and are **not interchangeable** with the new center fields. Missing or wrong-size CT companion files cause an error; retain all three files for each time index. For plotting `jz`, use `xjz, yjz`, not the cell-center coordinate arrays.

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
