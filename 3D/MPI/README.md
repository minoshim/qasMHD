## qasMHD/3D/MPI

MPI parallel codes for the following three-dimensional problems are available:

- `KHI`... Kelvin-Helmholtz instability;
- `MRX`... Magnetic reconnection;
- `OTvortex` ... Orszag-Tang vortex problem;
- `RTI` ... Rayleigh-Taylor instability;
- `blast` ... blast wave propagation problem;
- `loop` ... field loop advection problem[^1].

### Project structure

- `common/` contains the shared MPI `MHD3D`, dissipative `DMHD3D`, and gravity-aware `GMHD3D` classes and solvers, plus the build rules in `common/case.mk`.
- Each problem directory retains its driver (`main.cpp`), initial conditions, parameters, macros, and a small `Makefile`.
- RTI uses the shared `GMHD3D`; its initial state, potential profile, and special upper/lower energy boundary conditions are in `RTI/gmhd3d_init_.cpp`.
- The existing `merge.out`, `batch.py`, and `python/` links remain in each problem directory.

The sources in this directory's `common/` are compiled separately for each problem, using that problem's `mymacros.hpp`. They are not a precompiled class library. The repository-level `common/` and `mpi/` supply the numerical kernels and MPI routines in `libqasmhd.a` and `libmympi.a`.

MRX selects the dissipative sources with `DISSIPATION = 1`; RTI selects the gravity sources with `GRAVITY := 1`, as in `2D/MPI`. The case Makefiles include `common/case.mk` and list their local sources in `CASE_SRCS`.

`common/mhd3d_control.cpp` provides time-step initialization and backup validation for all cases, including RTI. `common/mhd3d_io.hpp` provides checked output helpers and MPI-wide error termination.

Solver workspaces are reused through `static thread_local std::vector` storage. The calling thread allocates shared field arrays and separate line buffers for each OpenMP worker before entering parallel regions. Storage is resized when the required size changes and retained until thread exit; same-thread reentrant solver calls are unsupported.

### Configuring a problem

- Edit `mhd3d_init_.cpp` for initial conditions, `dmhd3d_init_.cpp` for MRX (also defines `setdc()`), or `gmhd3d_init_.cpp` for RTI (also defines its special `bound()`).
- Edit `mhd3d_paras.cpp` for the domain and boundary conditions. `setup_grid(xmin, xmax, ymin, ymax, zmin, zmax)` takes **global** physical domain bounds excluding ghost cells; it sets local coordinates, mesh spacings, and the initial timestep using the existing MPI decomposition. Optional `xshift, yshift, zshift` arguments specify offsets in cell widths: all default to `0.5`; OTvortex uses `0.0, 0.0, 0.0`.
- The output directory defaults to `./dat/`. To change it, assign the `std::string` member `fildir` in `paras()` and create the directory before running.
- Edit `mymacros.hpp` for global mesh size, MPI process counts, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), `R_K` (1–3), and positive `MNP_X/Y/Z` are checked by `static_assert` in the shared class header. `CTW` controls multidimensional CT upwinding[^1]. Uneven MPI partitions are supported; the process counts need not divide the corresponding global mesh sizes exactly.

`GMHD3D` inherits the common execution loop and includes a prescribed, time-independent gravitational potential in the conserved total energy. Gravity-aware pressure conversion does not modify that stored energy, and momentum sources use density at the start of each RK stage. RTI's special energy boundary correction applies only to the actual conserved state at global Z boundaries; it is not a general prescription for gravity. MRX uses the diffusion substep safety coefficient of 9, the 3D counterpart of 6 in `2D/MPI`; its time-splitting order is unchanged.

KHI, MRX, and RTI use `rand_noise_mt()` with `std::mt19937`. Rank 0 broadcasts an unsigned time-based base seed when random perturbations are used. For reproducible runs, replace `std::time(nullptr)` in the case's initial-condition source with a fixed value such as `10u`, then rebuild. Each generator is initialized once per `init_()` call and reused while filling the perturbation arrays.

| Case and random amplitude | Seed inputs to `std::seed_seq` | Ranks sharing the same sequence |
| --- | --- | --- |
| KHI/MRX: `vy` perturbation | `seed, mpi_ranx, mpi_ranz` | Same X/Z subdomain, independent of Y rank |
| KHI: shear amplitude in `vx` | `seed, mpi_ranz` | Same Z subdomain, independent of X/Y rank |
| MRX: magnetic perturbation amplitude `b1amp` | `seed, mpi_ranz` | Same Z subdomain, independent of X/Y rank |
| RTI: `vz` perturbation | `seed, mpi_ranx, mpi_rany` | Same X/Y subdomain, independent of Z rank |

The table describes the random amplitudes; each case then applies its spatial profiles. `RANDOM=0` selects a sinusoidal `vy` perturbation in KHI, no `vy` perturbation in MRX, and a cosine-product `vz` perturbation in RTI. **KHI's `vx` and MRX's magnetic perturbation remain random even with `RANDOM=0`**, using separate Z-only generators unaffected by the `vy` generator.

Random sequences differ from the legacy generator. Reproducibility requires the same base seed, mesh, decomposition, and execution/call history; changing the MPI layout does not preserve the same global random field. Initial perturbation arrays use `std::vector`; allocation/initialization exceptions and failed seed broadcasts stop all ranks. These three cases use random numbers only for initialization, not during time integration.

### How to run the simulation

From `3D/MPI/`:

```sh
cd OTvortex/
make -C ../../../common
make -C ../../../mpi
make
mkdir -p dat
mpiexec -np 16 env OMP_NUM_THREADS=2 ./a.out
```

Use the MPI C++ compiler `MPICXX` and flags configured in the repository's `Makefile.inc`. The MPI process count must equal `MNP_X*MNP_Y*MNP_Z` in `mymacros.hpp`. The example uses sixteen MPI processes, each with two OpenMP threads; `OMP_NUM_THREADS=1` uses one thread per process.

The executable is `a.out`; object and dependency files are stored in each problem's `build/`. Results are written to `dat/`. Local macro and shared-header edits are tracked by dependency files: rerun `make` after editing. Rebuild the repository-level libraries after changing their sources; after changing compiler flags, clean and rebuild both the library objects and the case objects.

Rerunning the same executable and MPI layout loads a saved backup from the output directory when available. Preserve previous results and use a separate directory for a fresh run; see the restart requirements below.

`make clean` removes `a.out` and the current build's object/dependency files, but preserves results. It does not remove objects left directly in a case directory by the previous layout; those files are no longer linked. `make cdata` deletes `dat/*.dat`, including backups; use it deliberately.

### Time stepping and restart

Time integration follows `2D/MPI`'s bookkeeping. On a new run, the initial CFL timestep is aligned to the output interval with `nrec=max(1,ceil(dtrec/dt))` and `dt=dtrec/nrec`.

- `exec_(0)` keeps this fixed timestep, runs `nrec*N_OUT` steps, and records every `nrec` steps.
- `exec_(1)` retains the existing CFL recalculation every two steps, runs until `tmax=N_OUT*DTREC`, and shortens the final step to reach `tmax` exactly. The initial step-count estimate does not limit the run. Records are written when the time crosses the next output threshold; intermediate record times need not equal that threshold exactly.
- Nonfinite/nonpositive timesteps and unsupported step counts stop all ranks. Common initialization/restart checks are in `common/mhd3d_control.cpp`.
- MRX uses the same actual step length for its ideal and dissipative updates.

Backup load/save results are collected with `mpi_sync_status()`: all ranks must either find no checkpoint or load one successfully. Missing rank files, invalid metadata, inconsistent rank statuses/times, and I/O errors stop the entire MPI job. A restored timestep is preserved, not realigned; fixed-step checkpoints must already be consistent with the current output interval and counters. Use the same mesh, MPI decomposition, physical parameters, and time-stepping mode/settings for restart. The old format does not encode all of these settings, so validation cannot detect every incompatible configuration.

Output checks use the same path/open/write/close conventions as SERIAL, but invoke `MPI_Abort()` on failure. The MPI file formats are unchanged: normal fields and gravitational potential use single precision, while backup fields use double precision. A new checkpoint is saved only after every rank has completed its normal output.

Backups still overwrite the existing files without retaining generations. This is **not** an atomic multi-file checkpoint: forced termination during a save can leave mixed or incomplete data, and size checks cannot detect every such case. Preserve a known-good checkpoint externally when required.

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

Trailing slashes are optional; quote directory names containing spaces. Uneven MPI partitions are supported. The utility retains the coordinate formats and the eight-field float32 layout of `merge_outdat_XXXXX.dat`, but fields 4, 5, and 6 (`Bx`, `By`, `Bz`, zero-based) now contain **cell-center** values. Other fields, including total energy, are unchanged. It copies `t.dat` and `params.dat` when necessary. Raw rank files are not modified. Stop the simulation before merging; input files must not change during the merge.

When all rank-local `g_potential_XXXXX.dat` files exist, they are also merged into `merge_g_potential.dat` (float32, global `(nz, ny, nx)`, without ghost cells). `batch.py` reads this optional file and subtracts `rho*phi_g` when calculating pressure, without modifying stored total energy in `data[7]`. Missing potential on only some ranks, invalid/nonfinite potential, and incorrect binary sizes are errors. With no rank potentials, the utility reports that gravity correction is disabled; it cannot distinguish a nongravitating run from a gravitational run whose potential files were all lost. It refuses to proceed if an old merged potential exists in the destination but no input potentials exist; use a clean output directory instead of reusing unrelated results.

Malformed metadata, missing rank data, allocation failures and I/O failures return a nonzero exit status. Each output is written to a temporary file and renamed only after successful close. This is not a transaction across all output files or a power-loss guarantee: a failed/interrupted merge can leave a mixture of old and new complete files. Rerun successfully before plotting, preferably into a separate directory.

The simulation outputs only CT face values for the magnetic fields in its rank-local files. The merger computes cell-center fields and, for each output index, also writes the original CT face fields:

| File | float32 shape | Location |
| --- | --- | --- |
| `merge_outdat_XXXXX.dat` | `(8, nz, ny, nx)` | All fields at cell centers |
| `merge_bx_face_XXXXX.dat` | `(nz, ny, nx+1)` | X faces, including the global upper X face |
| `merge_by_face_XXXXX.dat` | `(nz, ny+1, nx)` | Y faces, including the global upper Y face |
| `merge_bz_face_XXXXX.dat` | `(nz+1, ny, nx)` | Z faces, including the global upper Z face |

Center fields are the averages of the two neighboring faces in their respective directions. Ghost cells in the rank-local data supply the final upper faces; at least one ghost cell per direction is required. Shared interfaces are taken from the rank on the positive side. Face values are copied without interpolation; no periodic wrapping is imposed at physical boundaries. The CT companion files add approximately 37.5% to the eight-field merged storage for large grids. Two-point interpolation is an output convention and need not reproduce the solver's higher-order interpolation.

In Python, **`bx, by, bz` are cell-center fields** loaded directly from `merge_outdat`. `bxct, byct, bzct` are the separately loaded staggered fields. `divb` is the second-order CT divergence at all cell centers, with shape `(nz,ny,nx)`. Currents are computed only at interior CT edges; outer edges require tangential ghost data that are not included:

| Variable | Shape | Coordinate arrays (X, Y, Z) |
| --- | --- | --- |
| `jx` | `(nz-1, ny-1, nx)` | `x, yf, zf` |
| `jy` | `(nz-1, ny, nx-1)` | `xf, y, zf` |
| `jz` | `(nz, ny-1, nx-1)` | `xf, yf, z` |

Here `xf, yf, zf` are the midpoints between adjacent cell coordinates. These diagnostics use float32 output and second-order differences, not necessarily the solver's full-precision/high-order operators.

Rebuild `merge.out` and rerun it on raw data before using the updated Python scripts. Old merged fields stored face values under the same eight-field layout and are **not interchangeable** with the new center fields. Missing or wrong-size CT companion files cause an error; retain all four files for each time index.

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
>>>plt3d.slice(x,y,z,ro,plane_orientation='y_axes',figure=fig) #Overplot in x-z plane
>>>plt3d.slice(x,y,z,ro,plane_orientation='z_axes',figure=fig) #Overplot in x-y plane
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../license/COPYING) file for details.

[^1]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
