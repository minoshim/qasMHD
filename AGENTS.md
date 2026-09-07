# Repository Guidelines

## Project Structure & Module Organization

`common/` contains the shared MHD numerics and builds `libqasmhd.a`; `mpi/` contains halo exchange and boundary routines and builds `libmympi.a`. Runnable cases live under `1D/SERIAL/`, `2D/SERIAL/`, `2D/MPI/`, and `3D/MPI/`. Each case keeps its initial conditions in `mhd[123]d_init_.cpp`, domain and boundary settings in `mhd[123]d_paras.cpp`, compile-time solver choices in `mymacros.hpp`, and its update loop in `mhd[123]d_solve.cpp`. Python post-processing scripts are in the dimension-specific `python/` directories; example figures are in `1D/imgs/`, `2D/imgs/`, and `3D/imgs/`. `references/` stores the papers behind the numerical methods. Runtime output belongs in each case's `dat/` directory.

## Build, Test, and Development Commands

- `make clean && make`: rebuild the two static libraries using settings in `Makefile.inc`.
- `cd 1D/SERIAL/shock && make && ./a.out`: build and run a serial example.
- `cd 2D/MPI/OTvortex && make`: build an MPI/OpenMP example.
- `mpiexec -np 4 env OMP_NUM_THREADS=2 ./a.out`: run that example; `-np` must equal `MNP_X*MNP_Y` (and `*MNP_Z` in 3D).
- `./merge.out dat/ dat/`: merge rank-local MPI output before using `python/batch.py`.
- `make clean`: remove local objects and executables. `make cdata` deletes case output, so use it deliberately.

## Coding Style & Naming Conventions

Follow the existing GNU C++ style: two-space indentation, braces on a new line for functions and on the same line for control statements, and short `snake_case` routine names. Classes use uppercase domain names such as `MHD2D` and `DMHD3D`. Keep solver macros uppercase (`RMN`, `ODR`, `R_K`, `CTW`). Preserve the flattened indexing convention (`ss=nx*j+i` in 2D). There is no configured formatter or linter; keep diffs consistent with neighboring code.

## Testing Guidelines

No automated test suite or coverage target exists. Validate shared-kernel changes with at least one smooth-wave case and one shock case. For multidimensional changes, also check a CT case such as `OTvortex` or `loop`, inspect pressure/density for nonphysical values, and verify that discrete `divB` remains near roundoff. Record compiler, grid, solver macros, MPI layout, and comparison results in the PR.

## Commit & Pull Request Guidelines

History favors short imperative summaries such as `add van Leer limiter` and `update readme`. Use a concise subject describing one logical change. PRs should explain the numerical motivation, list affected cases, document validation commands and results, and include plots when solution behavior changes. Do not commit generated `dat/*.dat`, executables, object files, archives, or Python cache files.

## Done when (added by TM)
- The code builds successfully, or the remaining failure is clearly explained.
- The diff is small and understandable.
- Any numerical, performance, or compatibility risks are explicitly noted.
