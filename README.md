# qasMHD

Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations with a finite-difference method.<br>
The qasMHD code has the characteristics of:

- up to 4th order accuracy in space and 3rd order accuracy in time;
- shock capturing by approximate Riemann solvers, including state-of-the-art low-dissipation all-speed solvers[^1][^2][^3];
- preservation of the solenoidal condition of the magnetic field by a Constrained Transport (CT) method[^4].

The qasMHD code is unique in that it can accurately solve MHD flows in wide-ranging Mach numbers, even when the flow is almost incompressible!<br>
The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

## System requirements

The following packages are required:

- Git to install and update the qasMHD code;
- a C++17 compiler with OpenMP support (the default is `g++`);
- an MPI library (MPICH or Open MPI) and its C++ compiler wrapper (`mpic++` by default) for MPI cases;
- GNU Make;
- Python 3 with NumPy and Matplotlib for data reading and visualization;
- optionally, [Mayavi](https://mayavi.readthedocs.io/ja/latest/index.html) for the 3D visualization examples.

The code has been tested on Linux (Ubuntu, Linux Mint, CentOS, including Windows Subsystem for Linux). SERIAL cases use one process with OpenMP support and do not require MPI.

## Installation

Clone the repository and enter its root directory:

```sh
git clone https://github.com/minoshim/qasMHD.git
cd qasMHD
```

Check [Makefile.inc](Makefile.inc) and adjust `CXX` and `MPICXX` for your environment. `CXXSTD` selects C++17; `OPTFLAGS`, `OMPFLAGS`, `WARNFLAGS`, and `LTOFLAGS` control optimization, OpenMP, warnings, and optional link-time optimization. These settings are shared by the libraries, cases, and merge utilities. `CXXFLAGS` combines the compilation flags; `LDFLAGS` includes the OpenMP and LTO flags needed at link time.

Rebuild the libraries for your compiler:

```sh
make clean
make
```

This builds `libqasmhd.a` from `common/` and `libmympi.a` from `mpi/`. The root Makefile builds the libraries only; cases and merge utilities are built separately. For SERIAL-only use without MPI, replace the commands above with:

```sh
make -C common clean
make -C common
```

## Composition

| Problem directory and guide | Execution | Shared classes |
| --- | --- | --- |
| [1D/SERIAL](1D/SERIAL/README.md) | One process with OpenMP | `MHD1D`, Hall-MHD `HMHD1D` |
| [2D/SERIAL](2D/SERIAL/README.md) | One process with OpenMP | `MHD2D`, dissipative `DMHD2D`, gravity-aware `GMHD2D` |
| [2D/MPI](2D/MPI/README.md) | MPI/OpenMP | `MHD2D`, `DMHD2D`, `GMHD2D` |
| [3D/MPI](3D/MPI/README.md) | MPI/OpenMP | `MHD3D`, `DMHD3D`, `GMHD3D` |

Each directory above has its own `common/` for shared classes, solvers, and `case.mk` build rules, and a `python/` directory for post-processing. Each problem retains its driver (`main.cpp`), initial conditions, parameters, `mymacros.hpp`, and a small `Makefile`. Case-specific boundary implementations remain with the problem, such as the derived `RMI2D` class and RTI's gravity boundary conditions.

These shared class sources are compiled separately for each problem, using that problem's `mymacros.hpp`; they are not a precompiled class library. Objects and dependency files are stored in each problem's `build/`, the executable is `a.out`, and the default output directory is `dat/`.

At the repository root:

- `common/` contains the numerical kernels and builds `libqasmhd.a`;
- `mpi/` contains MPI communication and boundary routines and builds `libmympi.a`;
- `license/` contains license documents;
- `Makefile` builds the libraries, and `Makefile.inc` defines shared build settings.

### Base classes

[common/mhd_class.hpp](common/mhd_class.hpp) and [mpi/mympi_class.hpp](mpi/mympi_class.hpp) define the base classes for MHD simulations and MPI execution. Main protected members of `MHD` include:

```cpp
  double gam=5.0/3.0;           // Specific heat ratio
  double *x,*y,*z;              // Spatial coordinates
  double *ro,*mx,*my,*mz,*en;   // Density, momentum, and total energy
  double *bx,*by,*bz;           // Staggered magnetic field (CT face values)
  double *vx,*vy,*vz,*pr;       // Velocity and pressure
  double *cx,*cy,*cz;           // Magnetic field @ cell center
  double *nu,*eta;              // Kinematic viscosity and resistivity
  double *phi_g;                // Gravitational potential
```

The MPI base class provides:

```cpp
  int mpi_num=1;                // Number of MPI processes
  int mpi_rank=0;               // MPI rank
```

## Configuring a problem

- Edit the case's initial-condition and parameter sources listed in its guide. `setup_grid()` takes physical domain bounds excluding ghost cells; MPI cases use global bounds.
- Edit `mymacros.hpp` for mesh sizes, output intervals, CFL, and solver choices. `RMN` (0–3), `ODR` (1–4), and `R_K` (1–3) are checked in the shared class headers. Multidimensional cases also use `CTW` for CT upwinding.
- Set MPI process counts with `MNP_X`, `MNP_Y`, and, in 3D, `MNP_Z`. Uneven partitions are supported; process counts need not divide the global mesh sizes exactly.
- To change the output directory, assign `fildir` in `paras()` and create the directory before running. See the case guide for random-perturbation seeds and reproducibility settings.

## How to run the simulation

After building the libraries, run a SERIAL example from the repository root:

```sh
cd 1D/SERIAL/shock
make
mkdir -p dat
OMP_NUM_THREADS=2 ./a.out
```

For an MPI example, also start from the repository root:

```sh
cd 2D/MPI/OTvortex
make
mkdir -p dat
mpiexec -np 4 env OMP_NUM_THREADS=2 ./a.out
```

The MPI process count must equal `MNP_X*MNP_Y` in 2D or `MNP_X*MNP_Y*MNP_Z` in 3D. The 2D OTvortex defaults use four ranks; adjust `-np` when changing the decomposition. `OMP_NUM_THREADS` sets threads per process; use `1` for one thread per process.

Preserve previous results before a fresh run. SERIAL runs can overwrite output; MPI runs load a saved backup from the output directory when available. Restart requires matching mesh, MPI decomposition, physical parameters, and time-stepping settings; see the MPI guides for checkpoint limitations.

## How to check the result

For MPI output, first build the shared merge utility and combine rank-local files from the problem directory after the simulation has stopped:

```sh
make -C ../common
./merge.out dat
```

The single directory is used for both input and output, equivalent to `./merge.out dat dat`. To write merged output elsewhere:

```sh
mkdir -p merged
./merge.out dat merged
```

The merger uses `CXX` and runs without an MPI launcher. The case's `merge.out` link points to the shared executable. SERIAL output does not need merging.

For either SERIAL or MPI results, start Python 3 from the problem directory:

```sh
python
```

Then execute the linked script in the interactive session:

```pycon
>>> exec(open("batch.py").read())
```

Enter the data directory (`dat`, or `merged` for the separate MPI destination) and select an output index. The `batch.py` and `python/` links make the scripts available from each problem directory; the interactive session keeps loaded arrays available for inspection. See the directory guides for loading all output times, wave spectra, and Mayavi visualization.

MPI merged output stores cell-center magnetic fields in `merge_outdat_XXXXX.dat`, with separate CT face-field files for current and divergence diagnostics. Keep these companion files together. Rebuild the merger and regenerate old merged output from raw rank files before using the updated Python scripts; the old and new magnetic-field layouts are not interchangeable.

For gravity cases, the scripts subtract `rho*phi_g` when the corresponding potential file is present: `g_potential.dat` for SERIAL, or `merge_g_potential.dat` for MPI. Keep the matching potential with each data set. File formats and diagnostic array locations are described in the directory guides.

## Updating and rebuilding

Use `git pull --ff-only` to update the checked-out branch from its configured upstream. Rebuild the root libraries after changing their sources, then rebuild the cases and MPI merge utilities you use. Case dependency files track local macro and shared-header edits.

After changing compilers or flags in `Makefile.inc`, clean and rebuild the libraries, cases, and merge utilities separately; changing flags alone does not rebuild every existing object. Root `make clean` cleans only the libraries and their objects, not the case or merge builds.

In a problem directory, `make clean` removes `a.out` and the current build's object/dependency files while preserving results. `make cdata` deletes `dat/*.dat`, including MPI backups; use it deliberately.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima, T. 2023, JPCS](https://iopscience.iop.org/article/10.1088/1742-6596/2742/1/012026)
[^4]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
