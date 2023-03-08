# qasMHD
Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations.

The qasMHD code has the capabilities of
- up to 4th order accuracy in space, and 3rd order accuracy in time,
- shock capturing by approximate Riemann solvers, including recent low-dissipation all-speed solvers[^1][^2],
- preservation of the solenoidal condition of the magnetic field by a sophisticated Contrained Transport (CT) method[^3].

The qasMHD code is unique in that it can accurately solve MHD flows in wide-ranging Mach numbers, even when the flow is neary incompressible!!

The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

## System requirements
Following packages are required to be installed on the system:
- Git for version management
- C++ compiler (GNU, Intel)
- MPI library (MPICH, OpenMPI) to use the MPI parallel code
- Python 3.X for data read and visualization
- These packages (excluding Intel compiler) can be installed by `apt-get` or `yum` commands

The qasMHD code is tested on Linux OSs (Ubuntu, Linux Mint, CentOS, including Windows Subsystem for Linux).

## Installation
1. Download the qasMHD code from Github via `>git clone hppts://`.
2. Move to the main directory `qasMHD/`.
3. Check `qasMHD/Makefile.inc` and modify environments variables `CC`, `CFLAGS`, and `MPICC` to meet your environments.
4. Execute `make clean` and `>make` commands to remake libraries
5. Since the code is updated without notice, you may update the code via `>git pull origin main`.

## Composition
- `1D/SERIAL/` contains serial codes for one-dimensional problems
- `2D/` contains serial and MPI parallel codes for two-dimensional problems
- `3D/MPI/` contains MPI parallel codes for three-dimensional problems
- `common/` contains central functions for MHD simulations
- `license/` contains license documents
- `mpi/` contains functions for MPI parallelization
- `Makefile` to make libraries there
- `Makefile.inc` to define environments variables
- `Readme.md` is this file
- `libmpi.a` is the library generated from `mpi/`
- `libqasmhd.a` is the library generated from `common/`

For information about `1D/`, `2D/`, and `3D/` problems, see `Readme.md` in each directory.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)