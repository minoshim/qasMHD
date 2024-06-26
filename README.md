# qasMHD
Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations with a finite-difference method.<br>
The qasMHD code has the characteristics of:
- up to 4th order accuracy in space and 3rd order accuracy in time;
- shock capturing by approximate Riemann solvers, including state-of-the-art low-dissipation all-speed solvers[^1][^2][^3];
- preservation of the solenoidal condition of the magnetic field by a well-designed Contrained Transport method[^4].

The qasMHD code is unique in that it can accurately solve MHD flows in wide-ranging Mach numbers, even when the flow is almost incompressible!<br>
The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

## System requirements
Following packages are required to be installed on the system:
- Git to install and update the qasMHD code;
- C++ compiler (GNU, Intel);
- MPI library (MPICH, OpenMPI) to use the MPI parallel code;
- Make
- Python 3.X with NumPy and matplotlib for data read and visualization (can be downloaded from [Anaconda](https://www.anaconda.com/products/distribution));
- (Optional) [Mayavi](https://mayavi.readthedocs.io/ja/latest/index.html) in Python for 3D data visualization (can be installed by `>pip install mayavi`).

The code is tested on Linux OSs (Ubuntu, Linux Mint, CentOS, including Windows Subsystem for Linux).<br>
Above packages (excluding Intel compiler) can be installed by `apt` or `yum` commands.

## Installation
1. Download the qasMHD code from GitHub by `>git clone https://github.com/minoshim/qasMHD`.
2. Move to the main directory `qasMHD/`.
3. Check `Makefile.inc` and edit environment variables `CC`, `CFLAGS`, and `MPICC` to meet users environment.
4. Execute `>make clean` and `>make` commands to remake libraries.

Since the code is updated without notice, users may need to update the code by `>git pull origin main`.

## Composition
- `1D/` contains serial codes for one-dimensional problems.
- `2D/` contains serial and MPI parallel codes for two-dimensional problems.
- `3D/` contains MPI parallel codes for three-dimensional problems.
- `common/` contains central functions for the simulation.
- `license/` contains license documents.
- `mpi/` contains functions for MPI parallelization.
- `Makefile` to make libraries.
- `Makefile.inc` to define environment variables.
- `README.md` is this file.
- `libmympi.a` is the library generated from `mpi/`.
- `libqasmhd.a` is the library generated from `common/`.

`common/mhd_class.hpp` and `mpi/mympi_class.hpp` define the base class of the MHD simulation and MPI parallelization.<br>
Member variables defined in `common/mhd_class.hpp` are:
```
  double gam=5.0/3.0;           // Specific heat ratio
  double *x,*y,*z;              // Spatial coordinate
  double *ro,*mx,*my,*mz,*en;   // Density, momentum, and total energy
  double *bx,*by,*bz;           // Magnetic field @ cell edge
  double *vx,*vy,*vz,*pr;       // Velocity and pressure
  double *cx,*cy,*cz;           // Magnetic field @ cell center
  double *nu,*eta;              // Kinematic viscosity and resistivity
  double *phi_g;                // Gravitational potential
```
and defined in `mpi/mympi_class.hpp` are:
```
  int mpi_num=1;                // Number of MPI processes
  int mpi_rank=0;               // MPI rank
```

For information about `1D/`, `2D/`, and `3D/` problems, see `README` in each directory.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
[^3]: [Minoshima, T. 2023, JPCS](https://iopscience.iop.org/article/10.1088/1742-6596/2742/1/012026)
[^4]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
