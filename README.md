# qasMHD
Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations.
The qasMHD code has the capabilities of
- up to 4th order accuracy in space, and 3rd order accuracy in time
- shock capturing by approximate Riemann solvers, including recent low-dissipation all-speed solvers
- preservation of the solenoidal condition of the magnetic field by the Contrained Transport (CT) method 
The qasMHD code is unique in that it can accurately solve MHD flows with very wide-ranging Mach numbers, even when the flow is neary incompressible!!
The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

## System requirements
Following packages are required to be installed on the system
- Git for version management
- C++ compiler (GNU, Intel)
- MPI library (MPICH, OpenMPI) to use the MPI parallel code
- Python 3.X for data read and visualization
- These packages (excluding Intel compiler) can be installed by 'apt-get' or 'yum install' commands
The packages are tested on Linux OSs (Ubuntu, Linut Mint, CentOS, including Windows Subsystem for Linux)

## Instllation

