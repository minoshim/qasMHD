# qasMHD
Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations.
The qasMHD code has the capabilities of
- up to 4th order accuracy in space, and 3rd order accuracy in time
- shock capturing by approximate Riemann solvers, including recent low-dissipation all-speed solvers
- preservation of the solenoidal condition of the magnetic field by the Contrained Transport (CT) method 
The qasMHD code is unique in that it can accurately solve MHD flows with very wide-ranging Mach numbers, even when the flow is neary incompressible!!
The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

## System requirements

## Instllation

