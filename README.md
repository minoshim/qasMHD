# qasMHD
Quasi all-speed magnetohydrodynamic (qasMHD) simulation code is a C++ package for solving compressible MHD equations.
The qasMHD code is unique in that it can accurately solve very wide-ranging Mach number MHD flow (even neary incompressible), by employing
- Low-dissipation version of the HLLD approximate Riemann solver
- Central Upwind Constrained Transport method

The current version supports 1D-3D space in Cartesian coordinates and hybrid MPI/OpenMP parallel computation.

