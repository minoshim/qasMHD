## Rayleigh-Taylor instability

When a heavy fluid is located above a light fluid, the system is unstable against the gravity, known as the Rayleigh-Taylor instability.<br>
The instability also occurs when a heavy fluid is accelerated toward a light fluid.<br>
The interface is set at *y=y<sub>0</sub>=&pm;56*, and the density is *&rho;=1.0* (for *y>|y<sub>0</sub>|*) and *&rho;=0.1* (for *y<|y<sub>0</sub>|*) where *-80<y<80* (the system is symmetric with respect to *y=0*).<br>
The gravitational potential is given as *&Phi;=g<sub>0</sub> &lambda;<sub>g</sub>log(cosh(y/&lambda;<sub>g</sub>))* where *g<sub>0</sub>=1.0* and *&lambda;<sub>g</sub>=8&Delta;y*.<br>
The pressure is determined to satisfy the hydrostatic equilibrium.<br>

Example density profiles at t=15.0,20.0 are shown below.

![RTI1](../../imgs/RTI/Figure_1.png)
![RTI2](../../imgs/RTI/Figure_2.png)

### Code organization

This case uses the shared gravity-aware `GMHD2D` implementation in `../common/`; its Makefile selects the gravity sources with `GRAVITY := 1`, as in SERIAL. `gmhd2d_init_.cpp` defines the initial state, prescribed potential, and `GMHD2D::bound()`. The latter contains the special upper/lower energy correction for this RTI setup, not a general-purpose gravitational boundary condition. Only ranks at global Y boundaries apply the correction, and only to the actual conserved-state arrays.

Edit `mhd2d_paras.cpp` for `setup_grid()` and ordinary boundary flags. The common MPI `MHD2D::exec_()` calls the overridden timestep, solver and output routines. Primitive conversion subtracts `rho*phi_g` locally without modifying stored total energy; gravity sources use the density at the start of each RK stage, matching SERIAL.

Output and backup formats are unchanged. The conserved total energy includes `rho*phi_g`; the static potential is written once to `dat/g_potential_XXXXX.dat` as binary floats with local shape `(ny, nx)`, including ghost cells, where `XXXXX` is the rank number. Keep these files with their corresponding fields and subtract `rho*phi_g` when deriving gas pressure. Restart regenerates the prescribed potential through `init_()` before loading conserved fields; use the same mesh, domain and physical parameters. Backups do not store the potential. Numerical results can differ from the old RTI solver because of the gravity-source correction and removal of energy round-trip rounding.

Run `make` in this directory; it compiles shared sources with the local `mymacros.hpp` and stores objects in `build/`. Build the merge utility with `make -C ../common`, then run `./merge.out dat` (equivalent to `./merge.out dat dat`), or specify a separate output directory as the second argument. This also creates `merge_g_potential.dat`; the common Python scripts read it and subtract `rho*phi_g` in pressure calculations. Keep the merged potential with the corresponding merged fields. If all potential files are absent, the utility cannot identify the data as gravitational, so check its status message before plotting RTI pressure. See [the MPI guide](../README.md) for build, execution, and result-merging instructions.

Random velocity perturbations use `rand_noise_mt()` with a common base seed and the X rank coordinate; ranks sharing an X subdomain use the same sequence across Y. The temporary `dvy` array is a checked `std::vector`. See [the MPI guide](../README.md) for fixed-seed runs and reproducibility limits.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.
