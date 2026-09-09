## Rayleigh-Taylor instability

When a heavy fluid is located above a light fluid, the system is unstable against the gravity, known as the Rayleigh-Taylor instability.<br>
The instability also occurs when a heavy fluid is accelerated toward a light fluid.<br>
The interface is set at *y=y<sub>0</sub>=&pm;56*, and the density is *&rho;=1.0* (for *y>|y<sub>0</sub>|*) and *&rho;=0.1* (for *y<|y<sub>0</sub>|*) where *-80<y<80* (the system is symmetric with respect to *y=0*).<br>
The gravitational potential is given as *&Phi;=g<sub>0</sub> &lambda;<sub>g</sub>log(cosh(y/&lambda;<sub>g</sub>))* where *g<sub>0</sub>=1.0* and *&lambda;<sub>g</sub>=8&Delta;y*.<br>
The pressure is determined to satisfy the hydrostatic equilibrium.<br>

Example density profiles at t=15.0,20.0 are shown below.

![RTI1](../../imgs/RTI/Figure_1.png)
![RTI2](../../imgs/RTI/Figure_2.png)

### Serial configuration

This case uses the shared gravity-aware `GMHD2D` implementation in `../common/`; its Makefile selects the gravity sources with `GRAVITY := 1`. `gmhd2d_init_.cpp` defines the initial state, prescribed potential, and `GMHD2D::bound()`. The latter contains the special upper/lower energy correction for this RTI setup, not a general-purpose gravitational boundary condition.

Edit `mhd2d_paras.cpp` for `setup_grid()` and ordinary boundary flags. `RANDOM=0` in `mymacros.hpp` uses the single-mode velocity perturbation; `RANDOM=1` uses `rand_noise_mt()`. Set a fixed `seed` in `gmhd2d_init_.cpp` for reproducibility; the default is time-based.

The conserved total energy includes `rho*phi_g`. The potential is written once to `dat/g_potential.dat` as binary doubles with shape `(ny, nx)`, including ghost cells. The shared `batch.py` and `batch_a.py` load this file when present and subtract `rho*phi_g` before deriving gas pressure and temperature. `batch_a.py` applies the same time-independent potential to every output, using each output's density. Missing files select the ordinary MHD formula; invalid file sizes or non-finite potential values raise an error. Keep `g_potential.dat` with the corresponding RTI output when copying or archiving results.

See the [serial build and visualization instructions](../README.md). Build with `make` in this directory; objects are kept in `build/`, the executable is `a.out`, and results are written to `dat/`.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.
