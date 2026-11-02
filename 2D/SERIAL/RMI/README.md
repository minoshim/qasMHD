## Richtmyer-Meshkov instability

When a shock collides with a corrugated contact discontinuity, the interface develops nonlinearly via the Richtmyer-Meshkov instability (RMI).<br>
The shock is initially located at *y=0* and the upstream variables at *y>0* are *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>,P)*=*(1,0,-1,0,0.0000346,0,0,0.00006)* so that the Mach number and the plasma beta are 100 and 10<sup>5</sup>[^1][^2].<br>
The downstream variables at *y<0* are set to satisfy the Rankine-Hugoniot condition.<br>
A corrugated contact discontinuity is imposed in the upstream region, and the density increases to *&rho;=10* behind the discontinuity.

Density profiles at *t=5.0,15.0* are shown below.<br>
This simulation was performed in a frame moving with *v<sub>y</sub>=-0.6* so that the structure of RMI remains at approximately *y=0*.

![RMI1](../../imgs/RMI/Figure_1.png)
![RMI3](../../imgs/RMI/Figure_3.png)

These figures were regenerated with `XMESH=128`, `YMESH=1024`, `RMN=2` (LHLLD), `ODR=2`, `R_K=3`, `CTW=1`, `CFL=0.4`, `DTREC=1`, `N_OUT=15`, `MAGNET=1`, and `RANDOM=0`, using GCC 16.2.0 (C++17, `-O2`) and `OMP_NUM_THREADS=2`. A fixed seed of `10` was used, although the perturbation is deterministic with `RANDOM=0`.

In `batch.py`, select output indices `5` or `15` and plot `ro` with `plt2d.image(x, y, ro, figsize=(4, 8), equal=1, xmin=0, xmax=1, ymin=-1, ymax=2, vmin=4, vmax=45)`. Output index `10` was also regenerated as [Figure_2.png](../../imgs/RMI/Figure_2.png). The titles use the actual times in `t.dat`, including the slight overshoot of the first requested time.

### Serial configuration

The local `RMI2D` class in `rmi2d_class.hpp/.cpp` derives from the shared `MHD2D` in `../common/`. Edit `rmi2d_init_.cpp` for the initial shock/contact state and inflow parameters, and `mhd2d_paras.cpp` for `setup_grid()` and ordinary boundary conditions.

`RANDOM=0` in `mymacros.hpp` selects the single-mode contact corrugation without density noise. `RANDOM=1` enables random phases and density perturbations. `rand_noise_mt()` uses one generator per instance, seeded in `rmi2d_init_.cpp`; replace the time-based `seed` with a fixed value for reproducibility.

`RMI2D::exec_()` refreshes the upper inflow once before each ideal-MHD timestep. The overridden `bound()` reuses the stored density values throughout the RK stages, without advancing the generator. The physical upper `By` face remains evolved by CT.

See the [serial build and visualization instructions](../README.md). Build with `make` in this directory; objects are kept in `build/`, the executable is `a.out`, and results are written to `dat/`.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
