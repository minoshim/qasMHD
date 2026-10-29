## Richtmyer-Meshkov instability

When a shock collides with a corrugated contact discontinuity, the interface develops nonlinearly via the Richtmyer-Meshkov instability (RMI).<br>
The shock is initially located at *y=0* and the upstream variables at *y>0* are *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>,P)*=*(1,0,-1,0,0.0000346,0,0,0.00006)* so that the Mach number and the plasma beta are 100 and 10<sup>5</sup>[^1][^2].<br>
The downstream variables at *y<0* are set to satisfy the Rankine-Hugoniot condition.<br>
A corrugated contact discontinuity is imposed in the upstream region, and the density increases to *&rho;=10* behind the discontinuity.

Density profiles at *t=5.0,15.0* are shown below.<br>
This simulation was performed in a frame moving with *v<sub>y</sub>=-0.6* so that the structure of RMI remains at approximately *y=0*.

![RMI1](../../imgs/RMI/Figure_1.png)
![RMI3](../../imgs/RMI/Figure_3.png)

### Code organization

The base `MHD2D` class and ideal solver are in `../common/`. As in SERIAL, the local `RMI2D` class (`rmi2d_class.hpp/.cpp`) owns the inflow parameters, density cache, and random generator. Initial conditions are in `rmi2d_init_.cpp`; domain and ordinary boundary settings remain in `mhd2d_paras.cpp`.

Only ranks on the global upper Y boundary apply the inflow. `RMI2D::exec_()` refreshes its density once per time step; the overridden `bound()` retains that realization through all RK stages. Auxiliary magnetic boundary calls do not change the conserved state or consume random numbers. The physical upper `By` face remains CT-evolved.

Run `make` in this directory; it compiles shared sources with the local `mymacros.hpp` and stores objects in `build/`. See [the MPI guide](../README.md) for build, execution, and result-merging instructions.

Random perturbations use `rand_noise_mt()`: density streams are seeded per rank, while interface phases are shared through rank-0 broadcasts. Each instance owns its stream; no function-local initialization flag remains. With `RANDOM=1`, results change from the former implementation that regenerated inflow during RK boundary calls. Generator state is not checkpointed, so random-inflow restarts are not exactly reproducible. See [the MPI guide](../README.md) for fixed-seed runs and reproducibility limits.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.

[^1]: [Minoshima T., Kitamura K., and Miyoshi T. 2020, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab8aee/meta)
[^2]: [Minoshima T. and Miyoshi T. 2021, JCP](https://www.sciencedirect.com/science/article/pii/S0021999121005349)
