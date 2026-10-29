## Magnetic reconnection

In highly conducting plasma with antiparallel magnetic fields, dissipation of electric current causes rapid change of magnetic topology and the convertion of magnetic energy into kinetic and thermal energies, known as magnetic reconnection.<br>
The initial condition is the Harris equilibrium, *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>)*=*(&rho;<sub>0</sub>/cosh<sup>2</sup>(y/&lambda;)+&rho;<sub>1</sub>,0,0,0,B<sub>0</sub>tanh(y/&lambda;),0,0)* where *&rho;<sub>0</sub>=B<sub>0</sub>=&lambda;=1,&rho;<sub>1</sub>=0.2*, and *-4<y<4*.<br>
The pressure is determined from the pressure balance *(P+B<sup>2</sup>/2=const)* and the plasma beta at the upstream region *(|y|>>&lambda;)* is 0.2.<br>
The reconnection is triggered by adding a magnetic perturbation around the center.

In MHD simualtions, the evolution of magnetic reconnection strongly depends on the resistivity model, which is given phenomenologically.<br>
Uniform resistivity leads to elongated current sheet and slow reconnection (Sweet-Paker model), while localized resistivity causes fast Petscheck-type reconnection.<br>
Following example is the out-of-plane current with localized resistivity (*&propto;exp(-(x<sup>2</sup>+y<sup>2</sup>)/&lambda;)*, defined in `dmhd2d_init_.cpp`), showing Petscheck-like bifurcated structure.

![MRX](../../imgs/MRX/Figure_1.png)

### Code organization

The MPI `MHD2D` and `DMHD2D` implementations are in `../common/`. This directory retains `dmhd2d_init_.cpp` (including `setdc()`), parameters, macros, and the driver. The Makefile selects the dissipative sources with `DISSIPATION = 1`.

Run `make` in this directory; it compiles shared sources with the local `mymacros.hpp` and stores objects in `build/`. See [the MPI guide](../README.md) for build, execution, and result-merging instructions.

Random velocity perturbations use `rand_noise_mt()` with a common base seed and the X rank coordinate; ranks sharing an X subdomain use the same sequence across Y. The temporary `dvy` array is a checked `std::vector`. See [the MPI guide](../README.md) for fixed-seed runs and reproducibility limits.

The diffusion substep safety coefficient is 6, matching SERIAL. The existing ideal/dissipative splitting order is retained.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.
