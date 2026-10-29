## Blast Wave Problem

This problem suggested by Balsara & Spicer[^1] demonstrates the propagation of MHD shocks in strongly magnetized meidum to assess the robustness of the code.<br>
The initial condition is *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>,P)*=*(1,0,0,0,10cos(&theta;),10sin(&theta;),0,1)*, and then a high pressure cylinder is imposed at the center of the domain, *P=100* for *&radic;(x<sup>2</sup>+y<sup>2</sup>)&le;0.125* where *-2<x,y<2*. The angle is measured from the x axis, as in `mhd2d_init_.cpp`.

The example result of the gas and magnetic pressures for *&theta;=30&deg;* is shown below (in logarithmic scale).

<img src="../../imgs/blast/Figure_1.png" alt="Gas pressure in blast wave" width="500px"> <img src="../../imgs/blast/Figure_2.png" alt="Magnetic pressure in blast wave" width="500px">

### Serial configuration

This case uses the shared `MHD2D` class and solver in `../common/`. Edit `mhd2d_init_.cpp` for the cylinder radius, densities, pressures, and magnetic-field angle. `mhd2d_paras.cpp` sets the domain with `setup_grid()` and defines the boundary conditions; `mymacros.hpp` controls the mesh, output, and solver choices.

See the [serial build and visualization instructions](../README.md). Build with `make` in this directory; objects are kept in `build/`, the executable is `a.out`, and results are written to `dat/`. Run the linked `python batch.py` from this directory to inspect the results.

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.

[^1]: [Balsara, D. S., and Spicer, D. S. 1999, JCP](https://www.sciencedirect.com/science/article/abs/pii/S0021999198961538?via%3Dihub)
