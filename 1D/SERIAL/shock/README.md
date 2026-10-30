## Shock Tube Problem
The shock tube problem is a standard benchmark test for compressible (M)HD simulations, in which fluids with two different thermodynamic states are separated by a diaphragm in a tube.<br>
The diaphragm is instantaneously removed at *t=0*, and then fluids evolve to produce a variety of shocks, rarefaction waves, and discontinuities depending on the initial condition.

Here, users can demonstrate five standard MHD shock tube problems adopted by Miyoshi and Kusano[^1], by editting the macro `NUM` in `mhd1d_init_.cpp`:
1. Dai & Woodward problem
2. Brio & Wu problem
3. Slow switch-off shock
4. Slow switch-off rarefaction
5. Super-fast expansion

Example results for each problem are shown below. These figures were regenerated with `XMESH=800`, `RMN=2` (LHLLD), `ODR=1`, `R_K=2`, `CFL=0.4`, and `DTREC=0.01`, using GCC 16.2.0 (C++17, `-O2`) and `OMP_NUM_THREADS=2`.

To reproduce them, set `NUM` to each problem number, rebuild and run as described in the [SERIAL guide](../README.md), then use `batch.py` to select output indices `20`, `10`, `20`, `20`, and `5`, respectively.

1. Dai & Woodward problem (`t=0.200`)
![Dai & Woodward problem](../../imgs/shock/Figure_1.png)
2. Brio & Wu problem (`t=0.100`)
![Brio & Wu problem](../../imgs/shock/Figure_2.png)
3. Slow switch-off shock (`t=0.200`)
![Slow switch-off shock](../../imgs/shock/Figure_3.png)
4. Slow switch-off rarefaction (`t=0.200`)
![Slow switch-off rarefaction](../../imgs/shock/Figure_4.png)
5. Super-fast expansion (`t=0.050`)
![Super-fast expansion](../../imgs/shock/Figure_5.png)

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.

[^1]:  [Miyoshi T. and Kusano K. 2005, JCP](https://www.sciencedirect.com/science/article/pii/S0021999105001142?via%3Dihub)
