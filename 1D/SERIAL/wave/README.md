## Wave Propagation Problem
This problem demonstrates the propagation of MHD waves in homogeneous plasma to confirm that the code successfully satisfies the linear dispersion.<br>
Random white noise is added to the magnetic field and pressure, and their perturabations propagate as the Alfven and magnetosonic waves.

Example results for high and low beta cases are shown below:<br>
(Top...spatial profiles, bottom left...spectrum of the magnetic field, bottom right...spectrum of the pressure. Users can create the spectrum plot with the scripts `batch_a.py` and `wk_plt.py`)

1. `beta=10` in `mhd1d_init_.cpp`. The pressure perturbation propagates as the fast magnetosonic wave.
![Wave profile for beta=10](../../imgs/wave/wave_b1e+1.png)
<img src="../../imgs/wave/wk_by_b1e+1.png" alt="Spectrum of by for beta=10" width="400px"><img src="../../imgs/wave/wk_pr_b1e+1.png" alt="Spectrum of pr for beta=10" width="400px">

2. `beta=0.1` in `mhd1d_init_.cpp`. The pressure perturbation propagates as the slow magnetosonic wave.
![Wave profile for beta=0.1](../../imgs/wave/wave_b1e-1.png)
<img src="../../imgs/wave/wk_by_b1e-1.png" alt="Spectrum of by for beta=0.1" width="400px"><img src="../../imgs/wave/wk_pr_b1e-1.png" alt="Spectrum of pr for beta=0.1" width="400px">

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.
