## Wave Propagation Problem
This problem demonstrates the propagation of MHD waves in homogeneous plasma to confirm that the code successfully satisfies the linear dispersion.<br>
Random white noise is added to the magnetic field and pressure, and their perturabations propagate as the Alfven and magnetosonic waves.

Example results for high and low beta cases are shown below:<br>
(Top...spatial profiles, bottom left...spectrum of the magnetic field, bottom right...spectrum of the pressure. Users can create the spectrum plot with the scripts `batch_a.py` and `wk_plt.py`)

These figures were regenerated with `XMESH=256`, `RMN=2` (LHLLD), `ODR=4`, `R_K=3`, `CFL=0.4`, `DTREC=0.01`, and `N_OUT=4095`, using GCC 16.2.0 (C++17, `-O2`) and `OMP_NUM_THREADS=2`. Initialization uses `rand_noise_mt()` with `std::mt19937` and seed `10`; the random realization differs from the legacy generator used for the old figures.

To reproduce them, set `beta` to `10.0` or `0.1` in `mhd1d_init_.cpp`, then rebuild and run as described in the [SERIAL guide](../README.md). In an interactive Python session, execute:

```pycon
>>> exec(open("batch_a.py").read())
>>> exec(open("wk_plt.py").read())
```

For `batch_a.py`, enter `dat` and select output index `1024` for the spatial profiles at `t=10.240`. The spectrum uses all 4096 snapshots from `t=0` to `40.95`. In `wk_plt.py`, use `val=by` for the magnetic-field spectrum; replace that assignment with `val=pr` and rerun the script for the pressure spectrum. The color scale is `log10(abs(FFT))`.

1. `beta=10` in `mhd1d_init_.cpp`. The pressure perturbation propagates as the fast magnetosonic wave.
![Wave profile for beta=10](../../imgs/wave/wave_b1e+1.png)
<img src="../../imgs/wave/wk_by_b1e+1.png" alt="Spectrum of by for beta=10" width="400px"><img src="../../imgs/wave/wk_pr_b1e+1.png" alt="Spectrum of pr for beta=10" width="400px">

2. `beta=0.1` in `mhd1d_init_.cpp`. The pressure perturbation propagates as the slow magnetosonic wave.
![Wave profile for beta=0.1](../../imgs/wave/wave_b1e-1.png)
<img src="../../imgs/wave/wk_by_b1e-1.png" alt="Spectrum of by for beta=0.1" width="400px"><img src="../../imgs/wave/wk_pr_b1e-1.png" alt="Spectrum of pr for beta=0.1" width="400px">

## License

This project is licensed under the GNU General Public License v3.0 - see the [license](../../../license/COPYING) file for details.
