## Wave Propagation Problem
This problem demonstrates the propagation of Hall-MHD waves in homogeneous plasma to confirm that the code successfully satisfies the linear dispersion.<br>
Random white noise is added to the magnetic field and pressure, and their perturbations propagate as the whistler and ion cyclotron,  and magnetosonic waves.

The result for *&beta;=10* and the ion inertia length *d<sub>i</sub>=4&Delta;x* is shown below.<br>
Dispersion relation of the whistler and ion cyclotron waves is observed in the spectrum plot of the magnetic field.
![Wave profile for beta=10](../../imgs/h-wave/wave_b1e+1.png)
![Spectrum of by for beta=10](../../imgs/h-wave/wk_by_b1e+1.png)

It is well known that the dispersive nature of the whistler waves (phase velocity proportional to wavenumber) is unfavorable condition for Hall-MHD simulations in terms of numerical stability.<br>
To relax this problem, the code can introduce the upper bound to the phase velocity of whistler waves through an "artificial" electron inertia, motivated by the study of Amano[^1].<br>

The result for *&beta;=10*, the ion inertia length *d<sub>i</sub>=4&Delta;x*, and the electron inertia length *d<sub>e</sub>=&Delta;x* is shown below.<br>
The dispersion relation of the whistler wave *&omega; = V<sub>a</sub>d<sub>i</sub>k<sup>2</sup>* is modified to *&omega; = V<sub>a</sub>d<sub>i</sub>k<sup>2</sup>/(1+d<sub>e</sub><sup>2</sup>k<sup>2</sup>)*, limiting the maximum phase velocity to *V<sub>a</sub>d<sub>i</sub>/d<sub>e</sub>*.
![Wave profile for beta=10](../../imgs/h-wave/wave_b1e+1_e.png)
![Spectrum of by for beta=10](../../imgs/h-wave/wk_by_b1e+1_e.png)

[^1]: [Amano, T. 2015, JCP](https://www.sciencedirect.com/science/article/abs/pii/S0021999115004805?via%3Dihub)
