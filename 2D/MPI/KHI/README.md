## Kelvin-Helmholtz instability

The Kelvin-Helmholtz instability (KHI) is the hydrodynamic instability that occurs at the velocity shear.<br>
The initial condition is *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>,P)*=*(1,(1/2)tanh(y-10),0,0,cos(&theta;),0,sin(&theta;),P<sub>0</sub>)*, where *0<x,y<20*.<br>
The flow can amplify the in-plane magnetic field when *&theta; &ne; 90&deg;* via stretching, until the field becomes strong enough to stabilize the flow by the magnetic tension.

The instability is essentially incompressible phenomenon, and conventional shock-capturing schemes fails to resolve the low-Mach KHI due to excess numerical diffusion.
This issue is fixed by the all-speed scheme adopted here.

The following movies show the streamline for two example runs with *&theta; = 90&deg;* and *&theta; = 71.5651&deg;*.<br>
The ambient pressure *P<sub>0</sub>=500* so that the Mach number of the initial flow is 0.0158.<br>
Images are produced with the [Line Integral Convolution](https://lic.readthedocs.io/en/latest/).

<img src="../../imgs/KHI/khi_movie_perp.gif" alt="KH instability (out-of-plane B)" width="480px"> <img src="../../imgs/KHI/khi_movie_oblique.gif" alt="KH instability (in-plane B)" width="480px">
