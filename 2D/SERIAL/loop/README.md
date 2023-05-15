## Field Loop Advection Problem

This problem is first suggested by Gardiner & Stone[^1] to assess the capability to preserve the shape of a weak magnetic field loop during the advection.<br>
The initial condition is *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>z</sub>,P)*=*(1,&radic;5cos(&theta;),&radic;5sin(&theta;),0,0,P<sub>0</sub>)*, and *(B<sub>x</sub>, B<sub>y</sub>)* are calculated from the vector potential *A<sub>z</sub>=*max*(10<sup>-3</sup>(0.3-r),0)* where *r<sup>2</sup>=x<sup>2</sup>+y<sup>2</sup>* and *-1<x<1, -0.5<y<0.5*.<br>
The problem is stringent especially when the flow is well aligned to the axis *(&theta;~0 or &pi;/2)*[^2] and/or the flow is subsonic *(P<sub>0</sub>>>1)*[^3].

The example results of the magnetic pressure for *&theta;=0.01* and *P<sub>0</sub>=1* are shown below.

<img src="../../imgs/loop/Figure_1.png" alt="Loop advection t=1" width="500px"> <img src="../../imgs/loop/Figure_2.png" alt="Loop advection t=2" width="500px">

[^1]: [Gardiner T. A., and Stone J. M. 2005, JCP](https://www.sciencedirect.com/science/article/pii/S0021999104004784)
[^2]: [Lee, D. 2013, JCP](https://www.sciencedirect.com/science/article/pii/S0021999113001836?via%3Dihub)
[^3]: [Minoshima T., Miyoshi T., and Matsumoto Y. 2019, ApJS](https://iopscience.iop.org/article/10.3847/1538-4365/ab1a36/meta)
