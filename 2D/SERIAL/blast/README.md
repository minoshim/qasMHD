## Blast Wave Problem

This problem demonstrates the propagation of MHD shocks in strongly magnetized meidum to assess the robustness of the code.

The initial condition is *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>,P)*=*(1,0,0,0,10sin(&theta;),10cos(&theta;),1)*, and then a high pressure cylinder is imposed at the center of the domain, *P=100* for *&radic;(x<sup>2</sup>+y<sup>2</sup>)<0.125* where *-2<x,y<2*.

The example result of the gas and magnetic pressures for *&theta;=30&deg;* is shown below (in logarithmic scale).

<img src="../imgs/blast/Figure_1.png" alt="Gas pressure in blast wave" width="500px"> <img src="../imgs/blast/Figure_2.png" alt="Magnetic pressure in blast wave" width="500px">
