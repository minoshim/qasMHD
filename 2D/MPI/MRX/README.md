## Magnetic reconnection

In highly conducting plasma with antiparallel magnetic fields, dissipation of electric current may cause rapid change of magnetic topology and the convertion of magnetic energy into kinetic and thermal energies, known as magnetic reconnection.
The initial condition is the Harris equilibrium (standard setup for reconnection study), *(&rho;,v<sub>x</sub>,v<sub>y</sub>,v<sub>z</sub>,B<sub>x</sub>,B<sub>y</sub>,B<sub>z</sub>)*=*(&rho;<sub>0</sub>/cosh<sup>2</sup>(y/&lambda;)+&rho;<sub>1</sub>,0,0,0,B<sub>0</sub>tanh(y/&lambda;),0,0)* where *&rho;<sub>0</sub>=B<sub>0</sub>=&lambda;=1,&rho;<sub>1</sub>=0.2*, and *-4<y<4*.<br>
The pressure is determined from the pressure balance *(P+B<sup>2</sup>/2=const)* and the plasma beta at the upstream region *(|y|>>&lambda;)* is 0.2.<br>
The reconnection is triggered by adding a magnetic perturbation around the center.

In MHD simualtions, the evolution of magnetic reconnection strongly depends on the resistivity model, which is described phenomenologically.<br>
Uniform resistivity leads to elongated current sheet and slow reconnection (Sweet-Paker model), while localized resistivity leads to fast Petscheck-type reconnection.<br>
Following example shows the out-of-plane current for localized resistivity (*&propto;exp(-(x<sup>2</sup>+y<sup>2</sup>)/&lambda;)*, defined in `dmhd2d_init_.cpp`).

![MRX](../../imgs/MRX/Figure_1.png)
