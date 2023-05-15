#include "dmhd2d_class.hpp"

int main(){

  DMHD2D mhd2d;
  mhd2d.paras();		// Set parameters
  mhd2d.init_();		// Set initial condition
  mhd2d.setdt(1);		// Set time step to satisfy CFL (if flag=1). 
  mhd2d.exec_(1);		// Run. flag 1/O with/without modifies dt @ every step.

  return 0;
}
