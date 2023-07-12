#include "mhd1d_class.hpp"

int main(){

  MHD1D mhd1d;
  mhd1d.paras();		// Set parameters
  mhd1d.init_();		// Set initial condition
  mhd1d.setdt(1);		// Set time step to satisfy CFL (if flag=1). 
  mhd1d.exec_(0);		// Run. flag 1/O with/without modifies dt @ every step.

  return 0;
}
