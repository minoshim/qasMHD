#include "mhd3d_class.hpp"

int main(int argc, char* argv[]){

  MHD3D mhd3d(&argc,&argv,MNP_X*MNP_Y*MNP_Z);
  mhd3d.paras();		// Set parameters
  mhd3d.init_();		// Set initial condition
  mhd3d.setdt(1);		// Set time step to satisfy CFL (if flag=1). 
  mhd3d.exec_(1);		// Run. flag 1/O with/without modifies dt @ every step.

  return 0;
}
