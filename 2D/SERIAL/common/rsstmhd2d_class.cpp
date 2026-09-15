#include "rsstmhd2d_class.hpp"
#include "mhd2d_io.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace serial2d_io;

void RSSTMHD2D::check_rsst_parameters() const
{
  if (!std::isfinite(cs_bnd) || cs_bnd <= 0.0 ||
      !std::isfinite(ma_bnd) || ma_bnd <= 0.0){
    std::fprintf(stderr,
                 "Invalid RSST parameters: cs_bnd=%g, ma_bnd=%g "
                 "(both must be finite and positive)\n",
                 cs_bnd,ma_bnd);
    std::exit(EXIT_FAILURE);
  }
}

void RSSTMHD2D::rsst_(int i)
{				// Set RSST factor 1/\xi
  double cs=v_snd(i);		// Sound velocity
  double ca=v_alf(i);		// Alfven velocity
  double vv=sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]); // |V|
  double cup=std::max(vv/ma_bnd,std::max(2.0*ca,cs_bnd)); // Empirical safety margin from RSST-HLLD intermediate-state tests.
  ixi[i]=std::min(cup/cs,1.0);
}

void RSSTMHD2D::setdt(int flg)
{
  // Set time step to meet CFL, if flg is set
  if (flg){
    double vtmp=0.0,vmax=1.0;
    double vfas[3];
    for (int j=yoff;j<ny-yoff;j++){
      for (int i=xoff;i<nx-xoff;i++){
	int ss=nx*j+i;
	cx[ss]=0.5*(bx[ss]+bx[nx*j+(i+stxs[4])]);
	cy[ss]=0.5*(by[ss]+by[nx*(j+stys[5])+i]);
	prmtv(ss);
	rsst_(ss);
	vfas[0]=v_snd(ss);
	vfas[0]*=ixi[ss];	// RSST
	vfas[1]=v_alf(ss);
	vfas[2]=sqrt(vfas[0]*vfas[0]+vfas[1]*vfas[1]);
	vtmp=sqrt(vx[ss]*vx[ss]+vy[ss]*vy[ss])+vfas[2];
	if (!std::isfinite(vtmp)) abort_run("Nonfinite wave speed in setdt().");
	if (vtmp > vmax) vmax=vtmp;
      }
    }
    dt=cfl*dr/vmax;
  }
  initialize_dt();
}

RSSTMHD2D::RSSTMHD2D()
{
  // constructor  
  ixi=new double[nd];
  std::fill(ixi,ixi+nd,1.0);
  cs_bnd=1e15;
  ma_bnd=1e-15;
}

RSSTMHD2D::~RSSTMHD2D()
{
  // destructor  
  delete[] ixi;
}
