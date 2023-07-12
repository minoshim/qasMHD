#include "hmhd1d_class.hpp"

double HMHD1D::haldt()
{
  // Return time step for Hall term to satisfy CFL condition
  double vtmp=0.0,vmax=1.0;
  for (int i=xoff;i<nx-xoff;i++){
    hallv(i);
    vtmp=fabs(hx[i])+vphix*fabs(cx[i]/ro[i]);
    if (vtmp > vmax) vmax=vtmp;
  }
  return cfl*dx/vmax;
}

void HMHD1D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  if (n == 0) dout_();
  
  while( (flg)?(n++ < nmax && tim < tmax):(n++ < nmax) ){
    tim+=dt;

    bound(val,nm,dnxs);

    // Syb-cycling of Hall term
    int hmax=1+(int)(dt/haldt());
    for (int h=0;h<hmax;h++){
      hall_(dt/hmax);
      bound(val,nm,dnxs);
    }
    
    ideal(dt);
    setdt(flg*(n % 2 == 0));
    
    if ( (flg)?(tim >= trec):((n % nrec) == 0) ){
      cnt++;
      trec+=dtrec;
      dout_();
    }
  }
}

HMHD1D::HMHD1D()
{
  // constructor  
  hx=new double[nx];
  hy=new double[nx];
  hz=new double[nx];
}

HMHD1D::~HMHD1D()
{
  // destructor  
  delete[] hx;
  delete[] hy;
  delete[] hz;
}
