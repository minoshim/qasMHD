#include "hmhd1d_class.hpp"

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
    hall_(dt);
    
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
