#include "dmhd2d_class.hpp"

double DMHD2D::setdc()
{
  // Set dissipation coefficients and return their maximum.
  double dcmax=0.0,dtmp;
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      int ss=nx*j+i;
      
      nu[ss]=nu0;
      eta[ss]=eta0;

      // Localized dissipation
      // double r=sqrt(x[i]*x[i]+y[j]*y[j]);
      // nu[ss]=nu0*exp(-r/1.0);
      // eta[ss]=eta0*exp(-r/1.0);
      
      dtmp=max(2*nu[ss],eta[ss]); // Factor 2 is multiplied in viscous coef. for robust estimation
      if (dtmp > dcmax) dcmax=dtmp;
    }
  }
  return dcmax;
}

void DMHD2D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time
  // Dissipation solver called
  
  if (n == 0) dout_(0);
  
  clock_t stim=clock();
  while(n++ < nmax && tim < tmax){
    tim+=dt;

    bound(val,nm,stxs,dnxs,stys,dnys);
    ideal(dt);

    double dcmax=setdc();
    if (dcmax > 1e-15){
      int nc,ncmax=1+(int)(2*(dcmax*dt)/(dr*dr));
      for (nc=0;nc<ncmax;nc++){	// sub-cycling
	bound(val,nm,stxs,dnxs,stys,dnys);
	dsptv(dt/ncmax);
      }
    }
    
    setdt(flg*(n % 2 == 0));
    
    if (tim >= trec){
      cnt++;
      trec+=dtrec;
      dout_(1);
    }
  }
  clock_t etim=clock();
  printf("CPU time = %f sec.\n",(double)(etim-stim)/CLOCKS_PER_SEC);

}

DMHD2D::DMHD2D()
{
  // constructor  
  nu=new double[nd];
  eta=new double[nd];
}

DMHD2D::~DMHD2D()
{
  // destructor  
  delete[] nu;
  delete[] eta;
}
