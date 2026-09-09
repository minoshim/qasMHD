#include "dmhd2d_class.hpp"

void DMHD2D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time
  // Dissipation solver called
  
  if (n == 0) dout_(0);
  
  clock_t stim=clock();
  while(flg ? (tim < tmax) : (n < nmax)){
    ++n;
    const double remaining=tmax-tim;
    const bool last_step=flg && dt >= remaining;
    const double dt_step=last_step ? remaining : dt;
    if (flg){
      tim=last_step ? tmax : tim+dt_step;
    } else{
      tim=(n == nmax) ? tmax : n*dt;
    }

    bound(val,nm,stxs,dnxs,stys,dnys);
    ideal(dt_step);

    double dcmax=setdc();
    if (dcmax > 1e-15){
      int nc,ncmax=1+(int)(6*(dcmax*dt_step)/(dr*dr));
      for (nc=0;nc<ncmax;nc++){	// sub-cycling
	bound(val,nm,stxs,dnxs,stys,dnys);
	dsptv(dt_step/ncmax);
      }
    }
    
    setdt(flg*(n % 2 == 0));
    
    if (flg ? (tim >= trec) : ((n % nrec) == 0)){
      cnt++;
      trec=(cnt+1)*dtrec;
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
