#include "rmi2d_class.hpp"

#include <cstddef>
#include <stdexcept>

void RMI2D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

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

    // Refresh once per time step; ideal() reuses the inflow at every RK stage.
    bound(val,nm,stxs,dnxs,stys,dnys,true);
    ideal(dt_step);

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

void RMI2D::bound(double *values[], int nvars,
                  const int stxs[], const int dnxs[], const int stys[], const int dnys[])
{
  bound(values,nvars,stxs,dnxs,stys,dnys,false);
}

void RMI2D::bound(double *values[], int nvars,
                  const int stxs[], const int dnxs[], const int stys[], const int dnys[],
                  bool refresh_inflow)
{
  MHD2D::bound(values,nvars,stxs,dnxs,stys,dnys);

  // Identify the actual arrays, not just their number.
  bool is_state=(nvars == nm);
  if (is_state){
    for (int m=0;m<nm;m++){
      if (values[m] != val[m]) is_state=false;
    }
  }
  if (is_state){
    apply_inflow(refresh_inflow);
  } else if (nvars == 2 && values[0] == cx && values[1] == cy){
    apply_inflow_center();
  }
}

void RMI2D::apply_inflow(bool refresh_inflow)
{
  const std::size_t required_size=static_cast<std::size_t>(nx)*yoff;
  if (refresh_inflow){
    if (inflow_ro.size() != required_size) inflow_ro.resize(required_size);
    for (std::size_t s=0;s<required_size;s++){
      inflow_ro[s]=inflow.rho+inflow.rho_noise*((double)random()/RAND_MAX-0.5)*2.0;
    }
  } else if (inflow_ro.size() != required_size){
    throw std::logic_error("RMI2D inflow is not initialized; call init_() first");
  }

  for (int j=ny-yoff;j<ny;j++){
    const bool ct_face=(j == ny-yoff);
    for (int i=0;i<nx;i++){
      const int ss=nx*j+i;
      ro[ss]=inflow_ro[nx*(j-(ny-yoff))+i];
      vx[ss]=inflow.vx;
      vy[ss]=inflow.vy;
      vz[ss]=inflow.vz;
      pr[ss]=inflow.pressure;
      bx[ss]=inflow.bx;
      // The physical upper By face evolves through CT and must not be reset.
      if (!ct_face) by[ss]=inflow.by;
      bz[ss]=inflow.bz;
      cx[ss]=bx[ss];
      cy[ss]=by[ss];
      cnsvt(ss);
    }
  }
}

void RMI2D::apply_inflow_center()
{
  // Preserve the existing inflow magnetic closure without changing the state or RNG.
  for (int j=ny-yoff;j<ny;j++){
    for (int i=0;i<nx;i++){
      const int ss=nx*j+i;
      cx[ss]=bx[ss];
      cy[ss]=by[ss];
    }
  }
}
