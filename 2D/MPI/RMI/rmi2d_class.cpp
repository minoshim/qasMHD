#include "rmi2d_class.hpp"
#include "mhd2d_io.hpp"

#include <cstddef>

using namespace mpi2d_io;

void RMI2D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  prepare_run(flg); // Load backup data without reinitializing its dt
  if (n == 0) dout_(0);

  double stim=MPI_Wtime();
  while(flg ? (tim < tmax) : (n < nmax)){
    check_step();
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
      dout_(!mpi_rank);
      // Do not advance the checkpoint until every rank has finished its output.
      check_mpi(MPI_Barrier(MPI_COMM_WORLD),"Output barrier failed.");
      bkup_(1); // Save backup data; wait for all ranks to succeed
    }
  }
  double etim=MPI_Wtime();
  if (!mpi_rank) printf("Elapse time = %lu sec.\n",(unsigned long)(etim-stim));
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
  if (mpi_rany != mpi_numy-1) return; // Only the global upper Y boundary injects.

  const std::size_t required_size=static_cast<std::size_t>(nx)*yoff;
  if (refresh_inflow){
    resize_work(inflow_ro,required_size);
    const double density_params[2]={inflow.rho,inflow.rho_noise};
    for (std::size_t s=0;s<required_size;s++){
      inflow_ro[s]=rand_noise_mt(density_params,noise_engine);
    }
  } else if (inflow_ro.size() != required_size){
    abort_run("RMI2D inflow is not initialized; call init_() first.");
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
  if (mpi_rany != mpi_numy-1) return;

  // Preserve the existing inflow magnetic closure without changing the state or RNG.
  for (int j=ny-yoff;j<ny;j++){
    for (int i=0;i<nx;i++){
      const int ss=nx*j+i;
      cx[ss]=bx[ss];
      cy[ss]=by[ss];
    }
  }
}
