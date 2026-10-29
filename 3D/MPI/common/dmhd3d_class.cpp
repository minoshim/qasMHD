#include "dmhd3d_class.hpp"
#include "mhd3d_io.hpp"

#include <cmath>
#include <limits>

using namespace mpi3d_io;

void DMHD3D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time
  // Dissipation solver called
  
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

    bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
    ideal(dt_step);

    double dcmax=setdc();
    if (!std::isfinite(dcmax) || dcmax < 0.0) abort_run("Invalid dissipation coefficient.");
    if (dcmax > 1e-15){
      // 3D counterpart of the 2D safety coefficient 6; covers fourth-order differences.
      const double substeps=9.0*(dcmax*dt_step)/(dr*dr);
      if (!std::isfinite(substeps) || substeps >= std::numeric_limits<int>::max()){
        abort_run("Dissipation substep count exceeds the supported integer range.");
      }
      int nc,ncmax=1+(int)(substeps);
      for (nc=0;nc<ncmax;nc++){	// sub-cycling
	bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
	dsptv(dt_step/ncmax);
      }
    }
    
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

DMHD3D::DMHD3D(int* argc, char*** argv, int mnp) : MHD3D(argc,argv,mnp)
{
  // constructor  
  nu=new double[nd];
  eta=new double[nd];
}

DMHD3D::~DMHD3D()
{
  // destructor  
  delete[] nu;
  delete[] eta;
}
