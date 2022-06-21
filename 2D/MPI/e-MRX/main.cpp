#include "global.hpp"
using namespace global;

#include "new_delete.hpp"
#include "init.hpp"
#include "get.hpp"
#include "dataio.hpp"

int main(int argc, char* argv[])
{
  // MPI Send/Recv arguments
  int mpi_rank;
  int mpi_num;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
  if (mpi_num != mnp){
    if (mpi_rank == 0) puts("MPI number of process is not correct. Check mnp in global.hpp");
    MPI_Finalize();
    return 0;
  }

  int n=0,cnt=0;
  double tim=0,trec=dtrec;
  double stim,etim;

  new_delete();
  init_grid(mpi_rank);
  init_plasma(mpi_rank);
  double *p[]={ro,mx,my,mz,bx,by,bz,en};
  
  get_dt(&dt,cfl,dr,!(mpi_rank));

  bkup_load(p,8,nd,&n,&cnt,&tim,&dt,&trec,mpi_rank,fildir);
  if (n == 0) dataio(p,8,nd,n,cnt,tim,0,mpi_rank);

  stim=MPI_Wtime();
  // Time integration
  while(n++ < nmax && tim < tend){
    tim+=dt;

    // Boundary condition
    boundary(p,8,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);

    // MHD update
    mhd_fd2d(p,dt,dx,dy,de,mpme,8,nx,ny,xoff,yoff,gam,mpi_rank,mpi_numx,mpi_numy);
    
#if (DIFF)			// Diffusion
    double dcmax=get_dcoefs(nu,eta);
    int nc,ncmax=1+(int)((dcmax*dt)/(dr*dr)/0.5);
    for (nc=0;nc<ncmax;nc++){	// sub-cycling
      boundary(p,8,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);
      mhd_diff2d(p,dt/ncmax,dx,dy,nu,eta,8,nx,ny,xoff,yoff,gam,mpi_rank,mpi_numx,mpi_numy);
    }
#endif

    // Re-calculate time step
    if (n % 2 == 0) get_dt(&dt,cfl,dr,0);

    // Output
    if (tim >= trec){
      trec+=dtrec;
      cnt++;
      bkup_save(p,8,nd,n,cnt,tim,dt,trec,mpi_rank,fildir);
      dataio(p,8,nd,n,cnt,tim,!(mpi_rank),mpi_rank);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  etim=MPI_Wtime();
  if (mpi_rank == 0) printf("Elapse time = %lu sec.\n",(unsigned long)(etim-stim));

  new_delete();
  MPI_Finalize();
  return 0;
}
  