#include "boundary.h"
#include "common_mpi.h"
#include "mpi.h"
#include <stdlib.h>
#include "mhd_func.h"

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys,
	      int mpi_rank, int mpi_numx, int mpi_numy)
{
  int i;
  mpi_sdrv2d(p,nm,nx,ny,xoff,yoff,dnxs[0],dnys[0],mpi_rank,mpi_numx,mpi_numy);
  /* Note: If dnx(y)s[0] == 0, boundary is periodic*/
  for (i=0;i<nm;i++){
    mpi_xbc2d(p[i],nx,ny,xoff,yoff,stxs[i],dnxs[i],mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(p[i],nx,ny,xoff,yoff,stys[i],dnys[i],mpi_rank,mpi_numx,mpi_numy);
  }
}

void injection(double *ro, double *mx, double *my, double *mz,
	       double *bx, double *by, double *bz, double *en,
	       double ro_3, double dro3, double vx_1, double vy_1, double vz_1,
	       double bx_1, double by_1, double bz_1, double pr_1,
	       int nx, int ny, int xoff, int yoff, double gamma,
	       int mpi_rank, int mpi_numx, int mpi_numy)
// Set upper boundary for RMI
// ro_3: CD density @ upper boundary
// dro3: CD density perturbation @ upper boundary
// v{x,y,z}_1, b{x,y,z}_1, pr_1: upstream variables
{
  int i,j,ss;
  static int irflg=0;
  
  if (irflg == 0){
    double stim;			// Time at node 0, used for seed of random
    unsigned seed;
    if (mpi_rank == 0) stim=MPI_Wtime();
    MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    seed=(unsigned)(stim*(1+mpi_rank));
    srandom(seed);
    irflg=1;
  }

  if ((mpi_rank / mpi_numx) == (mpi_numy-1)){
    for (j=ny-yoff;j<ny;j++){
      int ctyflg=(j == (ny-yoff));
      for (i=0;i<nx;i++){
	double vx,vy,vz,cx,cy,cz,pr;

	ss=nx*j+i;
	
	ro[ss]=ro_3;
	ro[ss]+=dro3*((double)random()/RAND_MAX-0.5)*2.0; /* Random perturbation */

	vx=vx_1;
	vy=vy_1;
	vz=vz_1;

	bx[ss]=cx=bx_1;
	by[ss]=ctyflg?by[ss]:by_1; // Due to CT spacing
	cy=by_1;
	bz[ss]=cz=bz_1;

	pr=pr_1;

	mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gamma);      
      }
    }
  }
}
  
