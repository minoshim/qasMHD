#include "boundary.h"
#include "common_mpi.h"

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys,
	      double gamma, const double *phi_g,
	      int mpi_rank, int mpi_numx, int mpi_numy)
{
  int i;
  mpi_sdrv2d(p,nm,nx,ny,xoff,yoff,dnxs[0],dnys[0],mpi_rank,mpi_numx,mpi_numy);
  /* Note: If dnx(y)s[0] == 0, boundary is periodic*/
  for (i=0;i<nm;i++){
    mpi_xbc2d(p[i],nx,ny,xoff,yoff,stxs[i],dnxs[i],mpi_rank,mpi_numx,mpi_numy);
    mpi_ybc2d(p[i],nx,ny,xoff,yoff,stys[i],dnys[i],mpi_rank,mpi_numx,mpi_numy);
  }

  /* Manually set Y boundary (except energy) */
  /* If dny==+1, all boudary cells have the same value */
  /* If dny==-1, all boudary cells have zero value */
  int j,m;
  for (m=0;m<nm-1;m++){
    if (dnys[m] != 0){
      if (mpi_rank/mpi_numx == 0){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++){
	    p[m][nx*j+i]=0.5*(1.0+dnys[m])*p[m][nx*yoff+i];
	  }
	}
      }
      if (mpi_rank/mpi_numx == (mpi_numy-1)){
	for (j=ny-yoff+stys[m];j<ny;j++){
	  for (i=0;i<nx;i++){
	    p[m][nx*j+i]=0.5*(1.0+dnys[m])*p[m][nx*(ny-yoff+stys[m]-1)+i];
	  }
	}
      }
    }
  }

  /* Y boundary for energy needs special care */
  if (nm == 8){
    int ss,sb;
    double fac=(2.0-gamma)/(gamma-1.0);
    double *ro=p[0],*en=p[7];
    if (mpi_rank/mpi_numx == 0){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sb=nx*yoff+i;
	  en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	}
      }
    }
    if (mpi_rank/mpi_numx == (mpi_numy-1)){
      for (j=ny-yoff;j<ny;j++){
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sb=nx*(ny-yoff-1)+i;
	  en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	}
      }
    }
  }  
}
