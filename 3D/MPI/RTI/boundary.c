#include "boundary.h"
#include "common_mpi.h"

void boundary(double *p[], int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	      int *stxs, int *dnxs, int *stys, int *dnys, int *stzs, int *dnzs,
	      double gamma, const double *phi_g,
	      int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
{
  int i;
  mpi_sdrv3d(p,nm,nx,ny,nz,xoff,yoff,zoff,dnxs[0],dnys[0],dnzs[0],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  /* Note: If dnx(y,z)s[0] == 0, boundary is periodic*/
  for (i=0;i<nm;i++){
    mpi_xbc3d(p[i],nx,ny,nz,xoff,yoff,zoff,stxs[i],dnxs[i],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(p[i],nx,ny,nz,xoff,yoff,zoff,stys[i],dnys[i],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(p[i],nx,ny,nz,xoff,yoff,zoff,stzs[i],dnzs[i],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  }

  /* Z boundary for energy needs special care */
  if (nm == 8){
    int j,k,ss,sb;
    int m_xy=mpi_numx*mpi_numy;
    double fac=(2.0-gamma)/(gamma-1.0);
    double *ro=p[0],*en=p[7];
    if ((mpi_rank/m_xy) == 0){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    sb=nx*(ny*zoff+j)+i;
	    en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	  }
	}
      }
    }
    if ((mpi_rank/m_xy) == (mpi_numz-1)){
      for (k=nz-zoff;k<nz;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    sb=nx*(ny*(nz-zoff-1)+j)+i;
	    en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	  }
	}
      }
    }
  }
}
