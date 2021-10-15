#include "boundary.h"
#include "common_mpi.h"

void boundary(double *p[], int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	      int *stxs, int *dnxs, int *stys, int *dnys, int *stzs, int *dnzs,
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
}
