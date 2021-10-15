#include "boundary.h"
#include "common_mpi.h"

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
