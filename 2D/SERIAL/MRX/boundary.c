#include "boundary.h"
#include "common_func.h"

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys)
{
  int i;
  for (i=0;i<nm;i++){
    bc2d(p[i],nx,ny,xoff,yoff,stxs[i],dnxs[i],stys[i],dnys[i]);
  }
}
