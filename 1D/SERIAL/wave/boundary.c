#include "boundary.h"
#include "common_func.h"

void boundary(double *p[], int n, int nx, int xoff, int *dnxs)
{
  int i;
  for (i=0;i<n;i++){
    bc1d(p[i],nx,xoff,dnxs[i]);
  }
}
