#include "boundary.h"
#include "common_func.h"

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys,
	      double gamma, const double *phi_g)
{
  int i;
  for (i=0;i<nm;i++){
    bc2d(p[i],nx,ny,xoff,yoff,stxs[i],dnxs[i],stys[i],dnys[i]);
  }

  /* Y boundary for energy needs special care */
  if (nm == 8){
    int j,ss,sb;
    double fac=(2.0-gamma)/(gamma-1.0);
    double *ro=p[0],*en=p[7];
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	ss=nx*j+i;
	sb=nx*yoff+i;
	en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
      }
    }
    for (j=ny-yoff;j<ny;j++){
      for (i=0;i<nx;i++){
	ss=nx*j+i;
	sb=nx*(ny-yoff-1)+i;
	en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
      }
    }
  }  
}
