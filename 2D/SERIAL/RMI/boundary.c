#include "boundary.h"
#include "common_func.h"
#include <time.h>
#include <stdlib.h>
#include "mhd_func.h"

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys)
{
  int i;
  for (i=0;i<nm;i++){
    bc2d(p[i],nx,ny,xoff,yoff,stxs[i],dnxs[i],stys[i],dnys[i]);
  }
}

void injection(double *ro, double *mx, double *my, double *mz,
	       double *bx, double *by, double *bz, double *en,
	       double ro_3, double dro3, double vx_1, double vy_1, double vz_1,
	       double bx_1, double by_1, double bz_1, double pr_1,
	       int nx, int ny, int xoff, int yoff, double gamma)
// Set upper boundary for RMI
// ro_3: CD density @ upper boundary
// dro3: CD density perturbation @ upper boundary
// v{x,y,z}_1, b{x,y,z}_1, pr_1: upstream variables
{
  int i,j,ss;
  static int irflg=0;
  
  if (irflg == 0){
    unsigned seed;
    seed=(unsigned)time(NULL);
    srandom(seed);
    irflg=1;
  }

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
