#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

void boundary(double *ro, double *mx, double *my, double *mz,
	      double *by, double *bz, double *en,
	      double bx,  double dt, double dx,
	      int nx, int xoff, double gamma,
	      int istt, int iend);

#endif
