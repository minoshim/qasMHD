#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys);
void injection(double *ro, double *mx, double *my, double *mz,
	       double *bx, double *by, double *bz, double *en,
	       double ro_3, double dro3, double vx_1, double vy_1, double vz_1,
	       double bx_1, double by_1, double bz_1, double pr_1,
	       int nx, int ny, int xoff, int yoff, double gamma);

#endif
