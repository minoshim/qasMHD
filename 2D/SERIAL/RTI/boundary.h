#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

void boundary(double *p[], int nm, int nx, int ny, int xoff, int yoff,
	      int *stxs, int *dnxs, int *stys, int *dnys,
	      double gamma, const double *phi_g);

#endif
