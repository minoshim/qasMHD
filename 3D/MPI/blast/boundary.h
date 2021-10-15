#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

void boundary(double *p[], int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	      int *stxs, int *dnxs, int *stys, int *dnys, int *stzs, int *dnzs,
	      int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

#endif
