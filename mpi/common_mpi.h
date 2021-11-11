#ifndef _COMMON_MPI_H_
#define _COMMON_MPI_H_

/* MPI SendRecv for 2D variables */
/* Set dn=0 for Periodic boundary */
/* For other condition, call mpi_xbc2d and mpi_ybc2d later */
void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff,
		int dnx, int dny,
		int mpi_rank, int mpi_numx, int mpi_numy);

/* 2D X and Y BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);
void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy);

/* MPI SendRecv for 3D variables */
/* Set dn=0 for Periodic boundary */
/* For other condition, call mpi_x(y,z)bc3d later */
void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int dnx, int dny, int dnz,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

/* 3D X,Y,Z BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);
void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

#endif
