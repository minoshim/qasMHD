#ifndef _COMMON_FUNC_H_
#define _COMMON_FUNC_H_

double max(double a, double b);	/* return max(a,b) */

double min(double a, double b);	/* return min(a,b) */

double minmod(double a, double b); /* return minmod(a,b) */

double minmod3(double a, double b, double c); /* return minmod(a,b,c) */

double rand_noise(const double *params, unsigned seed); /* return uniform random distribution (params[0] +- paramas[1]) */

void cpy_array(double *a, const double *b, int n); /* copy b to a. n = length */

void conv_d2f(float *valo, const double *vali, int n); /* convert double (vali) to single (valo) presicion */

void conv_f2d(double *valo, const float *vali, int n); /* convert single (vali) to double (valo) presicion */

void rk_updt(double *f1, double f0, double df, double rkfac0, double rkfac1);
/* TVD Runge-Kutta update. f1=rkfac0*f0+rkfac1*(f1+df) */

void bc1d(double *f, int nx, int xoff, int dnx);
/* 1D boundary condition */
/* xoff = number of ghost cell */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */

void bc2d(double *f, int nx, int ny, int xoff, int yoff,
	  int stx, int dnx, int sty, int dny);
/* 2D Boundary condition */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */

void bc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	  int stx, int dnx, int sty, int dny, int stz, int dnz);
/* 3D Boundary condition */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */

void bkup_load(double *p[], int nm, int nd, int *n, int *cnt, double *tim, double *dt, double *trec, int mpi_rank, const char *fildir);
void bkup_save(double *p[], int nm, int nd, int n, int cnt, double tim, double dt, double trec, int mpi_rank, const char *fildir);
/* Load and save backup data for restart */
/* p: Array of pointer for variables (ro,mx,my,...) */
/* nm: Number of variables (7 for MHD1D, 8 for multi-D MHD) */
/* nd: Number of spatial grid (nx*ny*nz) */
/* n,cnt: Numbers of iteration and output  */
/* tim,dt,trec: Simulation time, time step, and next time for output */
/* mpi_rank: # of MPI process. Set 0 for serial calc. */
/* fildir: data directory */

#endif

