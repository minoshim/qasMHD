#ifndef _EMHD_FUNC2D_H_
#define _EMHD_FUNC2D_H_

/* Extended MHD conversion of Enew and Eorg, satisfying (ro-de^2 \nabla^2) Enew = ro*Eorg */
/* ro is number density and de is electron inertia length (both are normalized) */
void emhd_enew2eorg(double *enew, double *eorg, const double *ro,
		    double de, double dx, double dy,
		    int nx, int ny, int xoff, int yoff,
		    int stx, int dnx, int sty, int dny,
		    int mpi_rank, int mpi_numx, int mpi_numy); /* Enew => Eorg */
void emhd_eorg2enew(double *eorg, double *enew, const double *ro,
		    double de, double dx, double dy,
		    int nx, int ny, int xoff, int yoff,
		    int stx, int dnx, int sty, int dny,
		    int mpi_rank, int mpi_numx, int mpi_numy); /* Wrapper of Eorg => Enew functions */
int emhd_eorg2enew_cg(double *eorg, double *enew, const double *ro,
		      double de, double dx, double dy,
		      int nx, int ny, int xoff, int yoff,
		      int stx, int dnx, int sty, int dny,
		      int mpi_rank, int mpi_numx, int mpi_numy); /* Eorg => Enew by Conjugate Gradient method */
double calc_jcen(const double *bx, const double *by, double idx, double idy,
		 int xoffset, int yoffset);
double calc_dnvv(const double *ro,
		 const double *vpx, const double *vpy, const double *vpz,
		 const double *vex, const double *vey, const double *vez,
		 double idx, double idy, double idz,
		 int xoffset, int yoffset, int zoffset);

#endif
