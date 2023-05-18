#ifndef _MHD_FUNC_H_
#define _MHD_FUNC_H_

/* Equation of State for Ideal MHD */
double mhd_ieos(double ro, double vx, double vy, double vz,
		double bx, double by, double bz, double en,
		double gamma);

/* MHD total energy */
double mhd_energy(double ro, double vx, double vy, double vz,
		  double bx, double by, double bz, double pr,
		  double gamma);

/* MHD primitive variables */
void mhd_prmtv(double ro, double mx, double my, double mz,
	       double bx, double by, double bz, double en,
	       double *vx, double *vy, double *vz, double *pr,
	       double gamma);

/* MHD conservative variables */
void mhd_cnsvt(double ro, double vx, double vy, double vz,
	       double bx, double by, double bz, double pr,
	       double *mx, double *my, double *mz, double *en,
	       double gamma);

/* Weight for CUCT (Minoshima+19, Eq. (42)) */
/* ro,vx,vy @ cell center, bx and by @ cell face */
/* xoffset and yoffset are (1,nx) for x-y plane, (nx,nx*ny) for y-z plane, (nx*ny,1) for z-x plane */
double mhd_cuct_weight(const double *ro, const double *vx, const double *vy,
		       const double *bx, const double *by,
		       int xoffset, int yoffset);

/* Resistive E-field (eta*j) @ CT grid */
/* B are defined @ Bx(i-1/2,j,k), By(i,j-1/2,k), Bz(i,j,k-1/2) */
/* E are defined @ Ex(i,j-1/2,k-1/2), Ey(i-1/2,j,k-1/2), Ez(i-1/2,j-1/2,k) */
/* eta is resistivity defined @ cell center (i,j,k) */
/* idx=1.0/dx */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
void mhd_ct_eres(const double *bx, const double *by, const double *bz,
		 const double *eta,
		 double idx, double idy, double idz,
		 int xoffset, int yoffset, int zoffset,
		 double (*func_df)(const double *f),
		 double *ex, double *ey, double *ez);

/* MHD left and right states at cell interface */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
/* func_lr = interpolation function */
void mhd_lrstate(const double *ro, const double *vx, const double *vy, const double *vz,
		 const double *by, const double *bz, const double *pr,
		 double bx, double gamma, int offset,
		 void (*func_lr)(const double*, double*, double*),
		 double *vl, double *vr);

/* Left and right state of numerical flux of B, by*vx-bx*vy, required for CUCT  */
/* Bx @ cell face, By @ cell center */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
void mhd_lr_fb(const double *vx, const double *vy, const double *bx, const double *by,
	       int offset,
	       void (*func_lr)(const double*, double*, double*),
	       double *vl, double *vr);

/* Left and right state of single variable  */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
void mhd_lr_single(const double *f, int offset,
		   void (*func_lr)(const double*, double*, double*),
		   double *vl, double *vr);

/* Update 1D MHD variables */
/* func_df = 1st derivative */
void mhd_updt1d(double *val, double val0, const double *fx,
		double dtdx, const double rk_fac[2], int xoffset,
		double (*func_df)(const double*));

/* Update 2D MHD variables @ cell center */
void mhd_updt2d(double *val, double val0,
		const double *fx, const double *fy,
		double dtdx, double dtdy, const double rk_fac[2],
		int xoffset, int yoffset,
		double (*func_df)(const double*));

/* Update 3D MHD variables @ cell center */
void mhd_updt3d(double *val, double val0,
		const double *fx, const double *fy, const double *fz,
		double dtdx, double dtdy, double dtdz, const double rk_fac[2],
		int xoffset, int yoffset, int zoffset,
		double (*func_df)(const double*));

/* Update 2D MHD B @ cell face (bx or by) by CT method */
/* dtdy should include the direction. For exsample, it is negative when val=by */
void mhd_updt2d_ctb(double *val, double val0, const double *ez,
		    double dtdy, const double rk_fac[2], int offset,
		    double (*func_df)(const double*));

/* Update 3D MHD Bx @ cell face by CT method */
/* yoffset and zoffset are nx and nx*ny */
void mhd_updt3d_ctb(double *bx, double bx0, const double *ey, const double *ez,
		    double dtdy, double dtdz, const double rk_fac[2],
		    int yoffset, int zoffset,
		    double (*func_df)(const double*));



/* Below is obsolete and will be deleted */

/* Update 3D MHD variables @ cell center */
/* Flag = 1 or 0. Update ith-variable if flag[i]=1 */
/* xoffset,yoffset,zoffset are 1,nx,nx*ny */
void mhd_updt3d(double *val[], const double *val0,
		const double *fx, const double *fy, const double *fz,
		double dtdx, double dtdy, double dtdz, const double rk_fac[2], const double flag[],
		int xoffset, int yoffset, int zoffset,
		double (*func_df)(const double*));

#endif
