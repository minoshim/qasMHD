#ifndef _MHD_EIGEN_H_
#define _MHD_EIGEN_H_

/* Eigenvalues and eigenvectors for non-conservative form, based on ATHENA (Stone+2008) */
void mhd_eigen_st(double ro, double vx, double vy, double vz, double bx, double by, double bz, double pr, double gamma,
		  double *al1, double *al2, double *al3, double *al4, double *al5, double *al6, double *al7,
		  double *l1, double *l2, double *l3, double *l4, double *l5, double *l6, double *l7, 
		  double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7);

/* Approximate eigenvectors for non-conservative form, based on MHD subsystem (Miyoshi+2020) */
void mhd_eigen_mi(double ro, double vx, double vy, double vz, double bx, double by, double bz, double pr, double gamma,
		  double *al1, double *al2, double *al3, double *al4, double *al5, double *al6, double *al7,
		  double *l1, double *l2, double *l3, double *l4, double *l5, double *l6, double *l7, 
		  double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7);

/* Approximate eigenvectors for non-conservative form, based on Alfven characteristics (Minoshima+2020) */
void mhd_eigen_ae(double ro, double vx, double vy, double vz, double bx, double by, double bz, double pr, double gamma,
		  double *al1, double *al2, double *al3, double *al4, double *al5, double *al6, double *al7,
		  double *l1, double *l2, double *l3, double *l4, double *l5, double *l6, double *l7, 
		  double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7);

/* Reconstruction of characteristic variables */
void mhd_c_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
		   const double *by, const double *bz, const double *pr,
		   double bx, double gamma, const int s,
		   double *ul, double *ur,
		   void (*func_lr)(const double*, double*, double*));

/* Reconstruction of fast and Alfven characteristic variables (Miyoshi+2020) */
void mhd_m_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
		   const double *by, const double *bz, const double *pr,
		   double bx, double gamma, const int s,
		   double *ul, double *ur,
		   void (*func_lr)(const double*, double*, double*));

/* Reconstruction of Alfven characteristic variables (Minoshima+2020) */
void mhd_a_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
		   const double *by, const double *bz, const double *pr,
		   double bx, double gamma, const int s,
		   double *ul, double *ur,
		   void (*func_lr)(const double*, double*, double*));

/* Reconstruction of primitive variables */
void mhd_reconst(const double *ro, const double *vx, const double *vy, const double *vz,
		 const double *by, const double *bz, const double *pr,
		 double bx, double gamma, const int s,
		 double *ul, double *ur,
		 void (*func_lr)(const double*, double*, double*));

#endif
