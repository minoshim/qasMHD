#ifndef _MHD_FLUX_H_
#define _MHD_FLUX_H_

/* MLAU solver (Minoshima+2020) */
void calc_flux_mlau(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		    double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		    double bnc, double gamma, const double *dv,
		    double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen);

/* Roe solver */
void calc_flux_roe(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		   double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		   double bnc, double gamma, const double *dv,
		   double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen);

/* HLLD solver (Miyoshi+2005) */
void calc_flux_hlld(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		    double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		    double bnc, double gamma, const double *dv,
		    double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen);

/* LHLLD solver (Minoshima+2021) */
void calc_flux_lhlld(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		     double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		     double bnc, double gamma, const double *dv,
		     double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen);

/* Viscous flux */
/* v{x,y,z} are defined @ cell center (i,j,k) */
/* rnu = viscosity coefficient (ro*nu) at the interface i-1/2, idx=1.0/dx */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
/* fmx,fmy,fmz,fen: output fluxes at the interface for momentum and energy */
void calc_flux_viscous(const double *vx, const double *vy, const double *vz,
		       double rnu, double idx, double idy, double idz,
		       int xoffset, int yoffset, int zoffset,
		       double (*func_fc)(const double *f),
		       double (*func_df)(const double *f),
		       double *fmx, double *fmy, double *fmz, double *fen);

/* Resistive flux for CT grid spacing */
/* B are defined @ By(i,j-1/2,k), Bz(i,j,k-1/2) */
/* E are defined @ Ey(i-1/2,j,k-1/2), Ez(i-1/2,j-1/2,k), prepared through MHD_CT_ERES */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
/* fby,fbz,fen: output fluxes at the interface for B and energy */
void calc_flux_resistive(const double *by, const double *bz, const double *ey, const double *ez,
			 int xoffset, int yoffset, int zoffset,
			 double (*func_fc)(const double *f),
			 double *fby, double *fbz, double *fen);

#endif
