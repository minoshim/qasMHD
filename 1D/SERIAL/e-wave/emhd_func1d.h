#ifndef _EMHD_FUNC1D_H_
#define _EMHD_FUNC1D_H_

/* Electron MHD convertion of B and Q, satisfying (1-de^2 \nabla^2) B = Q, where de is electron inertia length */
void emhd_b2q(double *b, double *q, double de, double dx, int nx, int xoff, int dnx); /* B => Q */
void emhd_q2b(double *q, double *b, double de, double dx, int nx, int xoff, int dnx); /* Wrapper of Q => B functions */
int emhd_q2b_cg(double *q, double *b, double de, double dx, int nx, int xoff, int dnx); /* Q => B by Conjudate Gradient */
int emhd_q2b_gs(double *q, double *b, double de, double dx, int nx, int xoff, int dnx); /* Q => B by Gauss-Seidel */

#endif
