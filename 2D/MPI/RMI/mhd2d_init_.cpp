#include "mhd2d_class.hpp"

void MHD2D::init_()
{
  // RM instability, including injection boundary condition. 
  
  int i,j;
  static int flg=0;

  // Shock upstream paramters
  double ro_1=+1.0;	// Density
  static double vx_1=0.0;	// Velocity in X
  static double vy_1=-1.0;	// Velocity in Y
  static double vz_1=0.0;	// Velocity in Z
  double ma_u=10.0;	// Sound Mach number
  static double pr_1=(ro_1*vy_1*vy_1)/(gam*ma_u*ma_u); // Pressure
  double beta=1e8;				// Beta
  double b0_1=(MAGNET)*sqrt(2.0*pr_1/beta);	// Magnetic field
  static double bx_1=b0_1;
  static double by_1=0.0;
  static double bz_1=0.0;
  static double vref=-0.6;	// Velocity of reference. 0 = shock rest frame
  double dw=0.5*dy;	// Shock width
  //// Contact Discon parameters
  static double ro_3=ro_1*1e1;	// Density jump
  static double dro3=RANDOM*ro_3*0.1;	// Density perturbation (Available when RANDOM=1)
  double lambda=getlx();	// Wavelength of perturbation for position
  double psi=0.1*lambda;	// Amplitude of perturbation for position
  double mmax=8;		// Number of modes (Available when RANDOM=1)
  
  if (flg != 0){		// Start injection
    // Injection @ upper boundary
    if (mpi_rany == (mpi_numy-1)){
      for (j=ny-yoff;j<ny;j++){
	int ctyflg=(j == (ny-yoff));
	for (i=0;i<nx;i++){
	  int ss=nx*j+i;

	  ro[ss]=ro_3+dro3*((double)random()/RAND_MAX-0.5)*2.0; /* Random perturbation */
	  vx[ss]=vx_1;
	  vy[ss]=vy_1;
	  vy[ss]-=vref;
	  vz[ss]=vz_1;
	  pr[ss]=pr_1;
	  bx[ss]=bx_1;
	  by[ss]=ctyflg?by[ss]:by_1; // Due to CT spacing
	  bz[ss]=bz_1;
	  cx[ss]=bx[ss];
	  cy[ss]=by[ss];

	  cnsvt(ss);
	}
      }
    }
  }	 // End of injection
  else {			// Start initialization
    flg=1;
    int isum=0,jsum=0,m;
    for (m=0;m<mpi_ranx;m++){
      isum+=(XMESH+m)/mpi_numx;
    }
    for (m=0;m<mpi_rany;m++){
      jsum+=(YMESH+m)/mpi_numy;
    }
    for (i=0;i<nx;i++) x[i]=(i-xoff+isum+0.5)*dx+xmin;
    for (j=0;j<ny;j++) y[j]=(j-yoff+jsum+0.5)*dy+ymin;
    
    // Shock downstream parameters satisfies R-H relation
    double gp1=gam+1.0;
    double gm1=gam-1.0;
#if (MAGNET)
    // MHD shock
    double a=2.0*(2.0-gam)/beta;
    double b=gam*(gm1*ma_u*ma_u+2.0/beta+2.0);
    double c=-gam*gp1*ma_u*ma_u;
    // Compression ratio
    double r=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    // Pr2/Pr1
    double rr=gam*ma_u*ma_u*(1.0-1.0/r)-(r*r-1.0)/beta+1.0;
#else
    // HD shock
    double r=(gp1*ma_u*ma_u)/(gm1*ma_u*ma_u+2.0); // Compression ration
    double rr=(2.0*gam*ma_u*ma_u-gm1)/gp1;	// Pr2/Pr1
#endif
    // Donwstream paramters
    double ro_2=ro_1*r;
    double vy_2=vy_1/r;
    double bx_2=bx_1*r;
    double bz_2=bz_1*r;
    double pr_2=pr_1*rr;

    // Messaage of initial condition
    if (mpi_rank == 0){
      puts("Initial condition in shock rest frame");
      printf("U: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_1,vx_1,vy_1,vz_1,bx_1,by_1,bz_1,pr_1);
      printf("D: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_2,vx_1,vy_2,vz_1,bx_2,by_1,bz_2,pr_2);
      printf("Mach = %.9f\n",ma_u);
    }
  
    // Position of contact discon.
    double stim;			// Time at node 0, used for seed of rand_noise
    unsigned seed;
    if (mpi_rank == 0) stim=MPI_Wtime();
    MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    seed=(unsigned)(stim*(1+mpi_rank));
    srandom(seed);
    double ypos[nx];
    for (i=0;i<nx;i++){
      ypos[i]=1.0+psi*cos(2*M_PI*x[i]/lambda); // Single mode
    }
#if (RANDOM)
    // Multi-mode
    double para[2]={M_PI,M_PI};
    for (int m=2;m<=mmax;m++){
      double xphase;
      if (mpi_rank == 0) xphase=rand_noise(para,seed);
      MPI_Bcast(&xphase,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      for (i=0;i<nx;i++){
	ypos[i]+=(psi/m)*cos(2*M_PI*m*(x[i]-xphase)/lambda);
      }
    }
#endif
 
    for (j=0;j<ny;j++){
      double sfunc=0.5*(1.0+tanh(y[j]/dw)); // 0 (y<0), 1 (y>0)
      for (i=0;i<nx;i++){
	int ss=nx*j+i;
      
	ro[ss]=ro_2+(ro_1-ro_2)*sfunc;
	// Contact discon
	double rocd=ro_3+dro3*((double)random()/RAND_MAX-0.5)*2.0; // Random perturbation
	ro[ss]+=(rocd-ro_1)*0.5*(1.0+tanh((y[j]-ypos[i])/dw));

	vx[ss]=vx_1;
	vy[ss]=vy_2+(vy_1-vy_2)*sfunc;
	vy[ss]-=vref;
	vz[ss]=vz_1;
	pr[ss]=pr_2+(pr_1-pr_2)*sfunc;
	bx[ss]=bx_2+(bx_1-bx_2)*sfunc;
	by[ss]=by_1;
	bz[ss]=bz_2+(bz_1-bz_2)*sfunc;
	cx[ss]=bx[ss];
	cy[ss]=by[ss];

	// In 2D, cell center Bz is identical to cell edge Bz.
	// Thus cz is not explicitly initialized here.
	// They should be defined in constructer (see mhd2d_class.cpp).
	cnsvt(ss);
      }
    }

    // Boundary condition
    bound(val,nm,stxs,dnxs,stys,dnys);
  } // End of initialization

}

