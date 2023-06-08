#include "dmhd3d_class.hpp"

inline double harris_field(double x, const double *params);
inline double harris_density(double x, const double *params);

const double lambda=1.0;	// Current sheet thickness

void DMHD3D::init_()
{
  // Magnetic reconnection
  int i,j,k;
  // Initial condition parameters
  const double beta=0.2;	// Plasma beta @ lobe
  const double ro0=1.0;		// Density @ CS
  const double ro1=0.2;		// Density @ lobe
  const double b0=1.0;		// Mag field @ lobe
  const double b1=0.05;		// Mag field perturbation by Zenitani
  const double bg=0.0;		// Guide mag field along Z
  const double dv=0.01;		// Random noize perturbation to Vy (avaiable when RANDOM=1)
  const double para[2]={0,lambda};
  const double dvpara[2]={0,dv};
  const double dbpara[2]={b1,b1};

  double dvy[nz][nx];
  unsigned seed;
  double stim;			// Time at node 0, used for seed of rand_noise
  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //Random seed depending on x and z, used for vy perturbation
  seed=(unsigned)(stim*(1+mpi_ranx+(mpi_numx*mpi_ranz)));
  
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      dvy[k][i]=0.0;
#if (RANDOM)
      dvy[k][i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#endif
    }
  }

  //Reset random seed depending on z, used for b1 perturbation
  seed=(unsigned)(stim*(1+mpi_ranz));
  srandom(seed);
  
  for (k=0;k<nz;k++){
    double b1amp=rand_noise(dbpara,seed); // Perturbation to b1 is independent of x nd y
    double zm=z[k]-0.5*dz;
    for (j=0;j<ny;j++){
      double ym=y[j]-0.5*dy;
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double xm=x[i]-0.5*dx;
	
	// ro[ss]=(ro0-ro1)*harris_density(y[j],para)+ro1;
	ro[ss]=ro0*harris_density(y[j],para)+ro1;

	vx[ss]=0.0;
	vy[ss]=0.0;
	vz[ss]=0.0;
	// Perturbation to vy
	vy[ss]+=dvy[k][i]*exp(-(y[j]*y[j])/(4*lambda*lambda));

	bx[ss]=b0*harris_field(y[j],para);
	by[ss]=0.0;
	bz[ss]=bg;
	cx[ss]=bx[ss];
	cy[ss]=by[ss];
	cz[ss]=bz[ss];
	
	pr[ss]=0.5*((1.0+beta)*(b0*b0+bg*bg)-(cx[ss]*cx[ss]+cy[ss]*cy[ss]+cz[ss]*cz[ss]));
	
	// Mag field perturbation by Zenitani
	bx[ss]-=b1amp*(y[j]/lambda)*exp(-(xm*xm+y[j]*y[j])/(4*lambda*lambda));
	by[ss]+=b1amp*(x[i]/lambda)*exp(-(x[i]*x[i]+ym*ym)/(4*lambda*lambda));
	cx[ss]-=b1amp*(y[j]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));
	cy[ss]+=b1amp*(x[i]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));
	
	cnsvt(ss);
      }
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);

  // Set kinematic viscosity and resistivity coefficients
  double al=sqrt((b0*b0+bg*bg)/ro0);
  nu0=al*lambda/REV;
  eta0=al*lambda/REM;
  setdc();
  // Message
  if (mpi_rank == 0){
    printf("Kinematic viscosity coef. = %f\n",nu0);
    printf("Resistivity coef. = %f\n",eta0);
  }
}

double DMHD3D::setdc()
{
  // Set dissipation coefficients and return their maximum.
  double dcmax=0.0,dtmp;
  for (int k=0;k<nz;k++){
    for (int j=0;j<ny;j++){
      for (int i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	
	// nu[ss]=nu0;
	// eta[ss]=eta0;

	// Localized dissipation
	double r=sqrt(x[i]*x[i]+y[j]*y[j]);
	nu[ss]=nu0*exp(-r/lambda);
	eta[ss]=eta0*exp(-r/lambda);
      
	dtmp=max(2*nu[ss],eta[ss]); // Factor 2 is multiplied in viscous coef. for robust estimation
	if (dtmp > dcmax) dcmax=dtmp;
      }
    }
  }

  // MPI Allreduce
  double dcmax_a;
  MPI_Allreduce(&dcmax,&dcmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  dcmax=dcmax_a;
  
  return dcmax;
}

inline double harris_field(double x, const double *params)
{
  /* Harris magnetic field */
  double x0=params[0],width=params[1]+1e-15;
  return( tanh((x-x0)/width) );
}
inline double harris_density(double x, const double *params)
{
  /* Harris distribution of the density */
  double x0=params[0],width=params[1]+1e-15;
  double cosh1=cosh((x-x0)/width);
  return( 1.0/(cosh1*cosh1) );
}
