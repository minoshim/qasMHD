#include "dmhd2d_class.hpp"
#include "mhd2d_io.hpp"

#include <ctime>
#include <random>
#include <vector>

inline double harris_field(double x, const double *params);
inline double harris_density(double x, const double *params);

const double lambda=1.0;	// Current sheet thickness

void DMHD2D::init_()
try
{
  // Magnetic reconnection
  int i,j;
  // Initial condition parameters
  const double beta=0.2;	// Plasma beta @ lobe
  const double ro0=1.0;		// Density @ CS
  const double ro1=0.2;		// Density @ lobe
  const double b0=1.0;		// Mag field @ lobe
  const double b1=0.05;		// Mag field perturbation by Zenitani
  const double bg=0.0;		// Guide mag field along Z
  const double dv=0.01;		// Random noize perturbation to Vy (avaiable when RANDOM=1)
  const double para[2]={0,lambda};

  std::vector<double> dvy;
  resize_work(dvy,nx); // Allocation failure must stop all ranks.
#if (RANDOM)
  // For reproducible runs, replace std::time(nullptr) with a fixed seed.
  unsigned seed=0;
  if (mpi_rank == 0) seed=static_cast<unsigned>(std::time(nullptr));
  mpi2d_io::check_mpi(MPI_Bcast(&seed,1,MPI_UNSIGNED,0,MPI_COMM_WORLD),
                      "Random seed broadcast failed.");
  // The same X subdomain uses the same stream, independent of Y rank.
  std::seed_seq seeds{seed,static_cast<unsigned>(mpi_ranx)};
  std::mt19937 engine(seeds);
#endif

  for (i=0;i<nx;i++){
    dvy[i]=0.0;
#if (RANDOM)
    double dvpara[2]={0,dv};
    dvy[i]+=rand_noise_mt(dvpara,engine); // Multiple mode perturbation
#endif    
  }
  
  for (j=0;j<ny;j++){
    double ym=y[j]-0.5*dy;
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double xm=x[i]-0.5*dx;
      
      // ro[ss]=(ro0-ro1)*harris_density(y[j],para)+ro1;
      ro[ss]=ro0*harris_density(y[j],para)+ro1;

      vx[ss]=0.0;
      vy[ss]=0.0;
      vz[ss]=0.0;
      // Perturbation to vy
      vy[ss]+=dvy[i]*exp(-(y[j]*y[j])/(4*lambda*lambda));

      bx[ss]=b0*harris_field(y[j],para);
      by[ss]=0.0;
      bz[ss]=bg;
      cx[ss]=bx[ss];
      cy[ss]=by[ss];
      // In 2D, cell center Bz is identical to cell edge Bz.
      // Thus cz is not explicitly initialized here.
      // They should be defined in constructer (see mhd2d_class.cpp).

      pr[ss]=0.5*((1.0+beta)*(b0*b0+bg*bg)-(cx[ss]*cx[ss]+cy[ss]*cy[ss]+cz[ss]*cz[ss]));

      // Mag field perturbation by Zenitani
      bx[ss]-=b1*(y[j]/lambda)*exp(-(xm*xm+y[j]*y[j])/(4*lambda*lambda));
      by[ss]+=b1*(x[i]/lambda)*exp(-(x[i]*x[i]+ym*ym)/(4*lambda*lambda));
      cx[ss]-=b1*(y[j]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));
      cy[ss]+=b1*(x[i]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));

      cnsvt(ss);
    }
  }
  
  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys);

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
} catch (...){
  mpi2d_io::abort_run("Exception while initializing MRX state or random generator.");
}

double DMHD2D::setdc()
{
  // Set dissipation coefficients and return their maximum.
  double dcmax=0.0,dtmp;
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      int ss=nx*j+i;
      
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

  // MPI Allreduce
  double dcmax_a;
  mpi2d_io::check_mpi(MPI_Allreduce(&dcmax,&dcmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD),
                      "Dissipation coefficient reduction failed.");
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
