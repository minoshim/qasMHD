#include "dmhd3d_class.hpp"

double DMHD3D::setdc()
{
  // Set dissipation coefficients and return their maximum.
  double dcmax=0.0,dtmp;
  for (int k=0;k<nz;k++){
    for (int j=0;j<ny;j++){
      for (int i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	
	nu[ss]=nu0;
	eta[ss]=eta0;

	// Localized dissipation
	// double r=sqrt(x[i]*x[i]+y[j]*y[j]);
	// nu[ss]=nu0*exp(-r/1.0);
	// eta[ss]=eta0*exp(-r/1.0);
      
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

void DMHD3D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time
  // Dissipation solver called
  
  bkup_(0);			// Load backup data
  if (n == 0) dout_(0);
  
  double stim=MPI_Wtime();
  while(n++ < nmax && tim < tmax){
    tim+=dt;

    bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
    ideal(dt);

    double dcmax=setdc();
    if (dcmax > 1e-15){
      int nc,ncmax=1+(int)(2*(dcmax*dt)/(dr*dr));
      for (nc=0;nc<ncmax;nc++){	// sub-cycling
	bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
	dsptv(dt/ncmax);
      }
    }
    
    setdt(flg*(n % 2 == 0));
    
    if (tim >= trec){
      cnt++;
      trec+=dtrec;
      bkup_(1);			// Save backup data
      dout_(!mpi_rank);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  double etim=MPI_Wtime();
  if (!mpi_rank) printf("Elapse time = %lu sec.\n",(unsigned long)(etim-stim));
}

DMHD3D::DMHD3D(int* argc, char*** argv, int mnp) : MHD3D(argc,argv,mnp)
{
  // constructor  
  nu=new double[nd];
  eta=new double[nd];
}

DMHD3D::~DMHD3D()
{
  // destructor  
  delete[] nu;
  delete[] eta;
}
