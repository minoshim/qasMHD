#include "gmhd3d_class.hpp"
#include "mhd3d_io.hpp"

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace mpi3d_io;

GMHD3D::GMHD3D(int* argc, char*** argv, int mnp) : MHD3D(argc,argv,mnp)
{
  try{
    phi_g=new double[nd];
  } catch (...){
    abort_run("Cannot allocate gravitational potential.",MPI_ERR_NO_MEM);
  }
}

GMHD3D::~GMHD3D()
{
  delete[] phi_g;
}

void GMHD3D::setdt(int flg)
{
  // Set time step to meet CFL, if flg is set
  if (flg){
    double vtmp=0.0,vmax=1.0;
    for (int k=zoff;k<nz-zoff;k++){
      for (int j=yoff;j<ny-yoff;j++){
	for (int i=xoff;i<nx-xoff;i++){
	  int ss=nx*(ny*k+j)+i;
	  cx[ss]=0.5*(bx[ss]+bx[nx*(ny*k+j)+(i+stxs[4])]);
	  cy[ss]=0.5*(by[ss]+by[nx*(ny*k+(j+stys[5]))+i]);
	  cz[ss]=0.5*(bz[ss]+bz[nx*(ny*(k+stzs[6])+j)+i]);
	  GMHD3D::prmtv(ss);
	  vtmp=sqrt(vx[ss]*vx[ss]+vy[ss]*vy[ss]+vz[ss]*vz[ss])+vfast(ss);
	  if (!std::isfinite(vtmp)) abort_run("Nonfinite wave speed in setdt().");
	  if (vtmp > vmax) vmax=vtmp;
	}
      }
    }

    // MPI Allreduce
    double vmax_a;
    check_mpi(MPI_Allreduce(&vmax,&vmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD),
              "Wave-speed reduction failed.");
    vmax=vmax_a;

    dt=cfl*dr/vmax;
  }
  initialize_dt();
}

void GMHD3D::dout_(int msg)
try
{
  MHD3D::dout_(msg);
  if (n != 0) return;

  // Retain the rank-local, single-precision MPI potential format.
  std::ostringstream filename;
  filename << "g_potential_" << std::setfill('0') << std::setw(5) << mpi_rank << ".dat";
  const std::string path=output_path(fildir,filename.str());
  const std::size_t output_size=static_cast<std::size_t>(nd);
  std::vector<float> fval(output_size);
  conv_d2f(fval.data(),phi_g,nd);
  std::FILE *outfil=open_output(path,"wb");
  errno=0;
  if (std::fwrite(fval.data(),sizeof(float),output_size,outfil) != output_size){
    output_error("fwrite",path);
  }
  close_output(outfil,path);
} catch (...){
  abort_run("Exception while constructing output paths or allocating output buffers.");
}
