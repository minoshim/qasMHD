#include "gmhd2d_class.hpp"
#include "mhd2d_io.hpp"

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>

using namespace serial2d_io;

GMHD2D::GMHD2D()
{
  phi_g=new double[nd];
}

GMHD2D::~GMHD2D()
{
  delete[] phi_g;
}

void GMHD2D::setdt(int flg)
{
  if (flg){
    double vtmp=0.0,vmax=1.0;
    for (int j=yoff;j<ny-yoff;j++){
      for (int i=xoff;i<nx-xoff;i++){
        const int ss=nx*j+i;
        cx[ss]=0.5*(bx[ss]+bx[nx*j+(i+stxs[4])]);
        cy[ss]=0.5*(by[ss]+by[nx*(j+stys[5])+i]);
	GMHD2D::prmtv(ss);
        vtmp=sqrt(vx[ss]*vx[ss]+vy[ss]*vy[ss])+vfast(ss);
        if (vtmp > vmax) vmax=vtmp;
      }
    }
    dt=cfl*dr/vmax;
  }
  // Reuse output-interval initialization without the ordinary MHD CFL calculation.
  MHD2D::setdt(0);
}

void GMHD2D::dout_(int msg)
{
  MHD2D::dout_(msg);
  if (n != 0) return;

  const std::string path=output_path(fildir,"g_potential.dat");
  std::FILE *outfil=open_output(path,"wb");
  const std::size_t output_size=static_cast<std::size_t>(nd);
  errno=0;
  if (std::fwrite(phi_g,sizeof(*phi_g),output_size,outfil) != output_size){
    output_error("fwrite",path);
  }
  close_output(outfil,path);
}
