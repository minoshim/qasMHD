#include "gmhd2d_class.hpp"

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {

[[noreturn]] void output_error(const char *operation, const std::string& path)
{
  const int error_number=errno;
  const char *reason=(error_number != 0)?std::strerror(error_number):"unknown I/O error";
  std::fprintf(stderr,"Output error: %s failed for '%s': %s\n",
               operation,path.c_str(),reason);
  std::exit(EXIT_FAILURE);
}

}

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

  std::string path=fildir;
  if (!path.empty() && path.back() != '/') path.push_back('/');
  path+="g_potential.dat";
  errno=0;
  std::FILE *outfil=std::fopen(path.c_str(),"wb");
  if (outfil == nullptr) output_error("fopen",path);
  const std::size_t output_size=static_cast<std::size_t>(nd);
  errno=0;
  if (std::fwrite(phi_g,sizeof(*phi_g),output_size,outfil) != output_size){
    output_error("fwrite",path);
  }
  errno=0;
  if (std::fclose(outfil) != 0) output_error("fclose",path);
}
