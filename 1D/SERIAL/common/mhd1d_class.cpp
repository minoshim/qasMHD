#include "mhd1d_class.hpp"

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
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

std::string output_path(const std::string& directory, const std::string& filename)
{
  std::string path(directory);
  if (!path.empty() && path.back() != '/') path.push_back('/');
  return path+filename;
}

std::FILE *open_output(const std::string& path, const char *mode)
{
  std::FILE *outfil=std::fopen(path.c_str(),mode);
  if (outfil == nullptr) output_error("fopen",path);
  return outfil;
}

void close_output(std::FILE *outfil, const std::string& path)
{
  if (std::fclose(outfil) != 0) output_error("fclose",path);
}

}

void MHD1D::setup_grid(double xmin_value, double xmax_value)
{
  xmin=xmin_value;
  xmax=xmax_value;
  dx=(xmax-xmin)/XMESH;
  dt=cfl*dx;
  for (int i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
  }
}

void MHD1D::bound(double *val[], int nm, const int dnxs[])
{
  // Boundary condition
  for (int m=0;m<nm;m++){
    bc1d(val[m],nx,xoff,dnxs[m]);
  }
}

void MHD1D::setdt(int flg)
{
  // Set time step to meet CFL, if flg is set
  if (flg){
    double vtmp=0.0,vmax=1.0;
    for (int i=xoff;i<nx-xoff;i++){
      prmtv(i);
      vtmp=std::fabs(vx[i])+vfast(i);
      if (vtmp > vmax) vmax=vtmp;
    }
    dt=cfl*dx/vmax;
  }
  if (!dt_initialized){
    nrec=static_cast<int>(std::ceil(dtrec/dt)); // Step for output
    if (nrec < 1) nrec=1;
    dt=dtrec/static_cast<double>(nrec);
    nmax=nrec*nout;			// Maximum step
    dt_initialized=true;
    printf("Data output every %d steps (%f duration) \n",nrec,dtrec);
  }
}

void MHD1D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  if (n == 0) dout_();
  
  while( (flg)?(n++ < nmax && tim < tmax):(n++ < nmax) ){
    tim+=dt;

    bound(val,nm,dnxs);
    ideal(dt);
    setdt(flg*(n % 2 == 0));
    
    if ( (flg)?(tim >= trec):((n % nrec) == 0) ){
      cnt++;
      trec+=dtrec;
      dout_();
    }
  }
}

void MHD1D::dout_()
{
  // Output data
  if (n == 0){
    std::string path=output_path(fildir,"params.dat");
    std::FILE *outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%.12f\n",gam) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"xoff.dat");
    outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%d\n",xoff) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"t.dat");
    outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%.12f\n",tim) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"x.dat");
    outfil=open_output(path,"w");
    for (int i=0;i<nx;i++){
      if (std::fprintf(outfil,"%.12f\n",x[i]) < 0) output_error("fprintf",path);
    }
    close_output(outfil,path);
  } else{
    const std::string path=output_path(fildir,"t.dat");
    std::FILE *outfil=open_output(path,"a");
    if (std::fprintf(outfil,"%.12f\n",tim) < 0) output_error("fprintf",path);
    close_output(outfil,path);
  }

  std::ostringstream filename;
  filename << "outdat_" << std::setfill('0') << std::setw(5) << cnt << ".dat";
  const std::string path=output_path(fildir,filename.str());
  std::FILE *outfil=open_output(path,"wb");
  const std::size_t output_size=static_cast<std::size_t>(nx);
  for (int i=0;i<nm;i++){
    if (std::fwrite(val[i],sizeof(*val[i]),output_size,outfil) != output_size){
      output_error("fwrite",path);
    }
  }
  close_output(outfil,path);
}

MHD1D::MHD1D()
{
  // constructor  
  x=new double[nx];
  ro=new double[nx];
  mx=new double[nx];
  my=new double[nx];
  mz=new double[nx];
  bx=new double[nx];
  by=new double[nx];
  bz=new double[nx];
  en=new double[nx];
  vx=new double[nx];
  vy=new double[nx];
  vz=new double[nx];
  pr=new double[nx];
  // In 1D, cell center B is identical to cell edge B.
  cx=bx;
  cy=by;
  cz=bz;
  // Array of pointers for MHD variables.
  val[0]=ro;
  val[1]=mx;
  val[2]=my;
  val[3]=mz;
  val[4]=bx;
  val[5]=by;
  val[6]=bz;
  val[7]=en;
  // Staggered flag for ro,mx,my,mz,bx,by,bz,en (trivial in 1D).
  stxs[0]=0;
  stxs[1]=0;
  stxs[2]=0;
  stxs[3]=0;
  stxs[4]=1;
  stxs[5]=0;
  stxs[6]=0;
  stxs[7]=0;
}

MHD1D::~MHD1D()
{
  // destructor  
  delete[] x;
  delete[] ro;
  delete[] mx;
  delete[] my;
  delete[] mz;
  delete[] bx;
  delete[] by;
  delete[] bz;
  delete[] en;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] pr;
}
