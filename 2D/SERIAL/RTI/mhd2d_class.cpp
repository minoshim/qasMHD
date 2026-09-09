#include "mhd2d_class.hpp"

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
  errno=0;
  std::FILE *outfil=std::fopen(path.c_str(),mode);
  if (outfil == nullptr) output_error("fopen",path);
  return outfil;
}

void close_output(std::FILE *outfil, const std::string& path)
{
  errno=0;
  if (std::fclose(outfil) != 0) output_error("fclose",path);
}

}

void MHD2D::bound(double *val[], int nm, const int stxs[], const int dnxs[], const int stys[], const int dnys[])
{
  // Boundary condition
  for (int m=0;m<nm;m++){
    bc2d(val[m],nx,ny,xoff,yoff,stxs[m],dnxs[m],stys[m],dnys[m]);
  }

  // With finite gravity along Y, Y boundary for energy needs special care
  if (nm == 8){
    int i,j,ss,sb;
    double fac=(2.0-gam)/(gam-1.0);
    {
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sb=nx*yoff+i;
	  en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	}
      }
    }
    {
      for (j=ny-yoff;j<ny;j++){
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sb=nx*(ny-yoff-1)+i;
	  en[ss]=en[sb]-fac*ro[sb]*(phi_g[ss]-phi_g[sb]);
	}
      }
    }
  }
}

void MHD2D::setdt(int flg)
{
  // Set time step to meet CFL, if flg is set
  if (flg){
    double vtmp=0.0,vmax=1.0;
    for (int j=yoff;j<ny-yoff;j++){
      for (int i=xoff;i<nx-xoff;i++){
	int ss=nx*j+i;
	cx[ss]=0.5*(bx[ss]+bx[nx*j+(i+stxs[4])]);
	cy[ss]=0.5*(by[ss]+by[nx*(j+stys[5])+i]);
	en[ss]-=ro[ss]*phi_g[ss]; // Subtract G-potential before calling prmtv()
	prmtv(ss);
	en[ss]+=ro[ss]*phi_g[ss]; // Return G-potential after calling prmtv()
	vtmp=sqrt(vx[ss]*vx[ss]+vy[ss]*vy[ss])+vfast(ss);
	if (vtmp > vmax) vmax=vtmp;
      }
    }
    dt=cfl*dr/vmax;
  }
  if (!dt_initialized){
    nrec=static_cast<int>(std::ceil(dtrec/dt));
    if (nrec < 1) nrec=1;
    dt=dtrec/static_cast<double>(nrec);
    nmax=nrec*nout;
    dt_initialized=true;
  }
}

void MHD2D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  if (n == 0) dout_(0);
  
  clock_t stim=clock();
  while(flg ? (tim < tmax) : (n < nmax)){
    ++n;
    const double remaining=tmax-tim;
    const bool last_step=flg && dt >= remaining;
    const double dt_step=last_step ? remaining : dt;
    if (flg){
      tim=last_step ? tmax : tim+dt_step;
    } else{
      tim=(n == nmax) ? tmax : n*dt;
    }

    bound(val,nm,stxs,dnxs,stys,dnys);
    ideal(dt_step);

    setdt(flg*(n % 2 == 0));
    
    if (flg ? (tim >= trec) : ((n % nrec) == 0)){
      cnt++;
      trec=(cnt+1)*dtrec;
      dout_(1);
    }
  }
  clock_t etim=clock();
  printf("CPU time = %f sec.\n",(double)(etim-stim)/CLOCKS_PER_SEC);

}

void MHD2D::dout_(int msg)
{
  // Output data
  const std::size_t output_size=static_cast<std::size_t>(nd);
  if (n == 0){
    std::string path=output_path(fildir,"params.dat");
    std::FILE *outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%f\n",gam) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"offsets.dat");
    outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%d %d\n",xoff,yoff) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"t.dat");
    outfil=open_output(path,"w");
    if (std::fprintf(outfil,"%f\n",tim) < 0) output_error("fprintf",path);
    close_output(outfil,path);

    path=output_path(fildir,"x.dat");
    outfil=open_output(path,"w");
    for (int i=0;i<nx;i++){
      if (std::fprintf(outfil,"%.12f\n",x[i]) < 0) output_error("fprintf",path);
    }
    close_output(outfil,path);

    path=output_path(fildir,"y.dat");
    outfil=open_output(path,"w");
    for (int i=0;i<ny;i++){
      if (std::fprintf(outfil,"%.12f\n",y[i]) < 0) output_error("fprintf",path);
    }
    close_output(outfil,path);

    path=output_path(fildir,"g_potential.dat");
    outfil=open_output(path,"wb");
    if (std::fwrite(phi_g,sizeof(*phi_g),output_size,outfil) != output_size){
      output_error("fwrite",path);
    }
    close_output(outfil,path);
  } else{
    const std::string path=output_path(fildir,"t.dat");
    std::FILE *outfil=open_output(path,"a");
    if (std::fprintf(outfil,"%f\n",tim) < 0) output_error("fprintf",path);
    close_output(outfil,path);
  }

  std::ostringstream filename;
  filename << "outdat_" << std::setfill('0') << std::setw(5) << cnt << ".dat";
  const std::string path=output_path(fildir,filename.str());
  std::FILE *outfil=open_output(path,"wb");
  for (int i=0;i<nm;i++){
    if (std::fwrite(val[i],sizeof(*val[i]),output_size,outfil) != output_size){
      output_error("fwrite",path);
    }
  }
  close_output(outfil,path);
  if (msg){
    printf("Output data at t=%.4f (%d iterations).\n",tim,n);
  }
}

MHD2D::MHD2D()
{
  // constructor  
  x=new double[nx];
  y=new double[ny];
  ro=new double[nd];
  mx=new double[nd];
  my=new double[nd];
  mz=new double[nd];
  bx=new double[nd];
  by=new double[nd];
  bz=new double[nd];
  en=new double[nd];
  vx=new double[nd];
  vy=new double[nd];
  vz=new double[nd];
  pr=new double[nd];
  cx=new double[nd];
  cy=new double[nd];
  // In 2D, cell center Bz is identical to cell edge Bz.
  cz=bz;
  // Gravitational potential
  phi_g=new double[nd];
  // Array of pointers for MHD variables.
  val[0]=ro;
  val[1]=mx;
  val[2]=my;
  val[3]=mz;
  val[4]=bx;
  val[5]=by;
  val[6]=bz;
  val[7]=en;
  // Staggered flag for ro,mx,my,mz,bx,by,bz,en.
  stxs[0]=0;
  stxs[1]=0;
  stxs[2]=0;
  stxs[3]=0;
  stxs[4]=1;
  stxs[5]=0;
  stxs[6]=0;
  stxs[7]=0;

  stys[0]=0;
  stys[1]=0;
  stys[2]=0;
  stys[3]=0;
  stys[4]=0;
  stys[5]=1;
  stys[6]=0;
  stys[7]=0;
}

MHD2D::~MHD2D()
{
  // destructor  
  delete[] x;
  delete[] y;
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
  delete[] cx;
  delete[] cy;
  delete[] phi_g;
}
