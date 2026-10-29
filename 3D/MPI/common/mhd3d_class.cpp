#include "mhd3d_class.hpp"
#include "mhd3d_io.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

using namespace mpi3d_io;

void MHD3D::setup_grid(double xmin_value, double xmax_value,
                       double ymin_value, double ymax_value,
                       double zmin_value, double zmax_value,
                       double xshift, double yshift, double zshift)
{
  xmin=xmin_value;
  xmax=xmax_value;
  ymin=ymin_value;
  ymax=ymax_value;
  zmin=zmin_value;
  zmax=zmax_value;
  dx=(xmax-xmin)/XMESH;
  dy=(ymax-ymin)/YMESH;
  dz=(zmax-zmin)/ZMESH;
  dr=min(min(dx,dy),dz);
  dt=cfl*dr; // Recalculated by setdt() after initialization.
  int isum=0,jsum=0,ksum=0;
  for (int m=0;m<mpi_ranx;m++) isum+=(XMESH+m)/mpi_numx;
  for (int m=0;m<mpi_rany;m++) jsum+=(YMESH+m)/mpi_numy;
  for (int m=0;m<mpi_ranz;m++) ksum+=(ZMESH+m)/mpi_numz;
  for (int i=0;i<nx;i++) x[i]=(i-xoff+isum+xshift)*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-yoff+jsum+yshift)*dy+ymin;
  for (int k=0;k<nz;k++) z[k]=(k-zoff+ksum+zshift)*dz+zmin;
}

void MHD3D::bound(double *val[], int nm, const int stxs[], const int dnxs[], const int stys[], const int dnys[], const int stzs[], const int dnzs[])
{
  // Boundary condition
  mpi_sdrv3d(val,nm,nx,ny,nz,xoff,yoff,zoff,dnxs[0],dnys[0],dnzs[0],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  /* Note: If dnx(y,z)s[0] == 0, boundary is periodic*/
  for (int m=0;m<nm;m++){
    mpi_xbc3d(val[m],nx,ny,nz,xoff,yoff,zoff,stxs[m],dnxs[m],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_ybc3d(val[m],nx,ny,nz,xoff,yoff,zoff,stys[m],dnys[m],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
    mpi_zbc3d(val[m],nx,ny,nz,xoff,yoff,zoff,stzs[m],dnzs[m],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  }
}

void MHD3D::setdt(int flg)
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
	  prmtv(ss);
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

void MHD3D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  prepare_run(flg); // Load backup data without reinitializing its dt
  if (n == 0) dout_(0);
  
  double stim=MPI_Wtime();
  while(flg ? (tim < tmax) : (n < nmax)){
    check_step();
    ++n;
    const double remaining=tmax-tim;
    const bool last_step=flg && dt >= remaining;
    const double dt_step=last_step ? remaining : dt;
    if (flg){
      tim=last_step ? tmax : tim+dt_step;
    } else{
      tim=(n == nmax) ? tmax : n*dt;
    }

    bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
    ideal(dt_step);

    setdt(flg*(n % 2 == 0));
    
    if (flg ? (tim >= trec) : ((n % nrec) == 0)){
      cnt++;
      trec=(cnt+1)*dtrec;
      dout_(!mpi_rank);
      // Do not advance the checkpoint until every rank has finished its output.
      check_mpi(MPI_Barrier(MPI_COMM_WORLD),"Output barrier failed.");
      bkup_(1); // Save backup data; wait for all ranks to succeed
    }
  }
  double etim=MPI_Wtime();
  if (!mpi_rank) printf("Elapse time = %lu sec.\n",(unsigned long)(etim-stim));
}

void MHD3D::dout_(int msg)
try
{
  // Retain the MPI single-precision output format used by merge.out.
  const std::size_t output_size=static_cast<std::size_t>(nd);
  std::vector<float> fval(output_size);
  if (n == 0){
    if (mpi_rank == 0){
      std::string path=output_path(fildir,"params.dat");
      std::FILE *outfil=open_output(path,"w");
      if (std::fprintf(outfil,"%f\n",gam) < 0) output_error("fprintf",path);
      close_output(outfil,path);

      path=output_path(fildir,"offsets.dat");
      outfil=open_output(path,"w");
      if (std::fprintf(outfil,"%d %d %d\n",xoff,yoff,zoff) < 0) output_error("fprintf",path);
      close_output(outfil,path);

      path=output_path(fildir,"mpinum.dat");
      outfil=open_output(path,"w");
      if (std::fprintf(outfil,"%d %d %d\n",mpi_numx,mpi_numy,mpi_numz) < 0) output_error("fprintf",path);
      close_output(outfil,path);

      path=output_path(fildir,"t.dat");
      outfil=open_output(path,"w");
      if (std::fprintf(outfil,"%f\n",tim) < 0) output_error("fprintf",path);
      close_output(outfil,path);
    }
    if ((mpi_ranz == 0) && (mpi_rany == 0)){
      std::ostringstream filename;
      filename << "x_" << std::setfill('0') << std::setw(5) << mpi_ranx << ".dat";
      const std::string path=output_path(fildir,filename.str());
      std::FILE *outfil=open_output(path,"w");
      for (int i=0;i<nx;i++){
        if (std::fprintf(outfil,"%.12f\n",x[i]) < 0) output_error("fprintf",path);
      }
      close_output(outfil,path);
    }
    if ((mpi_ranz == 0) && (mpi_ranx == 0)){
      std::ostringstream filename;
      filename << "y_" << std::setfill('0') << std::setw(5) << mpi_rany << ".dat";
      const std::string path=output_path(fildir,filename.str());
      std::FILE *outfil=open_output(path,"w");
      for (int i=0;i<ny;i++){
        if (std::fprintf(outfil,"%.12f\n",y[i]) < 0) output_error("fprintf",path);
      }
      close_output(outfil,path);
    }
    if ((mpi_rank%m_xy) == 0){
      std::ostringstream filename;
      filename << "z_" << std::setfill('0') << std::setw(5) << mpi_ranz << ".dat";
      const std::string path=output_path(fildir,filename.str());
      std::FILE *outfil=open_output(path,"w");
      for (int i=0;i<nz;i++){
        if (std::fprintf(outfil,"%.12f\n",z[i]) < 0) output_error("fprintf",path);
      }
      close_output(outfil,path);
    }
  } else{
    if (mpi_rank == 0){
      const std::string path=output_path(fildir,"t.dat");
      std::FILE *outfil=open_output(path,"a");
      if (std::fprintf(outfil,"%f\n",tim) < 0) output_error("fprintf",path);
      close_output(outfil,path);
    }
  }
  std::ostringstream filename;
  filename << "outdat_" << std::setfill('0') << std::setw(5) << cnt
           << "_" << std::setw(5) << mpi_rank << ".dat";
  const std::string path=output_path(fildir,filename.str());
  std::FILE *outfil=open_output(path,"wb");
  for (int i=0;i<nm;i++){
    conv_d2f(fval.data(),val[i],nd);
    if (std::fwrite(fval.data(),sizeof(float),output_size,outfil) != output_size){
      output_error("fwrite",path);
    }
  }
  close_output(outfil,path);
  if (msg){
    printf("Output data at t=%.4f (%d iterations).\n",tim,n);
  }
} catch (...){
  abort_run("Exception while constructing output paths or allocating output buffers.");
}

MHD3D::MHD3D(int* argc, char*** argv, int mnp) : MYMPI(argc,argv,mnp)
{
  // constructor  
  x=new double[nx];
  y=new double[ny];
  z=new double[nz];
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
  cz=new double[nd];
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

  stzs[0]=0;
  stzs[1]=0;
  stzs[2]=0;
  stzs[3]=0;
  stzs[4]=0;
  stzs[5]=0;
  stzs[6]=1;
  stzs[7]=0;
}

MHD3D::~MHD3D()
{
  // destructor  
  delete[] x;
  delete[] y;
  delete[] z;
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
  delete[] cz;
}
