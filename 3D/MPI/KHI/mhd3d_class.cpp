#include "mhd3d_class.hpp"

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
  static int nflg=0;
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
	  if (vtmp > vmax) vmax=vtmp;
	}
      }
    }

    // MPI Allreduce
    double vmax_a;
    MPI_Allreduce(&vmax,&vmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    vmax=vmax_a;

    dt=cfl*dr/vmax;
  }
  if (nflg == 0){
    nmax=(int)(tmax/dt+0.5);
    nflg=1;
  }
}

void MHD3D::exec_(int flg)
{
  // Run simulation.
  // If flg=0, dt unchanged and output @ constant step
  // If flg=1, dt changed and output @ constant time

  bkup_(0);			// Load backup data
  if (n == 0) dout_(0);
  
  double stim=MPI_Wtime();
  while(n++ < nmax && tim < tmax){
    tim+=dt;

    bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
    ideal(dt);

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

void MHD3D::dout_(int msg)
{
  // Output data
  int i;
  char filname[200];
  FILE *outfil;
  float *fval;
  fval=new float[nd];
  if (n == 0){
    if (mpi_rank == 0){
      sprintf(filname,"%s/params.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f\n",gam);
      fclose(outfil);

      sprintf(filname,"%s/offsets.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%d %d %d\n",xoff,yoff,zoff);
      fclose(outfil);

      sprintf(filname,"%s/mpinum.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%d %d %d\n",mpi_numx,mpi_numy,mpi_numz);
      fclose(outfil);
      
      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
    if ((mpi_ranz == 0) && (mpi_rany == 0)){
      sprintf(filname,"%s/x_%05d.dat",fildir,mpi_ranx);
      outfil=fopen(filname,"w");
      for (i=0;i<nx;i++) fprintf(outfil,"%.12f\n",x[i]);
      fclose(outfil);
    }
    if ((mpi_ranz == 0) && (mpi_ranx == 0)){
      sprintf(filname,"%s/y_%05d.dat",fildir,mpi_rany);
      outfil=fopen(filname,"w");
      for (i=0;i<ny;i++) fprintf(outfil,"%.12f\n",y[i]);
      fclose(outfil);
    }
    if ((mpi_rank%m_xy) == 0){
      sprintf(filname,"%s/z_%05d.dat",fildir,mpi_ranz);
      outfil=fopen(filname,"w");
      for (i=0;i<nz;i++) fprintf(outfil,"%.12f\n",z[i]);
      fclose(outfil);
    }
  } else{
    if (mpi_rank == 0){
      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"a");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
  }
  sprintf(filname,"%s/outdat_%05d_%05d.dat",fildir,cnt,mpi_rank);
  outfil=fopen(filname,"wb");
  for (i=0;i<nm;i++){
    conv_d2f(fval,val[i],nd);
    fwrite(fval,sizeof(*fval),nd,outfil);
  }
  fclose(outfil);
  if (msg){
    printf("Output data at t=%.4f (%d iterations).\n",tim,n);
  }

  delete[] fval;
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
