#include "mhd1d_class.hpp"

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
  static int nflg=0;
  if (flg){
    double vtmp=0.0,vmax=1.0;
    for (int i=xoff;i<nx-xoff;i++){
      prmtv(i);
      vtmp=fabs(vx[i])+vfast(i);
      if (vtmp > vmax) vmax=vtmp;
    }
    dt=cfl*dx/vmax;
  }
  if (nflg == 0){
    nrec=(int)(dtrec/dt+0.5);		// Step for output
    nmax=nrec*nout;			// Maximum step
    nflg=1;
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
  int i;
  char filname[200];
  FILE *outfil;
  if (n == 0){
    sprintf(filname,"%s/params.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%.12f\n",gam);
    fclose(outfil);

    sprintf(filname,"%s/xoff.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d\n",xoff);
    fclose(outfil);

    sprintf(filname,"%s/t.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%.12f\n",tim);
    fclose(outfil);

    sprintf(filname,"%s/x.dat",fildir);  
    outfil=fopen(filname,"w");
    for (i=0;i<nx;i++){
      fprintf(outfil,"%.12f\n",x[i]);
    }
    fclose(outfil);
  } else{
    sprintf(filname,"%s/t.dat",fildir);  
    outfil=fopen(filname,"a");
    fprintf(outfil,"%.12f\n",tim);
    fclose(outfil);
  }
  sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
  outfil=fopen(filname,"wb");
  for (i=0;i<sizeof(val)/sizeof(val[0]);i++){
    fwrite(val[i],sizeof(*val[i]),nx,outfil);
  }
  fclose(outfil);
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
