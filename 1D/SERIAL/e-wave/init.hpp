void init_grid();
void init_plasma();

void init_grid()
// Define X coordinate
{
  int i;
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);
}

void init_plasma()
// Set initial condition
{
  int i;
  unsigned seed=10;
  // seed=(unsigned)time(NULL);
  double para1[]={0,0.01},para2[]={1,0.01};
  for (i=0;i<nx;i++){
    double vx,vy,vz,pr;
    ro[i]=ro0;
    vx=vx0;
    vy=0.0;
    vz=0.0;
    by[i]=rand_noise(para1,seed);
    bz[i]=rand_noise(para1,seed);
    pr=pr0*rand_noise(para2,seed);

    mhd_cnsvt(ro[i],vx,vy,vz,bx,by[i],bz[i],pr,&mx[i],&my[i],&mz[i],&en[i],gam);
  }
}
