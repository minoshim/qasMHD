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

  // Set normal magnetic field
#if (NUM == 2)
  // BrioWu
  bx=0.75;
  gam=2.0;
#elif (NUM == 3)
  // Slow switch-off shock
  bx=1.0;
#elif (NUM == 4)
  // Slow switch-off rarefaction
  bx=1.0;
#elif (NUM == 5)
  // Super-fast expansion
  bx=0.0;
#else  // Default
  // DaiWoodward
  bx=2.0/sqrt(4*pi);
#endif

  for (i=0;i<nx;i++){
    double vx,vy,vz,pr;
    char flag=(x[i] <= 0);

#if (NUM == 2)
    // BrioWu
    ro[i]=(flag)?(1.00):(0.125);
    vx=0;
    vy=0;
    vz=0;
    by[i]=(flag)?(1):(-1);
    bz[i]=0;
    pr=(flag)?(1.0):(0.1);
#elif (NUM ==3)
  // Slow switch-off shock
    ro[i]=(flag)?(1.368):(1.0);
    vx=(flag)?(0.269):(0.0);
    vy=(flag)?(1.0):(0.0);
    vz=(flag)?(0.0):(0.0);
    by[i]=(flag)?(0.0):(1.0);
    bz[i]=0.0;
    pr=(flag)?(1.769):(1.0);
#elif (NUM == 4)
  // Slow switch-off rarefaction
    ro[i]=(flag)?(1.0):(0.2);
    vx=(flag)?(0.0):(1.186);
    vy=(flag)?(0.0):(2.967);
    vz=(flag)?(0.0):(0.0);
    by[i]=(flag)?(0.0):(1.6405);
    bz[i]=0.0;
    pr=(flag)?(2.0):(0.1368);
#elif (NUM == 5)    
  // Super-fast expansion
    ro[i]=1.0;
    vx=(flag)?(-3.1):(3.1);
    vy=0.0;
    vz=0.0;
    by[i]=0.5;
    bz[i]=0.0;
    pr=0.45;
#else  // Default
    // DaiWoodward
    ro[i]=(flag)?(1.08):(1.0);
    vx=(flag)?(1.20):(0.0);
    vy=(flag)?(0.01):(0.0);
    vz=(flag)?(0.50):(0.0);
    by[i]=(flag)?(3.6/sqrt(4*pi)):(4.0/sqrt(4*pi));
    bz[i]=2.0/sqrt(4*pi);
    pr=(flag)?(0.95):(1.0);
#endif

    mhd_cnsvt(ro[i],vx,vy,vz,bx,by[i],bz[i],pr,&mx[i],&my[i],&mz[i],&en[i],gam);
  }
    
}
