void init_grid();
void init_plasma();

void init_grid()
// Define X and Y coordinates
{
  int i,j;
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (j=0;j<ny;j++){
    y[j]=(j-yoff+0.5)*dy+ymin;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
// Set initial condition
{
  // RMI
  int i,j,ss;
  double dw=0.5*dy;		// Shock width
  double gp1=gam+1.0;
  double gm1=gam-1.0;

#if (MAGNET)
  // MHD shock
  // R-H relation
  double a=2.0*(2.0-gam)/beta;
  double b=gam*(gm1*ma_u*ma_u+2.0/beta+2.0);
  double c=-gam*gp1*ma_u*ma_u;
  // Compression ratio
  double r=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
  // Pr2/Pr1
  double rr=gam*ma_u*ma_u*(1.0-1.0/r)-(r*r-1.0)/beta+1.0;
#else
  // HD shock
  double r=(gp1*ma_u*ma_u)/(gm1*ma_u*ma_u+2.0); // Compression ration
  double rr=(2.0*gam*ma_u*ma_u-gm1)/gp1;	// Pr2/Pr1
#endif
  // Donwstream paramters
  double ro_2=ro_1*r;
  double vy_2=vy_1/r;
  double bx_2=bx_1*r;
  double bz_2=bz_1*r;
  double pr_2=pr_1*rr;

  // Messaage of initial condition
  {
    puts("Initial condition in shock rest frame");
    printf("U: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_1,vx_1,vy_1,vz_1,bx_1,by_1,bz_1,pr_1);
    printf("D: %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",ro_2,vx_1,vy_2,vz_1,bx_2,by_1,bz_2,pr_2);
    printf("Mach = %.9f\n",ma_u);
  }

  // Random seed for contact discon.
  unsigned seed;
  seed=(unsigned)time(NULL);

  // Position of contact discon.
  double ypos[nx];
  for (i=0;i<nx;i++){
    ypos[i]=1.0+psi*cos(2*pi*x[i]/lambda); // Single mode
  }
#if (RANDOM)
  // Multi-mode
  double para[2]={pi,pi};
  for (int m=2;m<=mmax;m++){
    double xphase;
    xphase=rand_noise(para,seed);
    for (i=0;i<nx;i++){
      ypos[i]+=(psi/m)*cos(2*pi*m*(x[i]-xphase)/lambda);
    }
  }
#endif

  for (j=0;j<ny;j++){
    double sfunc=0.5*(1.0+tanh(y[j]/dw)); // 0 (y<0), 1 (y>0)
    for (i=0;i<nx;i++){
      double vx,vy,vz,cx,cy,cz,pr;

      ss=nx*j+i;

      ro[ss]=ro_2+(ro_1-ro_2)*sfunc;

      vx=vx_1;
      vy=vy_2+(vy_1-vy_2)*sfunc;
      vy-=vref;
      vz=vz_1;

      bx[ss]=cx=bx_2+(bx_1-bx_2)*sfunc;
      by[ss]=cy=by_1;
      bz[ss]=cz=bz_2+(bz_1-bz_2)*sfunc;
      
      pr=pr_2+(pr_1-pr_2)*sfunc;
      
      // Contact discon
      double rocd=ro_3;
#if (RANDOM)
      double rpara[2]={ro_3,dro3};
      rocd=rand_noise(rpara,seed);
#endif
      ro[ss]+=(rocd-ro_1)*0.5*(1.0+tanh((y[j]-ypos[i])/dw));

      mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);      
    }
  }
}
