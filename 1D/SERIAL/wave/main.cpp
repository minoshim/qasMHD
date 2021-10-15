#include "global.hpp"
using namespace global;

#include "new_delete.hpp"
#include "init.hpp"
#include "dataio.hpp"

int main(void)
{
  int i,n=0,cnt=0;
  double tim=0.0;

  new_delete();
  init_grid();
  init_plasma();

  // Re-calculate time step
  for (i=0;i<nx;i++){
    double vx,vy,vz,pr,cf,vtmp;
    mhd_prmtv(ro[i],mx[i],my[i],mz[i],bx,by[i],bz[i],en[i],&vx,&vy,&vz,&pr,gam);
    cf=sqrt((gam*pr+(bx*bx+by[i]*by[i]+bz[i]*bz[i]))/ro[i]);
    vtmp=fabs(vx)+cf;
    if (vtmp > vmax) vmax=vtmp;
  }
  dt=cfl*dx/vmax;
  nrec=(int)(dtrec/dt+0.5);	// Number of iterations for output
  nmax=nrec*nout;		// Number of maximum iteration
  printf("Data output every %d steps (%f duration) \n",nrec,dtrec);

  dataio(n,cnt,tim);
  
  /* Time integration */
  while(n++ < nmax){
    tim+=dt;

    // Boundary condition
    double *p[]={ro,mx,my,mz,by,bz,en};
    boundary(p,7,nx,xoff,dnxs);

    // MHD update
    mhd_fd1d(ro,mx,my,mz,by,bz,en,bx,dt,dx,nx,xoff,gam);

    /* Output */
    if ((n % nrec) == 0){
      cnt++;
      dataio(n,cnt,tim);
    }
  }
  
  new_delete();
  return 0;
}
