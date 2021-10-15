#include "global.hpp"
using namespace global;

#include "new_delete.hpp"
#include "init.hpp"
#include "get.hpp"
#include "dataio.hpp"

int main(void)
{
  int n=0,cnt=0;
  double tim=0,trec=dtrec;
  clock_t stim,etim;

  new_delete();
  init_grid();
  init_plasma();
  double *p[]={ro,mx,my,mz,bx,by,bz,en};
  
  get_dt(&dt,cfl,dr,1);
  
  dataio(p,8,nd,n,cnt,tim,0);

  stim=clock();
  // Time integration
  while(n++ < nmax && tim < tend){
    tim+=dt;

    // Boundary condition
    boundary(p,8,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys);
    
    // MHD update
    mhd_fd2d(p,dt,dx,dy,8,nx,ny,xoff,yoff,gam);

    // Re-calculate time step
    if (n % 2 == 0) get_dt(&dt,cfl,dr,0);

    // Output
    if (tim >= trec){
      trec+=dtrec;
      cnt++;
      dataio(p,8,nd,n,cnt,tim,1);
    }
  }
  etim=clock();
  printf("CPU time = %f sec.\n",(double)(etim-stim)/CLOCKS_PER_SEC);

  new_delete();
  return 0;
}
