void new_delete()
{
  static int nflg=0;
  if (nflg == 0){
    nflg=1;
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
    phi_g=new double[nd];
  } else{
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
    delete[] phi_g;
    nflg=0;
  }
}
