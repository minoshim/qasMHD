void new_delete()
{
  static int nflg=0;
  if (nflg == 0){
    nflg=1;
    x=new double[nx];
    ro=new double[nx];
    mx=new double[nx];
    my=new double[nx];
    mz=new double[nx];
    by=new double[nx];
    bz=new double[nx];
    en=new double[nx];
  } else{
    delete[] x;
    delete[] ro;
    delete[] mx;
    delete[] my;
    delete[] mz;
    delete[] by;
    delete[] bz;
    delete[] en;
    nflg=0;
  }
}
