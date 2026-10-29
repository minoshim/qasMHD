#include "mhd2d_class.hpp"
#include "mhd2d_io.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace mpi2d_io;

namespace {

int checked_nrec(double steps, int nout)
{
  const int max_int=std::numeric_limits<int>::max();
  if (!std::isfinite(steps) || steps < 1.0 || steps > max_int){
    abort_run("Steps per output exceed the supported integer range.");
  }
  const int nrec=static_cast<int>(steps);
  if (nout > max_int/nrec) abort_run("Total step count exceeds the supported integer range.");
  return nrec;
}

bool same_time(double a, double b)
{
  const double scale=std::max(1.0,std::max(std::fabs(a),std::fabs(b)));
  return std::fabs(a-b) <= 128.0*std::numeric_limits<double>::epsilon()*scale;
}

}

void MHD2D::check_time_parameters() const
{
  if (!std::isfinite(dtrec) || dtrec <= 0.0 || nout < 0 ||
      nout == std::numeric_limits<int>::max() || !std::isfinite(tmax) ||
      !std::isfinite(cfl) || cfl <= 0.0 || !std::isfinite(dr) || dr <= 0.0 ||
      !std::isfinite(dt) || dt <= 0.0){
    abort_run("Invalid output interval, CFL, grid spacing, or time step.");
  }
}

void MHD2D::initialize_dt()
{
  check_time_parameters();
  if (!dt_initialized){
    // Same initial output-interval alignment as 2D/SERIAL.
    const double steps=std::max(1.0,std::ceil(dtrec/dt));
    nrec=checked_nrec(steps,nout);
    dt=dtrec/static_cast<double>(nrec);
    nmax=nrec*nout;
    dt_initialized=true;
    check_time_parameters();
  }
}

void MHD2D::check_step() const
{
  check_time_parameters();
  if (n == std::numeric_limits<int>::max() ||
      !std::isfinite(tim) || tim < 0.0 ||
      tim+std::min(dt,tmax-tim) <= tim){
    abort_run("Time integration cannot advance (step counter or floating-point limit).");
  }
}

void MHD2D::check_restart_state() const
{
  check_time_parameters();
  if (n < 0 || cnt < 0 || cnt > nout || !std::isfinite(tim) ||
      tim < 0.0 || tim > tmax || !std::isfinite(trec) || trec <= 0.0 ||
      !same_time(trec,(cnt+1)*dtrec) || ((n == 0) != (tim == 0.0))){
    abort_run("Invalid or incompatible backup counters/times.");
  }

  // Every rank must resume the same step and follow the same collective sequence.
  const int counts[]={n,cnt};
  int counts_min[2],counts_max[2];
  const double times[]={tim,dt,trec};
  double times_min[3],times_max[3];
  check_mpi(MPI_Allreduce(counts,counts_min,2,MPI_INT,MPI_MIN,MPI_COMM_WORLD),
            "Backup counter reduction failed.");
  check_mpi(MPI_Allreduce(counts,counts_max,2,MPI_INT,MPI_MAX,MPI_COMM_WORLD),
            "Backup counter reduction failed.");
  check_mpi(MPI_Allreduce(times,times_min,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD),
            "Backup time reduction failed.");
  check_mpi(MPI_Allreduce(times,times_max,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD),
            "Backup time reduction failed.");
  for (int i=0;i<2;i++){
    if (counts_min[i] != counts_max[i]) abort_run("Backup counters differ between ranks.");
  }
  for (int i=0;i<3;i++){
    if (times_min[i] != times_max[i]) abort_run("Backup times differ between ranks.");
  }
}

int MHD2D::bkup_(int flg)
try
{
  if (flg != 0 && flg != 1) abort_run("Invalid backup operation.");
  if (flg == 1) check_restart_state();
  const int local_status=(flg == 0)?
    bkup_load(val,nm,nd,&n,&cnt,&tim,&dt,&trec,mpi_rank,fildir.c_str()):
    bkup_save(val,nm,nd,n,cnt,tim,dt,trec,mpi_rank,fildir.c_str());
  int global_status=BKUP_ERROR;
  check_mpi(mpi_sync_status(local_status,&global_status),
            "Backup status reduction failed or ranks disagree.");
  if (global_status == BKUP_ERROR || (flg == 1 && global_status != BKUP_OK)){
    abort_run("Backup I/O failed.");
  }
  return global_status;
} catch (...){
  abort_run("Exception while loading or saving backup data.");
}

void MHD2D::prepare_run(int flg)
{
  int global_mode;
  check_mpi(mpi_sync_status(flg != 0,&global_mode),
            "Time integration modes differ between ranks.");
  const int status=bkup_(0);
  if (status == BKUP_OK){
    check_restart_state();
    if (!flg){
      // Recover bookkeeping only. Never change a restored time step.
      nrec=checked_nrec(std::round(dtrec/dt),nout);
      nmax=nrec*nout;
      if (!same_time(dtrec/dt,static_cast<double>(nrec)) ||
          n > nmax || cnt != n/nrec || !same_time(tim,n*dt)){
        abort_run("Fixed-step backup is incompatible with the output interval; dt was not changed.");
      }
    }
    dt_initialized=true;
  } else{
    initialize_dt();
  }
}
