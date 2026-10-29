#include "mhd2d_class.hpp"
#include "mhd2d_io.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace serial2d_io;

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
    // Same initial output-interval alignment as 2D/MPI.
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
