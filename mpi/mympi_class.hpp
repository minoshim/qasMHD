#ifndef _CLASS_MYMPI_
#define _CLASS_MYMPI_

#define ABORT_INIT (1)

#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include "common_mpi.h"

class MYMPI{

private:
  [[noreturn]] static void abort_all(int error_code) noexcept
  {
    (void)MPI_Abort(MPI_COMM_WORLD,error_code);
    std::abort();
  }

public:
  MYMPI(int* argc, char*** argv, int mnp)
  {
    // Constructor
    int provided=MPI_THREAD_SINGLE;
    int ret=MPI_Init_thread(argc,argv,MPI_THREAD_SERIALIZED,&provided);
    if (ret != MPI_SUCCESS){
      std::fprintf(stderr,"MPI_Init_thread failed with error code %d.\n",ret);
      std::abort();
    }
    mpi_initialized=true;

    if (ret == MPI_SUCCESS){
      ret=MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
    }
    if (ret == MPI_SUCCESS){
      ret=MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    }
    if ((ret != MPI_SUCCESS) || (provided < MPI_THREAD_SERIALIZED) || (mpi_num != mnp)){
      abort_all((ret == MPI_SUCCESS)?ABORT_INIT:ret);
    } else{
      if (!mpi_rank) std::puts("MPI initialized successsfully.");
    }
  }

  MYMPI(const MYMPI&)=delete;
  MYMPI& operator=(const MYMPI&)=delete;
  MYMPI(MYMPI&&)=delete;
  MYMPI& operator=(MYMPI&&)=delete;

  virtual ~MYMPI() noexcept
  {
    // Destructor
    if (!mpi_initialized) return;

    int finalized=0;
    int ret=MPI_Finalized(&finalized);
    if (ret != MPI_SUCCESS) abort_all(ret);
    if (finalized){
      mpi_initialized=false;
      return;
    }

    ret=MPI_Finalize();
    if (ret != MPI_SUCCESS){
      std::fprintf(stderr,"MPI_Finalize failed with error code %d.\n",ret);
      abort_all(ret);
    }
    mpi_initialized=false;
    if (!mpi_rank) std::puts("MPI finalized successsfully.");
  }
  int get_mrank()
  {
    // Get my MPI rank
    return mpi_rank;
  }
  int get_mnum()
  {
    // Get number of MPI processes
    return mpi_num;
  }
  
protected:
  int mpi_num=1;		// Number of MPI processes
  int mpi_rank=0;		// MPI rank

private:
  bool mpi_initialized=false;
};

#endif
