#ifndef _CLASS_MYMPI_
#define _CLASS_MYMPI_

#define ABORT_INIT (1)

#include <iostream>
#include <mpi.h>
#include "common_mpi.h"

class MYMPI{

public:
  MYMPI(int* argc, char*** argv, int mnp)
  {
    // Constructor
    int ret;
    ret=MPI_Init(argc,argv);
    if (ret == MPI_SUCCESS){
      ret=MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
    }
    if (ret == MPI_SUCCESS){
      ret=MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    }
    if ((ret != MPI_SUCCESS) || (mpi_num != mnp)){
      MPI_Abort(MPI_COMM_WORLD,ABORT_INIT);
    } else{
      if (!mpi_rank) puts("MPI initialized successsfully.");
    }
  }
  virtual ~MYMPI()
  {
    // Destructor
    int ret=MPI_Finalize();
    if (ret == MPI_SUCCESS && !mpi_rank) puts("MPI finalized successsfully.");
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
};

#endif
