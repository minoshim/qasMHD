#ifndef QASMHD_MPI2D_IO_HPP
#define QASMHD_MPI2D_IO_HPP

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <mpi.h>

namespace mpi2d_io {

[[noreturn]] inline void abort_run(const char *message, int error_code=EXIT_FAILURE)
{
  int rank=-1;
  (void)MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::fprintf(stderr,"Rank %d: %s\n",rank,message);
  (void)MPI_Abort(MPI_COMM_WORLD,error_code);
  std::abort();
}

inline void check_mpi(int error_code, const char *operation)
{
  if (error_code != MPI_SUCCESS) abort_run(operation,error_code);
}

// Same path/open/close checks as SERIAL, but failures stop all MPI ranks.
[[noreturn]] inline void output_error(const char *operation, const std::string& path)
{
  const int error_number=errno;
  const char *reason=(error_number != 0)?std::strerror(error_number):"unknown I/O error";
  int rank=-1;
  (void)MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::fprintf(stderr,"Rank %d: Output error: %s failed for '%s': %s\n",
               rank,operation,path.c_str(),reason);
  (void)MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  std::abort();
}

inline std::string output_path(const std::string& directory, const std::string& filename)
{
  std::string path(directory);
  if (!path.empty() && path.back() != '/') path.push_back('/');
  return path+filename;
}

inline std::FILE *open_output(const std::string& path, const char *mode)
{
  errno=0;
  std::FILE *outfil=std::fopen(path.c_str(),mode);
  if (outfil == nullptr) output_error("fopen",path);
  return outfil;
}

inline void close_output(std::FILE *outfil, const std::string& path)
{
  errno=0;
  if (std::fclose(outfil) != 0) output_error("fclose",path);
}

}
#endif

