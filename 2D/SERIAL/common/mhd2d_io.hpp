#ifndef QASMHD_SERIAL2D_IO_HPP
#define QASMHD_SERIAL2D_IO_HPP

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace serial2d_io {

// Same path/open/close checks as MPI, but failures stop this serial process.
[[noreturn]] inline void output_error(const char *operation, const std::string& path)
{
  const int error_number=errno;
  const char *reason=(error_number != 0)?std::strerror(error_number):"unknown I/O error";
  std::fprintf(stderr,"Output error: %s failed for '%s': %s\n",
               operation,path.c_str(),reason);
  std::exit(EXIT_FAILURE);
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
