#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>

namespace {

static_assert(sizeof(float) == 4, "Merged data require 32-bit floats");

[[noreturn]] void fail(const std::string& path, const std::string& message)
{
  throw std::runtime_error(path+": "+message);
}

std::string join(const std::string& directory, const std::string& name)
{
  return directory+(directory.back() == '/' ? "" : "/")+name;
}

std::string index_name(const std::string& prefix, std::size_t index)
{
  std::ostringstream out;
  out << prefix << std::setfill('0') << std::setw(5) << index;
  return out.str();
}

bool file_info(const std::string& path, struct stat& info)
{
  if (stat(path.c_str(),&info) == 0) return true;
  if (errno == ENOENT) return false;
  fail(path,std::strerror(errno));
}

std::size_t file_size(const std::string& path)
{
  struct stat info;
  if (!file_info(path,info)) fail(path,"file not found");
  if (!S_ISREG(info.st_mode) || info.st_size < 0) fail(path,"not a regular file");
  if (static_cast<unsigned long long>(info.st_size) > std::numeric_limits<std::size_t>::max()){
    fail(path,"file too large");
  }
  return static_cast<std::size_t>(info.st_size);
}

std::size_t product(std::size_t a, std::size_t b)
{
  if (b && a > std::numeric_limits<std::size_t>::max()/b){
    throw std::runtime_error("array size overflow");
  }
  return a*b;
}

template<class T> std::vector<T> read_text(const std::string& path)
{
  (void)file_size(path);
  std::ifstream input(path);
  if (!input) fail(path,"cannot open for reading");
  std::vector<T> values;
  std::string token;
  while (input >> token){
    std::istringstream parser(token);
    T value;
    if (!(parser >> value) || !parser.eof() || !std::isfinite(static_cast<double>(value))){
      fail(path,"invalid number '"+token+"'");
    }
    values.push_back(value);
  }
  if (input.bad() || !input.eof()) fail(path,"read failed");
  input.clear();
  input.close();
  if (input.fail()) fail(path,"close failed");
  if (values.empty()) fail(path,"empty numeric file");
  return values;
}

template<class T> std::vector<T> read_binary(const std::string& path, std::size_t count)
{
  const std::size_t bytes=product(count,sizeof(T));
  if (file_size(path) != bytes) fail(path,"unexpected file size; expected "+std::to_string(bytes)+" bytes");
  if (bytes > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())){
    fail(path,"read exceeds stream size limit");
  }
  std::ifstream input(path,std::ios::binary);
  if (!input) fail(path,"cannot open for reading");
  std::vector<T> values(count);
  if (bytes && !input.read(reinterpret_cast<char*>(values.data()),bytes)) fail(path,"short/failed read");
  if (input.peek() != std::char_traits<char>::eof() || input.bad()) fail(path,"file changed or read failed");
  input.clear();
  input.close();
  if (input.fail()) fail(path,"close failed");
  return values;
}

// POSIX temporary file in the destination directory: publish only after close succeeds.
// This protects individual files, not the entire multi-file output set.
void write_atomic(const std::string& path, const void* data, std::size_t bytes)
{
  std::string temporary=path+".tmp.XXXXXX";
  std::vector<char> name(temporary.begin(),temporary.end());
  name.push_back('\0');
  const int fd=mkstemp(name.data());
  if (fd < 0) fail(path,std::strerror(errno));
  std::FILE* output=fdopen(fd,"wb");
  if (!output){
    const int error=errno;
    close(fd);
    std::remove(name.data());
    fail(path,std::strerror(error));
  }
  const bool written=(bytes == 0 || std::fwrite(data,1,bytes,output) == bytes);
  const bool closed=(std::fclose(output) == 0);
  if (!written || !closed){
    std::remove(name.data());
    fail(path,"write/close failed");
  }
  if (std::rename(name.data(),path.c_str()) != 0){
    const int error=errno;
    std::remove(name.data());
    fail(path,"rename failed: "+std::string(std::strerror(error)));
  }
}

void copy_file(const std::string& source, const std::string& destination)
{
  struct stat src,dst;
  if (!file_info(source,src)) fail(source,"file not found");
  if (file_info(destination,dst) && src.st_dev == dst.st_dev && src.st_ino == dst.st_ino) return;
  const auto bytes=read_binary<char>(source,file_size(source));
  write_atomic(destination,bytes.data(),bytes.size());
}

struct Axis{
  std::size_t offset;
  std::vector<std::size_t> sizes,starts;
  std::vector<double> coordinates;
};

Axis read_axis(const std::string& directory, const std::string& axis,
               std::size_t ranks, std::size_t offset)
{
  Axis result;
  result.offset=offset;
  for (std::size_t r=0;r<ranks;r++){
    const std::string path=join(directory,index_name(axis+"_",r)+".dat");
    const auto values=read_text<double>(path);
    if (values.size() <= product(2,offset)) fail(path,"no physical cells after removing ghosts");
    for (std::size_t i=1;i<values.size();i++){
      if (values[i] <= values[i-1]) fail(path,"coordinates must increase strictly");
    }
    if (!result.coordinates.empty() && values[offset] <= result.coordinates.back()){
      fail(path,"overlapping or unordered rank coordinates");
    }
    result.sizes.push_back(values.size());
    result.starts.push_back(result.coordinates.size());
    result.coordinates.insert(result.coordinates.end(),values.begin()+offset,values.end()-offset);
  }
  return result;
}

void write_axis(const std::string& path, const Axis& axis)
{
  std::ostringstream text;
  text << std::fixed << std::setprecision(6); // Preserve the original coordinate format.
  for (double value : axis.coordinates) text << value << '\n';
  const std::string data=text.str();
  write_atomic(path,data.data(),data.size());
}

struct FaceFields{
  std::vector<float> bx,by; // (nyall,nxall+1), (nyall+1,nxall)
};

// prefix is either outdat_NNNNN_ (eight fields) or g_potential_ (one field).
// For MHD, retain the CT faces separately and put two-face averages in fields 4,5.
std::vector<float> merge_fields(const std::string& directory, const std::string& prefix,
                                const Axis& x, const Axis& y, std::size_t fields,
                                bool require_finite, FaceFields* faces=nullptr)
{
  const std::size_t nxall=x.coordinates.size(),nyall=y.coordinates.size();
  const std::size_t global_cells=product(nxall,nyall);
  std::vector<float> merged(product(fields,global_cells));
  if (faces){
    faces->bx.resize(product(nxall+1,nyall));
    faces->by.resize(product(nxall,nyall+1));
  }
  for (std::size_t ry=0;ry<y.sizes.size();ry++){
    for (std::size_t rx=0;rx<x.sizes.size();rx++){
      const std::size_t rank=ry*x.sizes.size()+rx;
      const std::string path=join(directory,index_name(prefix,rank)+".dat");
      const std::size_t nx=x.sizes[rx],ny=y.sizes[ry],local_cells=product(nx,ny);
      const auto values=read_binary<float>(path,product(fields,local_cells));
      if (require_finite){
        for (float value : values){
          if (!std::isfinite(value)) fail(path,"nonfinite gravitational potential");
        }
      }
      for (std::size_t s=0;s<fields;s++){
        for (std::size_t j=y.offset;j<ny-y.offset;j++){
          const std::size_t source=s*local_cells+nx*j+x.offset;
          const std::size_t dest=s*global_cells+nxall*(y.starts[ry]+j-y.offset)+x.starts[rx];
          std::copy_n(values.data()+source,nx-2*x.offset,merged.data()+dest);
        }
      }
      if (faces){
        const std::size_t width=nx-2*x.offset,height=ny-2*y.offset;
        // Shared interfaces belong to the rank on their positive side. Only
        // the final rank contributes the global right/top face from its halo.
        const std::size_t xfaces=width+(rx+1 == x.sizes.size());
        const std::size_t yfaces=height+(ry+1 == y.sizes.size());
        for (std::size_t j=0;j<height;j++){
          std::copy_n(values.data()+4*local_cells+nx*(j+y.offset)+x.offset,xfaces,
                      faces->bx.data()+(nxall+1)*(y.starts[ry]+j)+x.starts[rx]);
        }
        for (std::size_t j=0;j<yfaces;j++){
          std::copy_n(values.data()+5*local_cells+nx*(j+y.offset)+x.offset,width,
                      faces->by.data()+nxall*(y.starts[ry]+j)+x.starts[rx]);
        }
      }
    }
  }
  if (faces){
    for (std::size_t j=0;j<nyall;j++){
      for (std::size_t i=0;i<nxall;i++){
        const std::size_t cell=nxall*j+i,bxf=(nxall+1)*j+i,byf=cell;
        // Average in double precision; keep the existing float32 output format.
        merged[4*global_cells+cell]=0.5*faces->bx[bxf]+0.5*faces->bx[bxf+1];
        merged[5*global_cells+cell]=0.5*faces->by[byf]+0.5*faces->by[byf+nxall];
      }
    }
  }
  return merged;
}

void merge(const std::string& source, const std::string& destination)
{
  for (const auto& directory : {source,destination}){
    struct stat info;
    if (directory.empty() || !file_info(directory,info) || !S_ISDIR(info.st_mode)){
      fail(directory,"directory does not exist");
    }
  }
  const auto offsets=read_text<long long>(join(source,"offsets.dat"));
  const auto ranks=read_text<long long>(join(source,"mpinum.dat"));
  if (offsets.size() != 2 || offsets[0] < 0 || offsets[1] < 0 ||
      offsets[0] > std::numeric_limits<int>::max() || offsets[1] > std::numeric_limits<int>::max()){
    fail(join(source,"offsets.dat"),"expected two nonnegative int offsets");
  }
  if (offsets[0] == 0 || offsets[1] == 0){
    fail(join(source,"offsets.dat"),"CT output requires at least one ghost cell in each direction");
  }
  if (ranks.size() != 2 || ranks[0] <= 0 || ranks[1] <= 0 ||
      ranks[0] > std::numeric_limits<int>::max() || ranks[1] > std::numeric_limits<int>::max()){
    fail(join(source,"mpinum.dat"),"expected two positive int rank counts");
  }
  const std::size_t nranks=product(ranks[0],ranks[1]);
  if (nranks > static_cast<std::size_t>(std::numeric_limits<int>::max())){
    fail(join(source,"mpinum.dat"),"rank count exceeds int range");
  }
  const Axis x=read_axis(source,"x",ranks[0],offsets[0]);
  const Axis y=read_axis(source,"y",ranks[1],offsets[1]);
  const auto times=read_text<double>(join(source,"t.dat"));
  if (times.front() < 0 || !std::is_sorted(times.begin(),times.end())){
    fail(join(source,"t.dat"),"expected nonnegative, nondecreasing times");
  }
  const auto parameters=read_text<double>(join(source,"params.dat"));
  if (parameters[0] <= 1) fail(join(source,"params.dat"),"gamma must exceed one");

  std::size_t potentials=0;
  for (std::size_t rank=0;rank<nranks;rank++){
    struct stat info;
    if (file_info(join(source,index_name("g_potential_",rank)+".dat"),info)) ++potentials;
  }
  if (potentials != 0 && potentials != nranks) fail(source,"gravitational potential is missing for some ranks");
  const std::string potential_path=join(destination,"merge_g_potential.dat");
  struct stat old_potential;
  if (!potentials && file_info(potential_path,old_potential)){
    fail(potential_path,"stale potential: input has no rank potentials; use a clean output directory");
  }
  // Preflight every expected binary size before replacing any existing output.
  for (std::size_t rank=0;rank<nranks;rank++){
    const std::size_t cells=product(x.sizes[rank%ranks[0]],y.sizes[rank/ranks[0]]);
    for (std::size_t n=0;n<times.size();n++){
      const std::string path=join(source,index_name(index_name("outdat_",n)+"_",rank)+".dat");
      if (file_size(path) != product(product(8,cells),sizeof(float))) fail(path,"unexpected MHD file size");
    }
  }
  if (potentials){
    const auto phi=merge_fields(source,"g_potential_",x,y,1,true);
    write_atomic(potential_path,phi.data(),product(phi.size(),sizeof(float)));
    std::cout << "Merged gravitational potential.\n";
  } else{
    std::cout << "No rank potentials: pressure will not be gravity-corrected.\n";
  }
  write_axis(join(destination,"merge_x.dat"),x);
  write_axis(join(destination,"merge_y.dat"),y);
  copy_file(join(source,"t.dat"),join(destination,"t.dat"));
  copy_file(join(source,"params.dat"),join(destination,"params.dat"));
  for (std::size_t n=0;n<times.size();n++){
    FaceFields faces;
    const auto values=merge_fields(source,index_name("outdat_",n)+"_",x,y,8,false,&faces);
    write_atomic(join(destination,index_name("merge_bx_face_",n)+".dat"),
                 faces.bx.data(),product(faces.bx.size(),sizeof(float)));
    write_atomic(join(destination,index_name("merge_by_face_",n)+".dat"),
                 faces.by.data(),product(faces.by.size(),sizeof(float)));
    write_atomic(join(destination,index_name("merge_outdat_",n)+".dat"),
                 values.data(),product(values.size(),sizeof(float)));
    std::cout << "Merged time index " << n << '\n';
  }
}

}

int main(int argc, char* argv[])
{
  if (argc != 2 && argc != 3){
    std::cerr << "Usage: " << argv[0] << " INPUT_DIRECTORY [OUTPUT_DIRECTORY]\n"
              << "If omitted, OUTPUT_DIRECTORY is INPUT_DIRECTORY.\n";
    return EXIT_FAILURE;
  }
  try{
    merge(argv[1],argc == 3 ? argv[2] : argv[1]);
    return EXIT_SUCCESS;
  } catch (const std::exception& error){
    std::cerr << "Merge failed: " << error.what() << '\n';
  } catch (...){
    std::cerr << "Merge failed: unknown error.\n";
  }
  return EXIT_FAILURE;
}
