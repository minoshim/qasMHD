CXX = g++
MPICXX = mpic++

CXXSTD = -std=c++11
OPTFLAGS = -O2
#OPTFLAGS = -O3 -march=native ## Aggressive optimization
OMPFLAGS = -fopenmp
WARNFLAGS =
DEPFLAGS = -MMD -MP
LTOFLAGS =
#LTOFLAGS = -flto ## Link-time optimization

CPPFLAGS =
CXXFLAGS = $(CXXSTD) $(OPTFLAGS) $(OMPFLAGS) $(WARNFLAGS) $(DEPFLAGS) $(LTOFLAGS)
LDFLAGS = $(OMPFLAGS) $(LTOFLAGS)

AR = ar
ARFLAGS = rcs

LIBNAME = qasmhd
LIBMPI = mympi
