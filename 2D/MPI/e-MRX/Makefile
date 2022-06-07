DIR0 = ../../..
include $(DIR0)/Makefile.inc
DIRI = $(DIR0)/common
DIRM = $(DIR0)/mpi
DIRL = $(DIR0)
INCDIR  = -I$(DIRI) -I$(DIRM)
LDFLAGS = -L$(DIRL)
LIBS = -lm -l$(LIBNAME) -l$(LIBMPI)
LIBRARY = $(DIRL)/lib$(LIBNAME).a $(DIRL)/lib$(LIBMPI).a 
TARGET = a.out
OBJS = main.o mhd_fd2d.o emhd_func2d.o boundary.o
DATADIR = dat

.PHONY : all
.PHONY : clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(MPICC) $^ $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@

%.o: %.c
	$(MPICC) -c $< $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@
%.o: %.cpp
	$(MPICC) -c $< $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@

$(LIBRARY):
	cd $(DIRL); make

clean:
	rm -f $(TARGET) $(OBJS)

cdata:
	rm -f $(DATADIR)/*.dat

main.o: $(LIBRARY) myfunc.h mhd_fd2d.o mhd_fd2d.h emhd_func2d.o emhd_func2d.h boundary.o boundary.h dataio.hpp get.hpp global.hpp init.hpp new_delete.hpp
mhd_fd2d.o: $(LIBRARY) mhd_fd2d.h emhd_func2d.o emhd_func2d.h boundary.o boundary.h
emhd_func2d.o: $(LIBRARY) emhd_func2d.h mhd_fd2d.h boundary.o boundary.h
boundary.o: $(LIBRARY) boundary.h
