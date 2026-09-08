DIR0 = ../../..
include $(DIR0)/Makefile.inc
DIRI = $(DIR0)/common
DIRM = $(DIR0)/mpi
DIRL = $(DIR0)
CPPFLAGS += -I$(DIRI) -I$(DIRM)
LDFLAGS += -L$(DIRL)
LDLIBS += -lm -l$(LIBNAME) -l$(LIBMPI)
LIBRARIES = $(DIRL)/lib$(LIBNAME).a $(DIRL)/lib$(LIBMPI).a
TARGET = a.out
OBJS = main.o dmhd3d_class.o dmhd3d_init_.o dmhd3d_solve.o mhd3d_class.o mhd3d_paras.o mhd3d_solve.o
DEPS = $(OBJS:.o=.d)
DATADIR = dat

.PHONY: all clean cdata

all: $(TARGET)

$(TARGET): $(OBJS) $(LIBRARIES)
	$(MPICXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

%.o: %.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(DIRL)/lib$(LIBNAME).a:
	$(MAKE) -C $(DIRL)/common

$(DIRL)/lib$(LIBMPI).a:
	$(MAKE) -C $(DIRL)/mpi

clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

cdata:
	$(RM) $(DATADIR)/*.dat

-include $(DEPS)
