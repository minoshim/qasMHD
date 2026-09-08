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
OBJS = main.o mhd2d_class.o mhd2d_init_.o mhd2d_paras.o mhd2d_solve.o
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
