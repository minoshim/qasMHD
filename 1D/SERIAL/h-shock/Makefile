DIR0 = ../../..
include $(DIR0)/Makefile.inc
DIRI = $(DIR0)/common
DIRL = $(DIR0)
INCDIR  = -I$(DIRI)
LDFLAGS = -L$(DIRL)
LIBS = -lm -l$(LIBNAME)
LIBRARY = $(DIRL)/lib$(LIBNAME).a
BASECLASS = $(DIRI)/mhd_class.hpp
TARGET = a.out
OBJS = main.o hmhd1d_class.o hmhd1d_paras.o hmhd1d_solve.o mhd1d_class.o mhd1d_init_.o mhd1d_solve.o
DATADIR = dat

.PHONY : all
.PHONY : clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $^ $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@

%.o: %.c
	$(CC) -c $< $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@
%.o: %.cpp
	$(CC) -c $< $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $@

$(LIBRARY):
	cd $(DIRL); make

clean:
	rm -f $(TARGET) $(OBJS)

cdata:
	rm -f $(DATADIR)/*.dat

main.o: $(LIBRARY) $(BASECLASS) hmhd1d_class.o hmhd1d_paras.o hmhd1d_solve.o mhd1d_class.o mhd1d_init_.o mhd1d_solve.o hmhd1d_class.hpp mhd1d_class.hpp mymacros.hpp
hmhd1d_class.o: $(LIBRARY) $(BASECLASS) mhd1d_class.o hmhd1d_class.hpp mhd1d_class.hpp mymacros.hpp
hmhd1d_solve.o: $(LIBRARY) $(BASECLASS) hmhd1d_class.o mhd1d_class.o hmhd1d_class.hpp mhd1d_class.hpp mymacros.hpp
hmhd1d_paras.o: $(BASECLASS) hmhd1d_class.o mhd1d_class.o hmhd1d_class.hpp mhd1d_class.hpp mymacros.hpp
mhd1d_class.o: $(LIBRARY) $(BASECLASS) mhd1d_class.hpp mymacros.hpp
mhd1d_init_.o: $(LIBRARY) $(BASECLASS) mhd1d_class.o mhd1d_class.hpp mymacros.hpp
mhd1d_solve.o: $(LIBRARY) $(BASECLASS) mhd1d_class.o mhd1d_class.hpp mymacros.hpp
