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
OBJS = main.o dmhd2d_class.o dmhd2d_init_.o dmhd2d_solve.o mhd2d_class.o mhd2d_paras.o mhd2d_solve.o
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

main.o: $(LIBRARY) $(BASECLASS) dmhd2d_class.o dmhd2d_init_.o dmhd2d_solve.o mhd2d_class.o mhd2d_paras.o mhd2d_solve.o dmhd2d_class.hpp mhd2d_class.hpp mymacros.hpp
dmhd2d_class.o: $(LIBRARY) $(BASECLASS) mhd2d_class.o dmhd2d_class.hpp mhd2d_class.hpp mymacros.hpp
dmhd2d_init_.o: $(LIBRARY) $(BASECLASS) dmhd2d_class.o mhd2d_class.o dmhd2d_class.hpp mhd2d_class.hpp mymacros.hpp
dmhd2d_solve.o: $(LIBRARY) $(BASECLASS) dmhd2d_class.o mhd2d_class.o dmhd2d_class.hpp mhd2d_class.hpp mymacros.hpp
mhd2d_class.o: $(LIBRARY) $(BASECLASS) mhd2d_class.hpp mymacros.hpp
mhd2d_paras.o: $(BASECLASS) mhd2d_class.o mhd2d_class.hpp mymacros.hpp
mhd2d_solve.o: $(LIBRARY) $(BASECLASS) mhd2d_class.o mhd2d_class.hpp mymacros.hpp
