DIR0 = ../../..
include $(DIR0)/Makefile.inc
DIRI = $(DIR0)/common
DIRL = $(DIR0)
INCDIR  = -I$(DIRI)
LDFLAGS = -L$(DIRL)
LIBS = -lm -l$(LIBNAME)
LIBRARY = $(DIRL)/lib$(LIBNAME).a
TARGET = a.out
OBJS = main.o mhd_fd2d.o boundary.o
DATADIR = dat

.PHONY : all
.PHONY : clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INCDIR) -o $(TARGET)

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

main.o: $(LIBRARY) myfunc.h mhd_fd2d.o mhd_fd2d.h boundary.o boundary.h dataio.hpp get.hpp global.hpp init.hpp new_delete.hpp
mhd_fd2d.o: $(LIBRARY) mhd_fd2d.h boundary.o boundary.h
boundary.o: $(LIBRARY) boundary.h
