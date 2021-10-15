DIR0 = ..
include $(DIR0)/Makefile.inc
LIBS = -lm
TARGET = $(DIR0)/lib$(LIBNAME).a
OBJS = common_func.o common_interp.o mhd_eigen.o mhd_flux.o mhd_func.o

.PHONY : all
.PHONY : clean

all: $(TARGET)

$(TARGET): $(OBJS)
	ar -r $(TARGET) $(OBJS)

%.o: %.c
	$(CC) -c $< $(CFLAGS) $(LIBS) -o $@

clean:
	rm -f $(OBJS)

common_func.o: common_func.h
common_interp.o: common_interp.h common_func.o common_func.h
mhd_eigen.o: mhd_eigen.h common_func.o common_func.h
mhd_flux.o: mhd_flux.h common_func.o common_func.h
mhd_func.o: mhd_func.h common_func.o common_func.h mhd_eigen.o mhd_eigen.h
