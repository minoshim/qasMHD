DIR0 = ..
include $(DIR0)/Makefile.inc
TARGET = $(DIR0)/lib$(LIBMPI).a
OBJS = common_mpi.o

.PHONY : all
.PHONY : clean

all: $(TARGET)

$(TARGET): $(OBJS)
	ar -r $@ $^

%.o: %.c
	$(MPICC) -c $< $(CFLAGS) -o $@

clean:
	rm -f $(OBJS)

common_mpi.o: common_mpi.h
