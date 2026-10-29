include Makefile.inc

.PHONY: all common mpi clean

all: common mpi

common:
	$(MAKE) -C common

mpi:
	$(MAKE) -C mpi

lib$(LIBNAME).a:
	$(MAKE) -C common

lib$(LIBMPI).a:
	$(MAKE) -C mpi

clean:
	$(RM) *.a
	$(MAKE) -C common clean
	$(MAKE) -C mpi clean
