include Makefile.inc

all:
	cd common; make
	cd mpi; make
lib$(LIBNAME).a:
	cd common; make
lib$(LIBMPI).a:
	cd mpi; make
clean:
	rm -f *.a
	cd common; make clean
	cd mpi; make clean

#test
