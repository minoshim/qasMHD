all:
	cd common; make
	cd mpi; make
clean:
	rm -f *.a
	cd common; make clean
	cd mpi; make clean
