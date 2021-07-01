EXECS=mpi_hot_potato
MPICC?=mpicc
all:${EXECS}
mpi_hot_potato:mpi_hot_potato.c
	${MPICC} -o mpi_hot_potato mpi_hot_potato.c
clean:
	rm -f ${EXECS}
