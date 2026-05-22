#Makefile to compile MD code
# ml foss/2025b
.SUFFIXES : .o .f90

FC_SERIAL = gfortran
FC_MPI    = mpif90

FFLAGS = -O3

TARGETS = exec exec_p2p exec_collective

all: $(TARGETS)

exec: modules_serial.o md_serial.o
	$(FC_SERIAL) $(FFLAGS) -o $@ $^

exec_p2p: modules_mpi.o md_mpi_p2p.o
	$(FC_MPI) $(FFLAGS) -o $@ $^

exec_collective: modules_mpi.o md_mpi_collective.o
	$(FC_MPI) $(FFLAGS) -o $@ $^

modules_serial.o: modules_serial.f90
	$(FC_SERIAL) $(FFLAGS) -c modules_serial.f90

md_serial.o: md_serial.f90 modules_serial.o
	$(FC_SERIAL) $(FFLAGS) -c md_serial.f90

modules_mpi.o: modules_mpi.f90
	$(FC_MPI) $(FFLAGS) -c modules_mpi.f90

md_mpi_p2p.o: md_mpi_p2p.f90 modules_mpi.o
	$(FC_MPI) $(FFLAGS) -c md_mpi_p2p.f90

md_mpi_collective.o: md_mpi_collective.f90 modules_mpi.o
	$(FC_MPI) $(FFLAGS) -c md_mpi_collective.f90

clean:
	rm -f *.o *.mod K*.xyz exec exec_p2p exec_collective
