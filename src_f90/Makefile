#############################################################################
# Fortran compiler options and directives

KROME_F90_DIR ?= ./krome/chem_net_ramses_rt/C_2_O_4

VPATH = $(KROME_F90_DIR):./

F90	= mpif90 -ffree-line-length-none -frecord-marker=4 -O3 -g -fbacktrace -fbounds-check -fallow-argument-mismatch
FFLAGS  = -O3 -cpp
LDFLAGS = -lmpi


OBJS = module_parallel_mpi.o opkda2.o opkda1.o opkdmain.o krome_user_commons.o krome_all.o module_constants.o module_spectra.o module_ramses.o module_krome.o module_cell.o module_file.o module_worker.o module_master.o
#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#rule for f
%.o:%.f
	$(F90) -ffree-line-length-none -w -c $^ -o $@

#############################################################################
ion_fractions: $(OBJS) ion_fractions.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) ion_fractions.o

clean:
	rm -f $(OBJS) *.o *.mod
#############################################################################
