FC=gfortran
FFLAG=-fbounds-check -lfftw3 -lm
objects=do_3COM_AUCC.o filter.o do_ncc.o do_pcc.o pcc.o analytic.o smooth.o do_norm.o do_agc.o taperf.o do_whiten.o do_water.o smoothf.o globe_data.o sacio.o usage.o
all:sacio.mod globe_data.mod ../bin/do_3COM_AUCC
%.o:%.f90
	$(FC) $^ -c $(FFLAG)
sacio.mod:sacio.f90
	$(FC) $^ -c $(FFLAG)
globe_data.mod:globe_data.f90
	$(FC) $^ -c $(FFLAG)
../bin/do_3COM_AUCC:$(objects)
	$(FC) $^ -o $@ $(FFLAG)
clean:
	-rm *.o *.mod
