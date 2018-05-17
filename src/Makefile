# Fortran compilation

# The following 3 variables must be set in order to compile qumvia.
#HOME=/home/armino
#MKL_HOME=/opt/intel/composerxe/mkl
#LIO_HOME=$(HOME)/progs/lio

BINDIR=../bin
OBJDIR=../obj

FC=ifort
FFLAGS= -O3 -ip -xHost -FR -I${MKL_HOME}/include 
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME) 
LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(LIO_HOME)/g2g -lg2g -L$(LIO_HOME)/lioamber -llio-g2g
LIOFLAGS= -I$(LIO_HOME)/lioamber -I$(LIO_HOME)/lioamber/obj


cpu: clean mkdirs qumvia.cpu mvdirs

lio: clean mkdirs qumvia.lio mvdirs

mkdirs:
	if [ ! -d "$(BINDIR)" ]; then mkdir $(BINDIR); fi
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

qumvia.cpu: M_kracken.o qvamod_common.o qvamod_cpu.o qvamain_cpu.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)

qumvia.lio: M_kracken.o qvamod_common.o qvamod_lioexcl.o qvamod_lio.o qvamain_lio.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(LIOFLAGS) $(LIOLIBS)

#----------------------------------------------------------------------------------------
# KRACKEN COMMAND LINE PARSER
#----------------------------------------------------------------------------------------
M_kracken.o: M_kracken.f90
	$(FC) $(FFLAGS) -c M_kracken.f90

#----------------------------------------------------------------------------------------
# COMMON COMPILATION RULES
#----------------------------------------------------------------------------------------
qvamod_common.o: qvamod_common.f90 
	$(FC) $(FFLAGS) -c qvamod_common.f90 $(FLIBS)

#----------------------------------------------------------------------------------------
# CPU COMPILATION RULES
#----------------------------------------------------------------------------------------
qvamod_cpu.o: qvamod_common.mod qvamod_alt.f90 
	$(FC) -fpp -Dqvacpu $(FFLAGS) -c qvamod_alt.f90 $(FLIBS) -o qvamod_cpu.o

qvamain_cpu.o: qvamod_common.mod qvamod_cpu.mod qvamain.f90 
	$(FC) -fpp -Dqvacpu $(FFLAGS) -c qvamain.f90 $(FLIBS) -o qvamain_cpu.o

#----------------------------------------------------------------------------------------
# LIO COMPILATION RULES
#----------------------------------------------------------------------------------------
qvamod_lioexcl.o: qvamod_common.mod qvamod_lioexcl.f90
	$(FC) -fpp -Dqvalio $(FFLAGS) -c qvamod_lioexcl.f90 $(FLIBS) $(LIOLIBS) -o qvamod_lioexcl.o

qvamod_lio.o: qvamod_common.mod qvamod_lioexcl.mod qvamod_alt.f90 
	$(FC) -fpp -Dqvalio $(FFLAGS) -c qvamod_alt.f90 $(FLIBS) $(LIOLIBS) \
		-o qvamod_lio.o

qvamain_lio.o: qvamod_common.mod qvamod_lioexcl.mod qvamod_lio.mod qvamain.f90 
	$(FC) -fpp -Dqvalio $(FFLAGS) $(LIOFLAGS) -c qvamain.f90 $(FLIBS) $(LIOLIBS) \
		-o qvamain_lio.o

#----------------------------------------------------------------------------------------

mvdirs:
	mv *.o *.mod $(OBJDIR)
	mv qumvia.* $(BINDIR)

clean:
	rm -rf *.o *.mod qumvia.* 
