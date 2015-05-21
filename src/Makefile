# Fortran compilation

HOME=/home/diegoa
#MKL_HOME=/opt/intel/composerxe/mkl
MKL_HOME=/home/diegoa/intel/composer_xe_2013.0.079/mkl
MKL=$(MKL_HOME)
MKL_PROCESSOR=intel64
BASEDIR=$(HOME)/dev/QUMVIA
SRCDIR=$(BASEDIR)/src
BINDIR=../bin
OBJDIR=../obj

FC=ifort
FFLAGS= -O3 -ip -xHost -FR -I${MKL_HOME}/include 
#FFLAGS= -O3 -ipo -xHost -FR -I${MKL_HOME}/include 
#FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_thread -L$(MKL_HOME) 
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME) 

LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(HOME)/progs/lio/g2g/ -lg2g -L$(HOME)/progs/lio/lioamber -llio-g2g

cpu: clean mkdirs qumvia.cpu mvdirs

lio: clean mkdirs qumvia.lio mvdirs

mkdirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

qumvia.cpu: M_kracken.o qvamod_common.o qvamod_cpu.o qvamain_cpu.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)

qumvia.lio: M_kracken.o qvamod_common.o qvamod_lio.o qvamain_cpu.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(LIOLIBS)

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
	$(FC) -fpp -Dqvalio $(FFLAGS) -c qvamod_lioexcl.f90 $(FLIBS) $(LIOLIBS) \
		-o qvamod_lioexcl.o

qvamod_lio.o: qvamod_common.mod qvamod_lioexcl.mod M_kracken.mod qvamod_alt.f90 
	$(FC) -fpp -Dqvalio $(FFLAGS) -c qvamod_alt.f90 $(FLIBS) $(LIOLIBS) \
		-o qvamod_lio.o

qvamain_lio.o: qvamod_common.mod qvamod_lioexcl.mod qvamod_lio.mod M_kracken.mod qvamain.f90 
	$(FC) -fpp -Dqvalio $(FFLAGS) -c qvamain.f90 $(FLIBS) $(LIOLIBS) \
		-o qvamain_lio.o

#----------------------------------------------------------------------------------------

mvdirs:
	mv *.o *.mod $(OBJDIR)
	mv qumvia.* $(BINDIR)

clean:
	rm -rf *.o *.mod $(OBJDIR) $(BINDIR)
