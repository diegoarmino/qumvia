# Fortran compilation

HOME=/home/armino
MKL_HOME=/opt/intel/composerxe/mkl
BASEDIR=$(HOME)/dev/QUMVIA
SRCDIR=$(BASEDIR)/src
BINDIR=../bin
OBJDIR=../obj

FC=ifort
#FFLAGS= -O3 -ip -xHost -FR -I${MKL_HOME}/include 
FFLAGS= -O3 -ipo -xHost -FR -I${MKL_HOME}/include 
#FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_thread -L$(MKL_HOME) 
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME) 

MKL=/opt/intel/composerxe/mkl
MKL_PROCESSOR=intel64
LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(HOME)/progs/lio/g2g/ -lg2g -L$(HOME)/progs/lio/lioamber -llio-g2g

cpu: clean mkdirs qumvia.cpu mvdirs

lio: clean mkdirs qumvia.lio mvdirs

mkdirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

qumvia.cpu: qvamod_common.o qvamod_cpu.o qvamain_cpu.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)

qumvia.lio: qvamod_common.o qvamod_lio.o qvamain_lio.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(LIOLIBS)

qvamod_common.o: qvamod_common.f90 
	$(FC) $(FFLAGS) -c qvamod_common.f90 $(FLIBS)

qvamod_cpu.o: qvamod_common.mod qvamod_cpu.f90 
	$(FC) $(FFLAGS) -c qvamod_cpu.f90 $(FLIBS) 

qvamod_lio.o: qvamod_common.mod qvamod_lio.f90 
	$(FC) $(FFLAGS) -c qvamod_lio.f90 $(FLIBS) $(LIOLIBS)

qvamain_cpu.o: qvamod_common.mod qvamod_cpu.mod qvamain.f90 
	$(FC) -fpp $(FFLAGS) -Dqvacpu -c qvamain.f90 $(FLIBS) -o qvamain_cpu.o

qvamain_lio.o: qvamod_common.mod qvamod_lio.mod qvamain.f90 
	$(FC) -fpp $(FFLAGS) -Dqvalio -c qvamain.f90 $(FLIBS) $(LIOLIBS) -o qvamain_lio.o

mvdirs:
	mv *.o *.mod $(OBJDIR)
	mv qumvia.* $(BINDIR)

clean:
	rm -rf *.o *.mod $(OBJDIR) $(BINDIR)
