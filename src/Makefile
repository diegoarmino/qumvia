# Fortran compilation

# The following 3 variables must be set in order to compile qumvia.
#HOME=/home/armino
#MKL_HOME=/opt/intel/composerxe/mkl
#LIO_HOME=$(HOME)/progs/lio

BINDIR=../bin
OBJDIR=../obj

FC=ifort
<<<<<<< 4e69ba7e2d5240265cd90d98ec9b794eb73bd52a
<<<<<<< 3d02e19adaf35b093b863312b609077861e2d1cd
<<<<<<< 8f89ea86e25185c5d553f4bcb7a1cc626cb9fd1b
<<<<<<< cb52330281e84ee344e85760f1b200e1508f82b8
FFLAGS= -O3 -ip -xHost -FR -I${MKL_HOME}/include
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME)
=======
FFLAGS= -O3 -ip -xHost -FR -I${MKL_HOME}/include 
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME) 
LIOFLAGS= -shared
>>>>>>> Added Resonant Raman Intesities code to qvamod_lioexcl
LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(LIO_HOME)/g2g -lg2g -L$(LIO_HOME)/lioamber -llio-g2g
LIOFLAGS= -I$(LIO_HOME)/lioamber -I$(LIO_HOME)/lioamber/obj

=======
FFLAGS= -O3 -ip -xHost -FR -fPIC -I${MKL_HOME}/include 
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME) 
LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(LIO_HOME)/g2g -lg2g -L$(LIO_HOME)/lioamber -llio-g2g
LIOINC=$(LIO_HOME)/lioamber/obj/garcha_g2g_obj/mods
<<<<<<< 6c5d4a25bd13af90767d575565abfcd127171896
LIOFLAGS= -shared -I$(LIOINC) -I$(LIO_HOME)/lioamber
>>>>>>> Debugged. Qumvia now compiles. Dont know if it works though
=======
#LIOFLAGS= -shared -I$(LIOINC) -I$(LIO_HOME)/lioamber
LIOFLAGS= -I$(LIOINC) -I$(LIO_HOME)/lioamber
>>>>>>> Fixed bug in Makefile: ditched -shared flag
=======
FFLAGS= -ip -xHost -FR -I$(MKL_HOME)/include 
=======
FFLAGS= -ip -xHost -FR -I$(MKL_HOME)/include
>>>>>>> Elimine el clean obligatorio y cositas menores en main
FLIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L$(MKL_HOME)/lib/intel64
LIOLIBS=-L/usr/lib -L/usr/lib64 -L$(LIO_HOME)/g2g -lg2g -L$(LIO_HOME)/lioamber -llio-g2g
LIOINC=$(LIO_HOME)/lioamber/obj
#LIOFLAGS= -shared -I$(LIOINC) -I$(LIO_HOME)/lioamber
#LIOFLAGS= -I$(LIOINC) -I$(LIO_HOME)/lioamber
LIOFLAGS= -I$(LIO_HOME)/lioamber -I$(LIO_HOME)/lioamber/obj
>>>>>>> Added feature: Resonant Raman intensities

cpu: clean mkdirs qumvia.cpu mvdirs

lio: mkdirs qumvia.lio mvdirs

mkdirs:
	if [ ! -d "$(BINDIR)" ]; then mkdir $(BINDIR); fi
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

qumvia.cpu: M_kracken.o qvamod_common.o qvamod_cpu.o qvamain_cpu.o
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS)

qumvia.lio: M_kracken.o qvamod_common.o qvamod_lioexcl.o qvamod_lio.o qvamain_lio.o
<<<<<<< 8f89ea86e25185c5d553f4bcb7a1cc626cb9fd1b
<<<<<<< cb52330281e84ee344e85760f1b200e1508f82b8
	$(FC) $(FFLAGS) -o $@ $^ $(FLIBS) $(LIOFLAGS) $(LIOLIBS)
=======
	$(FC) $(FFLAGS) $(LIOFLAGS) -o $(LIO_HOME)/lioamber/obj/mods/garcha_mod.o $@ $^ $(FLIBS) $(LIOLIBS)
>>>>>>> Added Resonant Raman Intesities code to qvamod_lioexcl
=======
	$(FC) $(FFLAGS) $(LIOFLAGS) -o $@ $^ $(LIOINC)/garcha_mod.o $(FLIBS) $(LIOLIBS)
<<<<<<< 3d02e19adaf35b093b863312b609077861e2d1cd
>>>>>>> Debugged. Qumvia now compiles. Dont know if it works though
=======
#	$(FC) $(FFLAGS) $(LIOFLAGS) -o $@ $^ $(FLIBS) $(LIOLIBS)
>>>>>>> Added feature: Resonant Raman intensities

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
<<<<<<< 8f89ea86e25185c5d553f4bcb7a1cc626cb9fd1b
	$(FC) -fpp -Dqvalio $(FFLAGS) $(LIOFLAGS) -c qvamod_lioexcl.f90 $(FLIBS) $(LIOLIBS) -o qvamod_lioexcl.o
=======
	$(FC) -fpp -Dqvalio $(FFLAGS) $(LIOFLAGS) -c qvamod_lioexcl.f90 $(FLIBS) \
<<<<<<< 4e69ba7e2d5240265cd90d98ec9b794eb73bd52a
<<<<<<< 3d02e19adaf35b093b863312b609077861e2d1cd
		$(LIOLIBS) -o qvamod_lioexcl.o $(LIOINC)/garcha_mod.o 
>>>>>>> Debugged. Qumvia now compiles. Dont know if it works though
=======
		$(LIOLIBS) -o qvamod_lioexcl.o 
#		$(LIOLIBS) -o qvamod_lioexcl.o $(LIOINC)/garcha_mod.o 
>>>>>>> Added feature: Resonant Raman intensities
=======
		$(LIOLIBS) -o qvamod_lioexcl.o
>>>>>>> Elimine el clean obligatorio y cositas menores en main

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
