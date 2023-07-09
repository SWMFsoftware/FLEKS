include ../Makefile.def
-include Makefile.DEPEND

# source files
SRCS := \
	Domain.cpp \
	ParticleTracker.cpp \
	Particles.cpp \
	TestParticles.cpp \
	LinearSolver.cpp \
	Utility.cpp \
	PlotWriter.cpp \
	TimeCtr.cpp \
	BC.cpp \
	FluidInterface.cpp \
	Pic.cpp \
	PicIO.cpp \
	Grid.cpp

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

OBJECTS_EXE = \
	main.o

FLEKSEXE = ${BINDIR}/FLEKS.exe

CONVERTER_OBJECTS = \
	Converter.o

CONVERTEREXE = ../bin/converter.exe

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_
# FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_ -D_AMR_DEV_


# If hdf5.h exists, assume AMReX is built with HDF5 support.
ifneq ("$(wildcard $(UTILDIR)/AMREX/InstallDir/include/hdf5.h)","")
    FLAGC_EXTRA := ${FLAGC_EXTRA} -D_USE_HDF5_
endif

LIBFLEKS  = libFLEKS.a

include ../Makefile.conf

default: EXE

EXE:
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	${LINK.f90} -o ${FLEKSEXE} ${OBJECTS_EXE} \
		-L -lFLEKS ${LflagMpi}

CONVERTER:	
	$(MAKE) ${CONVERTER_OBJECTS}
	${LINK.f90} ${NOMAINFLAG} -o ${CONVERTEREXE} ${CONVERTER_OBJECTS} ${LflagMpi}
	@echo
	@echo ${CONVERTEREXE} has been brought up to date.
	@echo

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} ${OBJECTS_EXE}

LIB: DEPEND
	$(MAKE) ${LIBFLEKS}
	@echo
	@echo ${LIBFLEKS} has been brought up to date.
	@echo
	$(MAKE) CONVERTER

${LIBFLEKS}:${OBJECTS}	
	${AR} ${LIBFLEKS} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${OBJECTS_EXE} ${LIBFLEKS})
	-@(rm -rf ${CONVERTEREXE} ${CONVERTER_OBJECTS})
	-@(rm -rf .d)

distclean: clean


FLAGCC_EXCEPTION =

ifeq (${C_COMPILER_NAME},icc)
	FLAGCC_EXCEPTION = -no-ip -no-vec -fno-inline
endif

# OH-PT simulations may fail with -O2 but without -no-ip -no-vec -fno-inline
# on Pleiades (icc). It seems the failure is caused by the split_particles()
# method. The crash can not be re-produced with -O1 or -O0 or with gcc. It
# can be an icc bug.
Particles.o: Particles.cpp
	${COMPILE.mpicxx} ${FLAGCC} ${FLAGCC_EXCEPTION}  Particles.cpp
