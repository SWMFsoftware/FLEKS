-include ../Makefile.def
include ../Makefile.conf

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
VERBOSE_GMRES = 
FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_ ${VERBOSE_GMRES}

# Detect if AMReX is built with CUDA support
AMREX_CONFIG_H = ${UTILDIR}/AMREX/InstallDir/include/AMReX_Config.H
AMREX_HAS_CUDA := $(shell grep -E '^\#define\s+AMREX_USE_CUDA\s+1' $(AMREX_CONFIG_H) 2>/dev/null)

ifneq ($(AMREX_HAS_CUDA),)
  CUDA_ARCH ?= 89
  COMPILE.mpicxx := OMPI_CXX=nvcc mpicxx
  FLAGCC := $(FLAGCC) --forward-unknown-to-host-compiler -ccbin=g++ --expt-relaxed-constexpr --expt-extended-lambda --relocatable-device-code=true -gencode arch=compute_$(CUDA_ARCH),code=[sm_$(CUDA_ARCH),compute_$(CUDA_ARCH)] -x cu
  CUDALIB = -lcudart -lcuda -lcurand -lcusparse
  DLINK_OBJ = dlink.o
  DLINK_CONVERTER_OBJ = dlink_converter.o
endif

-include Makefile.DEPEND

LIBFLEKS = libFLEKS.a

# source files
SRCS := \
	Domain.cpp \
	ParticleTracker.cpp \
	Particles.cpp \
	ParticlesInit.cpp \
	ParticlesBC.cpp \
	ParticlesMoments.cpp \
	ParticlesMassMatrix.cpp \
	ParticlesMover.cpp \
	ParticlesResample.cpp \
	ParticlesReactions.cpp \
	TestParticles.cpp \
	LinearSolver.cpp \
	GridUtility.cpp \
	PlotWriter.cpp \
	TimeCtr.cpp \
	BC.cpp \
	WaveBC.cpp \
	FluidInterface.cpp \
	Pic.cpp \
	PicParam.cpp \
	PicBC.cpp \
	PicHybrid.cpp \
	PicFieldSolver.cpp \
	PicDivE.cpp \
	PicIO.cpp \
	Grid.cpp  \
	DataContainer.cpp \
	VisitWriter.cpp  \
	FleksDistributionMap.cpp \
	InitialCondition.cpp \
	ic/WaveIC.cpp \
	ic/BeamIC.cpp \
	ic/TopHatIC.cpp \
	ic/FadeevIC.cpp \
	ic/RegisterAll.cpp

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

# Keep converter independent of libFLEKS.a: the full archive includes PT/OH
# coupling objects that require SWMF-only OH symbols.
CONVERTER_OBJECTS = \
	Converter.o \
	DataContainer.o \
	Grid.o \
	VisitWriter.o \
	FleksDistributionMap.o

CONVERTEREXE = ../bin/converter.exe

default: EXE

EXE:
	$(MAKE) LIB
	$(MAKE) main.o
ifneq ($(AMREX_HAS_CUDA),)
	$(MAKE) dlink.o
endif
	mkdir -p ../bin
	${LINK.f90} -o ../bin/FLEKS.exe main.o ${DLINK_OBJ} \
		-L. -lFLEKS -L${LIBDIR} -lSHARE ${LflagMpi} ${CUDALIB} ${NOMAINFLAG}

ifneq ($(AMREX_HAS_CUDA),)
dlink.o: ${OBJECTS} main.o
	nvcc -dlink -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) main.o ${OBJECTS} -L${UTILDIR}/AMREX/InstallDir/lib -lamrex -o dlink.o

dlink_converter.o: ${CONVERTER_OBJECTS}
	nvcc -dlink -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) ${CONVERTER_OBJECTS} -L${UTILDIR}/AMREX/InstallDir/lib -lamrex -o dlink_converter.o
endif

CONVERTER:	${CONVERTER_OBJECTS}
ifneq ($(AMREX_HAS_CUDA),)
	$(MAKE) dlink_converter.o
endif
	${LINK.f90} -o ${CONVERTEREXE} ${CONVERTER_OBJECTS} ${DLINK_CONVERTER_OBJ} \
		-L${LIBDIR} -lSHARE ${LflagMpi} ${CUDALIB} ${NOMAINFLAG}

	@echo
	@echo ${CONVERTEREXE} has been brought up to date.
	@echo

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} main.o ${CONVERTER_OBJECTS}

LIB: DEPEND
	$(MAKE) ${LIBFLEKS}
	@echo
	@echo ${LIBFLEKS} has been brought up to date.
	@echo

${LIBFLEKS}:${OBJECTS}	
	${AR} ${LIBFLEKS} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} main.o ${LIBFLEKS} dlink.o dlink_converter.o)
	-@(rm -rf ${CONVERTEREXE} ${CONVERTER_OBJECTS})
	-@(rm -rf .d)

distclean: clean


#FLAGCC_EXCEPTION =

#ifeq (${C_COMPILER_NAME},icc)
#	FLAGCC_EXCEPTION = -no-ip -no-vec -fno-inline
#endif

# OH-PT simulations may fail with -O2 but without -no-ip -no-vec -fno-inline
# on Pleiades (icc). It seems the failure is caused by the split_particles()
# method. The crash can not be re-produced with -O1 or -O0 or with gcc. It
# can be an icc bug.
#Particles.o: Particles.cpp
#	${COMPILE.mpicxx} ${FLAGCC} ${FLAGCC_EXCEPTION}  Particles.cpp
