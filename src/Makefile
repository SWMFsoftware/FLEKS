FLEKS_DIR ?= $(abspath ..)

ifeq ($(STANDALONE),)
  STANDALONE := $(if $(wildcard $(FLEKS_DIR)/../../share/Library/src),,YES)
endif

ifeq ($(STANDALONE), YES)
  include ../Makefile.conf
  # The parent Makefile overrides these when FLEKS is inside an SWMF tree.
  SHARE_ROOT ?= $(FLEKS_DIR)/share
  UTIL_ROOT  ?= $(FLEKS_DIR)/util
  LIB_ROOT   ?= $(FLEKS_DIR)/lib
  SHARE_DIR ?= $(SHARE_ROOT)/Library/src
  AMREX_DIR ?= $(UTIL_ROOT)/AMREX/InstallDir
  FLEKS_SCRIPTDIR ?= $(SHARE_ROOT)/Scripts
  LDFS ?= -L$(LIB_ROOT) -lSHARE
  LDFLAGS ?= -L$(AMREX_DIR)/lib -lamrex -lgfortran -lm
  SEARCH_C = -I../include -I${SHARE_DIR} -I${AMREX_DIR}/include
  FLAGC_EXTRA = -D_PC_COMPONENT_ -DFLEKS_STANDALONE
  BINDIR = ../bin
  ifeq ($(origin CXX), default)
    CXX = mpicxx
  endif
  CXXFLAGS ?= -O3 -std=c++17
  COMPILE.mpicxx = $(CXX) $(CXXFLAGS) $(SEARCH_C) $(FLAGC_EXTRA)


  %.o: %.cpp
	$(COMPILE.mpicxx) -c $< -o $@
else
  include ../Makefile.def
  -include Makefile.DEPEND
  include ../Makefile.conf
  SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
  FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_
  LDFS ?= -L${LIBDIR} -lSHARE
endif

LIBFLEKS = libFLEKS.a

# source files
SRCS := \
	Domain.cpp \
	ParticleTracker.cpp \
	Particles.cpp \
	TestParticles.cpp \
	LinearSolver.cpp \
	GridUtility.cpp \
	PlotWriter.cpp \
	TimeCtr.cpp \
	BC.cpp \
	FluidInterface.cpp \
	Pic.cpp \
	PicIO.cpp \
	Grid.cpp  \
	DataContainer.cpp \
	VisitWriter.cpp  \
	FleksDistributionMap.cpp

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

OBJECTS_EXE = \
	main.o

FLEKSEXE = ${BINDIR}/FLEKS.exe

CONVERTER_OBJECTS = Converter.o

CONVERTEREXE = ../bin/converter.exe

default: EXE

EXE:
ifeq ($(STANDALONE), YES)
	$(MAKE) LIB STANDALONE=YES
	$(MAKE) main.o STANDALONE=YES
	mkdir -p ${BINDIR}
	$(COMPILE.mpicxx) -o ${FLEKSEXE} main.o \
		-L. -lFLEKS ${LDFS} $(LDFLAGS)
else
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	${LINK.f90} -o ${FLEKSEXE} ${OBJECTS_EXE} \
		-L. -lFLEKS ${LflagMpi}
endif

CONVERTER:	${CONVERTER_OBJECTS}	
	$(MAKE) ${LIBFLEKS}
ifeq ($(STANDALONE), YES)
	mkdir -p ${BINDIR}
	$(COMPILE.mpicxx) -o ${CONVERTEREXE} ${CONVERTER_OBJECTS} \
		-L. -lFLEKS ${LDFS} $(LDFLAGS)
else
	${LINK.f90} ${NOMAINFLAG} -o ${CONVERTEREXE} ${CONVERTER_OBJECTS} \
		-L. -lFLEKS ${LDFS} ${LflagMpi}
endif

	@echo
	@echo ${CONVERTEREXE} has been brought up to date.
	@echo

DEPEND:
ifeq ($(STANDALONE), YES)
	@perl ${FLEKS_SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} main.cpp ${CONVERTER_OBJECTS}
else
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} ${OBJECTS_EXE} ${CONVERTER_OBJECTS}
endif

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
