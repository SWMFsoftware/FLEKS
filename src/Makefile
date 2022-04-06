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
	PicIO.cpp	

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

OBJECTS_EXE = \
	main.o

FLEKSEXE = ${BINDIR}/FLEKS.exe

CONVERTER_OBJECTS = \
	Converter.o

CONVERTEREXE = ../bin/converter.exe

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

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

LIB:    DEPEND
	$(MAKE) ${LIBFLEKS}
	@echo
	@echo ${LIBFLEKS} has been brought up to date.
	@echo
	$(MAKE) CONVERTER

${LIBFLEKS}: ${OBJECTS} 
	${AR} ${LIBFLEKS} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${OBJECTS_EXE} ${LIBFLEKS})
	-@(rm -rf ${CONVERTEREXE} ${CONVERTER_OBJECTS})
	-@(rm -rf .d)

distclean: clean
