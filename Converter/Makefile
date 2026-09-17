# Must come first: ../Makefile.conf defines targets of its own (cleanfiles),
# and the first target of the file is make's default goal.
default: CONVERTER

-include ../Makefile.def
include ../Makefile.conf

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_

-include Makefile.DEPEND

CONVERTEREXE = ../bin/converter.exe

# Sources of the standalone AMReX/IDL -> VTK/TEC format converter. It is kept
# out of src/ on purpose: none of it is part of the FLEKS library, and it has
# its own main().
SRCS := \
	Converter.cpp \
	DataContainer.cpp \
	VisitWriter.cpp

# DataContainer::AMReXDataContainer derives from Grid, so the converter needs
# these two FLEKS sources. They are compiled here rather than taken from
# libFLEKS.a: the full archive also holds the PT/OH coupling objects, which
# reference OH symbols that only exist inside SWMF.
SHARED_SRCS := \
	../src/Grid.cpp \
	../src/FleksDistributionMap.cpp

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))
SHARED_OBJECTS := $(notdir $(patsubst %, %.o,$(basename $(SHARED_SRCS))))

default: CONVERTER

CONVERTER: ${CONVERTEREXE}

${CONVERTEREXE}: ${OBJECTS} ${SHARED_OBJECTS}
	mkdir -p ../bin
	${LINK.f90} -o ${CONVERTEREXE} ${OBJECTS} ${SHARED_OBJECTS} \
		-L${LIBDIR} -lSHARE ${LflagMpi} ${NOMAINFLAG}
	@echo
	@echo ${CONVERTEREXE} has been brought up to date.
	@echo

# The shared sources live in ../src; build their objects here so that they do
# not collide with the ones libFLEKS.a is built from.
Grid.o: ../src/Grid.cpp
	${COMPILE.mpicxx} ${FLAGCC} $< -o $@

FleksDistributionMap.o: ../src/FleksDistributionMap.cpp
	${COMPILE.mpicxx} ${FLAGCC} $< -o $@

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} ${SHARED_OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${SHARED_OBJECTS} ${CONVERTEREXE})
	-@(rm -rf .d)

distclean: clean
