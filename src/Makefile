# source files
SRCS := \
	Domain.cpp \
	Particles.cpp \
	LinearSolver.cpp \
	Utility.cpp \
	PlotWriter.cpp \
	TimeCtr.cpp \
	BC.cpp \
	FluidInterface.cpp \
	Pic.cpp \
	PicIO.cpp

include ../Makefile.conf.local
include ../Makefile.def

OBJECTS = ${OBJS}

OBJECTS_EXE = \
	main.o

EXE = ${BINDIR}/FLEKS.exe

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

FLAGC_EXTRA = ${FLAGC_LOCAL} 

include ../Makefile.conf

EXE: $(EXE)

default: $(EXE)

$(EXE):
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} \
		-L../lib/ -lFLEKS ${Lflag}

#

%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	${COMPILE.mpicxx} -o $@ -c ${DEPFLAGS} $(FLAGCC) $<

DEPEND:


LIBFLEKS  = ../lib/libFLEKS.a

LIB:    
	$(MAKE) ${LIBFLEKS}
	@echo
	@echo ${LIBFLEKS} has been brought up to date.
	@echo

${LIBFLEKS}: ${OBJECTS} 
	${AR} ${LIBFLEKS} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${OBJECTS_EXE} ${LIBFLEKS})
	-@(rm -rf .d)

distclean: clean
