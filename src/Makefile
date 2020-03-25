# source files
SRCS := \
	Domain.cpp \
	Particles.cpp \
	LinearSolver.cpp \
	Utility.cpp \
	PlotWriter.cpp \
	TimeCtr.cpp \
	DomainIO.cpp \
	BC.cpp \
	FluidInterface.cpp

include ../Makefile.conf.local
include ../Makefile.def

OBJECTS = ${OBJS}

OBJECTS_EXE = \
	main.o

EXE = ${BINDIR}/FLEKS.exe

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

FLAGC_EXTRA = ${FLAGC_LOCAL} 

# LDFLAGS = -L$(MYDIR)/${AMREX_INSTALL_DIR}/lib
# LIBRARIES = -lamrex  -lgfortran

include ../Makefile.conf

EXE: $(EXE)

default: $(EXE)

$(EXE):
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	${COMPILE.mpicxx} -o ${EXE} ${MY_LIB} ${OBJECTS_EXE} \
		-L../libFLEKS/ -lFLEKS ${CPPLIB} ${MPILIB} \
		$(LDFLAGS) $(LIBRARIES)


#

%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	${COMPILE.mpicxx} -o $@ -c ${DEPFLAGS} $(FLAGCC) $<

DEPEND:


MY_LIB  = ../libFLEKS/libFLEKS.a

LIB:    
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} 
	${AR} ${MY_LIB} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${OBJECTS_EXE} ${MY_LIB})
	-@(rm -rf .d)

distclean: clean
