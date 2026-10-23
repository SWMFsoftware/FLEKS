-include ../Makefile.def
include ../Makefile.conf

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
VERBOSE_GMRES = 
FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_ ${VERBOSE_GMRES}

-include Makefile.DEPEND

LIBFLEKS = libFLEKS.a

# source files
#
# Everything in this directory and in ic/ belongs to libFLEKS.a, except the
# translation units that carry their own main() and are linked separately:
# main.cpp into ../bin/FLEKS.exe, Converter.cpp into ../bin/converter.exe and
# main_gpu.cpp into ../bin/main_gpu.exe.
# $(sort) fixes the archive member order: $(wildcard) returns directory order,
# which is not reproducible across machines.
SRCS := $(sort $(filter-out main.cpp main_gpu.cpp Converter.cpp, \
                 $(wildcard *.cpp) $(wildcard ic/*.cpp)))

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

# The AMReX/IDL -> VTK/TEC format converter is no longer part of this library;
# it is built from ../Converter with its own Makefile.

default: EXE

EXE:
	$(MAKE) LIB
	$(MAKE) main.o
	mkdir -p ../bin
	${LINK.f90} -o ../bin/FLEKS.exe main.o \
		-L. -lFLEKS -L${LIBDIR} -lSHARE ${LflagMpi} ${NOMAINFLAG}

# AMReX GPU test program. It does not use libFLEKS.a. Configure with
# ./Config.pl -install -compiler=nvfortran,nvcc so that main_gpu.cpp is
# compiled by nvcc and linked with the CUDA AMReX library.
GPUTEST: main_gpu.o
	mkdir -p ../bin
	${LINK.cpp} -o ../bin/main_gpu.exe main_gpu.o ${AMREXLIB} ${CPPLIB}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} main.o

LIB: DEPEND
	$(MAKE) ${LIBFLEKS}
	@echo
	@echo ${LIBFLEKS} has been brought up to date.
	@echo

${LIBFLEKS}:${OBJECTS}	
	${AR} ${LIBFLEKS} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} main.o main_gpu.o ${LIBFLEKS})
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
