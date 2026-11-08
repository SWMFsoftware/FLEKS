-include ../Makefile.def
include ../Makefile.conf

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
VERBOSE_GMRES = 
FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_ ${VERBOSE_GMRES}

# GPU build: extra nvcc warning/error flags and device linking when Makefile.conf
# uses share/build/Makefile.nvcc (./Config.pl -cuda).
ifeq (${C_COMPILER_NAME},nvcc)
  CUDA_HOME ?= $(patsubst %/bin/nvcc,%,$(shell which nvcc))
  CUDA_LIBS = -L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs -lcudart
  FLAGCC += -Wno-psabi --Werror cross-execution-space-call \
	--Werror ext-lambda-captures-this --display-error-number \
	--diag-error 20092 -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
  DLINK_OBJ = dlink.o
endif

-include Makefile.DEPEND

LIBFLEKS = libFLEKS.a

# source files
#
# Everything in this directory and in ic/ belongs to libFLEKS.a, except the two
# translation units that carry their own main() and are linked separately:
# main.cpp into ../bin/FLEKS.exe, Converter.cpp into ../bin/converter.exe.
# $(sort) fixes the archive member order: $(wildcard) returns directory order,
# which is not reproducible across machines.
SRCS := $(sort $(filter-out main.cpp Converter.cpp, \
                 $(wildcard *.cpp) $(wildcard ic/*.cpp)))

OBJECTS := $(patsubst %, %.o,$(basename $(SRCS)))

# The AMReX/IDL -> VTK/TEC format converter is no longer part of this library;
# it is built from ../Converter with its own Makefile.

default: EXE

EXE:
	$(MAKE) LIB
	$(MAKE) main.o
ifeq (${C_COMPILER_NAME},nvcc)
	$(MAKE) dlink.o
endif
	mkdir -p ../bin
	${LINK.f90} -o ../bin/FLEKS.exe main.o ${DLINK_OBJ} \
		-L. -lFLEKS -L${LIBDIR} -lSHARE ${LflagMpi} ${CUDA_LIBS} ${NOMAINFLAG}

# Device link of all relocatable device code (FLEKS and AMReX).
dlink.o: main.o ${LIBFLEKS}
	${LINK.cpp} -dlink main.o ${LIBFLEKS} ${AMREXLIB} -o $@

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
	-@(rm -f ${OBJECTS} main.o dlink.o ${LIBFLEKS})
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
