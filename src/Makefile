-include ../Makefile.def
include ../Makefile.conf

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}
VERBOSE_GMRES = 
FLAGC_EXTRA = -D_${COMPONENT}_COMPONENT_ ${VERBOSE_GMRES}

# GPU build: on when the linked AMReX library was built with CUDA
# (./Config.pl -cuda). Every source is then compiled by nvcc as CUDA with
# relocatable device code, and the device code is linked by nvcc -dlink before
# the host link. The target architectures and nvcc's host compiler are taken
# from the AMReX build so that the two always match.
AMREX_INSTALL = ${UTILDIR}/AMREX/InstallDir
USE_CUDA := $(if $(shell grep -s '^\#define AMREX_USE_CUDA 1' \
               ${AMREX_INSTALL}/include/AMReX_Config.H),TRUE,FALSE)

ifeq (${USE_CUDA},TRUE)
  AMREX_PC := ${AMREX_INSTALL}/lib/pkgconfig/amrex.pc
  NVCC_ARCH_FLAGS := $(shell grep -o -- '--generate-code [^ ]*' ${AMREX_PC})
  NVCC_CCBIN := $(patsubst -ccbin=%,%, \
                  $(firstword $(shell grep -o -- '-ccbin=[^ ]*' ${AMREX_PC})))
  CUDA_HOME ?= $(patsubst %/bin/nvcc,%,$(shell which nvcc))
  CUDA_LIBS := -L${CUDA_HOME}/lib64 -lcudart $(sort $(filter -lcu%, \
                 $(shell grep '^Libs:' ${AMREX_PC})))
  NVCC = nvcc -ccbin=${NVCC_CCBIN}
  COMPILE.mpicxx = ${NVCC} $(shell mpicxx --showme:compile)
  FLAGCC += -x cu --expt-relaxed-constexpr --expt-extended-lambda \
	--relocatable-device-code=true ${NVCC_ARCH_FLAGS} \
	--forward-unknown-to-host-compiler -Wno-deprecated-gpu-targets \
	--Werror cross-execution-space-call --Werror ext-lambda-captures-this \
	--display-error-number --diag-error 20092 \
	-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
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
ifeq (${USE_CUDA},TRUE)
	$(MAKE) dlink.o
endif
	mkdir -p ../bin
	${LINK.f90} -o ../bin/FLEKS.exe main.o ${DLINK_OBJ} \
		-L. -lFLEKS -L${LIBDIR} -lSHARE ${LflagMpi} ${CUDA_LIBS} ${NOMAINFLAG}

# Device link of all relocatable device code (FLEKS and AMReX).
dlink.o: main.o ${LIBFLEKS}
	${NVCC} -dlink ${NVCC_ARCH_FLAGS} -Wno-deprecated-gpu-targets \
		main.o ${LIBFLEKS} -L${AMREX_INSTALL}/lib -lamrex -o $@

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
