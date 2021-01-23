#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR} 

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

#FLAGC_EXTRA = -DBATSRUS

# source files
SRCS := \
	SWMFInterface.cpp

include ../Makefile.conf
include ../Makefile.conf.local

ifeq ($(C_COMPILER_NAME), nvc)
      DEPFLAGS = 
endif

ifeq ($(C_COMPILER_NAME), pgcc)
      DEPFLAGS = 
endif

# OBJECT FILE FOR THE WRAPPER

OBJECTS = ${OBJS} PC_wrapper.o

PC_wrapper.o: ../include/SWMFInterface.h \
	SWMFInterface.cpp \
	${LIBRARYDIR}/CON_comp_info.o

LIBFLEKS = ../src/libFLEKS.a

LIBPC    = ${LIBDIR}/libPC.a

LIB:    ${LIBFLEKS} 
	make ${LIBPC}
	@echo
	@echo ${LIBPC} has been brought up to date.
	@echo

${LIBPC}: ${OBJECTS} ${LIBFLEKS}
	cp -f ${LIBFLEKS} ${LIBPC}
	${AR} ${LIBPC} ${OBJECTS}


%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	${COMPILE.mpicxx} -o $@ -c ${DEPFLAGS} $(FLAGCC) $<


clean:
	-@(rm -f ${OBJECTS} ${LIBPC})
	-@(rm -rf .d)


distclean: clean
