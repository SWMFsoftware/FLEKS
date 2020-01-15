#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR} 

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

#FLAGC_EXTRA = -DBATSRUS

include ../Makefile.conf
include Makefile.DEPEND

# OBJECT FILE FOR THE WRAPPER

OBJECTS = ${OBJS} PC_wrapper.o

PC_wrapper.o: ../include/SWMFInterface.h \
	SWMFInterface.cpp \
	${LIBRARYDIR}/CON_comp_info.o

LIBFLEKS = ../libFLEKS/libFLEKS.a

MY_LIB    = ${LIBDIR}/libPC.a

LIB:    ${LIBFLEKS} 
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} ${LIBFLEKS}
	cp -f ${LIBFLEKS} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}


%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	${COMPILE.mpicxx} -o $@ -c ${DEPFLAGS} $(FLAGCC) $<


clean:
	-@(rm -f ${OBJECTS} ${MY_LIB})
	-@(rm -rf .d)


distclean: clean
