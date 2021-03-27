include ../Makefile.def
-include Makefile.DEPEND

SEARCH_C = -I../include -I${SHAREDIR} ${AMREXSEARCH}

# source files
SRCS := \
	SWMFInterface.cpp

OBJS := $(patsubst %, %.o,$(basename $(SRCS)))
OBJECTS = ${OBJS} PC_wrapper.o

include ../Makefile.conf

PC_wrapper.o: ../include/SWMFInterface.h \
	SWMFInterface.o \
	${LIBRARYDIR}/CON_comp_info.o

LIBFLEKS = ../src/libFLEKS.a

LIBPC    = ${LIBDIR}/libPC.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJS} ${OBJECTS_EXE}

LIB:    ${LIBFLEKS} DEPEND
	make ${LIBPC}
	@echo
	@echo ${LIBPC} has been brought up to date.
	@echo

${LIBPC}: ${OBJECTS} ${LIBFLEKS}
	cp -f ${LIBFLEKS} ${LIBPC}
	${AR} ${LIBPC} ${OBJECTS}

clean:
	-@(rm -f ${OBJECTS} ${LIBPC})
	-@(rm -rf .d)

distclean: clean
