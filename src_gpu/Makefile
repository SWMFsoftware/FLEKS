default: GPUTEST

-include ../Makefile.def
include ../Makefile.conf

.PHONY: default GPUTEST GPUTEST_REFERENCE clean distclean

SEARCH_C = -I. ${AMREXSEARCH}
FLAGC_EXTRA =
GPUDEMOEXE = ../bin/main_gpu.exe
REFERENCEEXE = ../bin/gpu_reference.exe
SRCS := $(sort $(filter-out ProductionReference.cpp,$(wildcard *.cpp)))
OBJECTS := $(SRCS:.cpp=.o)
REFERENCEOBJECTS = ProductionReference.o DemoConfig.o DemoGrid.o DemoFields.o DemoLinearSolver.o

GPUTEST: ${GPUDEMOEXE}

${GPUDEMOEXE}: ${OBJECTS} Makefile ../Makefile.conf ../Makefile.def
	mkdir -p ../bin
	${LINK.cpp} -o $@ ${OBJECTS} ${AMREXLIB} ${CPPLIB}

# nvcc and the supported host compilers emit Make-compatible dependencies.
# -MP supplies harmless rules for removed headers, avoiding stale .d failures.
%.o: %.cpp Makefile ../Makefile.conf ../Makefile.def
	${COMPILE.mpicxx} ${FLAGCC} -MMD -MP -MF $(@:.o=.d) $< -o $@

GPUTEST_REFERENCE: ${REFERENCEEXE}

ProductionOperators.h: validate_demo.py ../src/GridUtility.cpp
	python3 validate_demo.py --emit-production-reference $@

ProductionReference.o: ProductionOperators.h

${REFERENCEEXE}: ${REFERENCEOBJECTS} Makefile ../Makefile.conf ../Makefile.def
	mkdir -p ../bin
	${LINK.cpp} -o $@ ${REFERENCEOBJECTS} ${AMREXLIB} ${CPPLIB}

-include $(OBJECTS:.o=.d) ProductionReference.d

clean:
	rm -f *.o *.d ${GPUDEMOEXE} ${REFERENCEEXE} ProductionOperators.h

distclean: clean
