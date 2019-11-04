include Makefile.depend
include ../Makefile.def
include ../Makefile.def.local

OBJECTS = ${OBJS}

OBJECTS_EXE = \
	main.o

EXE = ${BINDIR}/mypic.exe

SEARCH_C = -Iinclude -I$(MYDIR)/${AMREX_INSTALL_DIR}/include -I${SHAREDIR} 

FLAGC_EXTRA = ${FLAGC_LOCAL} 

LDFLAGS = -L$(MYDIR)/${AMREX_INSTALL_DIR}/lib
LIBRARIES = -lamrex  -lgfortran

include ../Makefile.conf

exe: $(EXE)

default: $(EXE)

$(EXE):
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	${COMPILE.mpicxx} -o ${EXE} ${MY_LIB} ${OBJECTS_EXE} \
	-L../libIPIC3D/ -lIPIC3D ${CPPLIB} ${MPILIB} \
	$(LDFLAGS) $(LIBRARIES)


#

%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	${COMPILE.mpicxx} -o $@ -c ${DEPFLAGS} $(FLAGCC) $<

DEPEND:


MY_LIB  = ../libIPIC3D/libIPIC3D.a

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
	-@(rm -f Makefile.DEPEND Makefile.local Makefile.RULES)
