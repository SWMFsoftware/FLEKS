default: mypic

include Makefile.conf
include Makefile.def

help:
	@echo Makefile targets:
	@echo
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo


mypic:
	cd src; make exe

bin:
	mkdir bin

libFLEKS:
	mkdir libFLEKS


install: bin libFLEKS

LIB: 
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB


rundir:
	mkdir -p ${RUNDIR}/PC
	cd ${RUNDIR}/PC; \
		mkdir restartIN restartOUT plots;\
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${SCRIPTDIR}/pIDL .


clean:
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean

distclean:
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf *~ ./bin libFLEKS ${TESTDIR})
	-@(rm -f test*.diff)

