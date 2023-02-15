default: FLEKS

-include Makefile.conf
-include Makefile.def

help:
	@echo Makefile targets:
	@echo
	#@echo 'make FLEKS                  - compile standalone executable'
	@echo 'make LIB                    - compile libPC.a for SWMF'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

include/show_git_info.h: include/show_git_info.h.orig
	cp -f include/show_git_info.h.orig include/show_git_info.h

include/Constants.h: include/Constants.h.orig
	cp -f include/Constants.h.orig include/Constants.h

GITINFO: include/show_git_info.h
	${SCRIPTDIR}/gitall -r=c > include/show_git_info.h

FLEKS: GITINFO 
	cd src; ${MAKE} EXE

bin:
	mkdir bin

install: bin include/Constants.h

LIB: include/Constants.h srcBATL/Makefile
	cd srcBATL; $(MAKE) LIB
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

srcBATL/Makefile:
	rm -rf srcBATL; mkdir srcBATL 
	cd srcBATL_orig; cp BATL*.f90 Makefile* ../srcBATL; \
	cd ../srcBATL; ${SCRIPTDIR}/Methods.pl FL *.f90; \
	${SCRIPTDIR}/Rename.pl -w -r -common=FL *.f90; \
	rm -f *~
	cp -f srcBATL/BATL_size_orig.f90 srcBATL/BATL_size.f90

CONVERTER:
	cd src; $(MAKE) CONVERTER

rundir:
	mkdir -p ${RUNDIR}/${COMPONENT}
	cd ${RUNDIR}/${COMPONENT}; \
		mkdir restartIN restartOUT plots;\
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${SCRIPTDIR}/pIDL .
	cd ${RUNDIR}; \
		(if [ -f ${BINDIR}/FLEKS.exe ]; then ln -s ${BINDIR}/FLEKS.exe .; fi)


clean:
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean
	rm -rf bin/*
	rm -rf srcBATL

distclean:
	rm -rf srcBATL*
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf *~ ./bin lib ${TESTDIR} include/Constants.h)
	-@(rm -f test*.diff)

