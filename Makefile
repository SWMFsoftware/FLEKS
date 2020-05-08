default: FLEKS

include Makefile.conf
include Makefile.def

help:
	@echo Makefile targets:
	@echo
	#@echo 'make FLEKS                  - compile standalone executable'
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

include/show_git_info.h: include/show_git_info.h.orig
	cp -f include/show_git_info.h.orig include/show_git_info.h

GITINFO: include/show_git_info.h
	${SCRIPTDIR}/gitall -r=c > include/show_git_info.h

FLEKS: GITINFO 
	cd src; ${MAKE} EXE

bin:
	mkdir bin

install: bin

LIB: 
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB


rundir:
	mkdir -p ${RUNDIR}/PC
	cd ${RUNDIR}/PC; \
		mkdir restartIN restartOUT plots;\
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${SCRIPTDIR}/pIDL .
	cd ${RUNDIR}; \
		(if [ -f ${BINDIR}/FLEKS.exe ]; then ln -s ${BINDIR}/FLEKS.exe .; fi)


clean:
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean
	rm -rf bin/*

distclean:
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf *~ ./bin lib ${TESTDIR})
	-@(rm -f test*.diff)

