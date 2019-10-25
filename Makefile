default: mypic

include Makefile.conf
include Makefile.def

help:
	@echo Makefile targets:
	@echo
	@echo 'make                        - compile IPIC3D.exe'
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make PIDL                   - compile PostIDL.exe for post processing'
	@echo 'make run                    - create run directory'
	@echo 'make test                   - test IPIC3D in stand alone mode'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo


mypic:
	cd src; make exe

INSTALLFILES =  src/Makefile.DEPEND \
		srcInterface/Makefile.DEPEND

bin:
	mkdir bin

libIPIC3D:
	mkdir libIPIC3D


install: bin libIPIC3D
	touch ${INSTALLFILES}

# Switch to "coupled" mode
LIB:
	touch ${INSTALLFILES}
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB


clean:
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean

distclean:
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf Makefile.local *~ ./bin libIPIC3D ${TESTDIR})
	-@(rm test*diff )
