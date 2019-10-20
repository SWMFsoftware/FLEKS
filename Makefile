include Makefile.def

default: mypic


mypic:
	cd src; make

clean:
	cd src; make clean

distclean:
	cd src; make distclean

rundir:
	mkdir run
	cd run; ln -s ${INSTALL_DIR}/bin/mypic.exe .
