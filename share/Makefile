default:libamrex

include ../Makefile.def
include ../Makefile.def.local

libdir = ${LIBDIR}

libamrex:
	cd $(MYDIR)/$(amrexdir); make all; make install; \
	cp -f $(MYDIR)/$(AMREX_INSTALL_DIR)/lib/libamrex.a $(libdir)	 

clean:
	cd $(MYDIR)/$(amrexdir); make clean

distclean:
	cd $(MYDIR)/$(amrexdir); make distclean 

