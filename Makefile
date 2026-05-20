default: EXE

-include Makefile.conf
-include Makefile.def

help:
	@echo Makefile targets:
	@echo
	@echo 'make EXE                    - compile standalone executable'
	@echo 'make LIB                    - compile libPC.a for SWMF'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

include/show_git_info.h: include/show_git_info.h.orig
	cp -f include/show_git_info.h.orig include/show_git_info.h

include/Constants.h: include/Constants.h.orig
	cp -f include/Constants.h.orig include/Constants.h 

include/UserSource.h: include/UserSource.h.orig
	cp -f include/UserSource.h.orig include/UserSource.h

GITINFO: include/show_git_info.h
	${SCRIPTDIR}/gitall -r=c > include/show_git_info.h

# Detect whether FLEKS lives inside the SWMF tree or is a pure standalone checkout.
# If ../../share/Library/src exists we are inside SWMF; otherwise pure standalone.
_SWMF_SHARE := $(abspath $(CURDIR)/../../share)
ifneq ($(wildcard $(_SWMF_SHARE)/Library/src),)
  # Inside SWMF tree: reuse SWMF's share and lib directories.
  SHARE_ROOT = $(_SWMF_SHARE)
  SHARE_SRC  = $(SHARE_ROOT)/Library/src
  SHARE_LIB  = $(abspath $(CURDIR)/../../lib/libSHARE.a)
else
  # Pure standalone: use share/ cloned under the FLEKS directory.
  SHARE_ROOT = $(abspath $(CURDIR)/share)
  SHARE_SRC  = $(SHARE_ROOT)/Library/src
  SHARE_LIB  = $(abspath $(CURDIR)/lib/libSHARE.a)
endif

GITCLONE_SHARE  = git clone git@github.com:SWMFsoftware/share $(CURDIR)/share
BUILD_MODE_FILE = src/.build_mode

EXE: GITINFO compile_commands
	@if [ ! -d $(SHARE_ROOT)/Library/src ]; then \
		echo "--- Cloning SWMFsoftware/share (not found at $(SHARE_ROOT)) ---"; \
		$(GITCLONE_SHARE); \
	fi
	@if [ ! -f Makefile.conf ]; then \
		echo "--- Configuring FLEKS ---"; \
		./Config.pl -install; \
	fi
	@if [ ! -f $(SHARE_LIB) ]; then \
		echo "--- Building libSHARE.a ---"; \
		$(MAKE) -C $(SHARE_SRC) LIB; \
	fi
	@if [ -f $(BUILD_MODE_FILE) ] && [ "$$(cat $(BUILD_MODE_FILE))" != "STANDALONE" ]; then \
		echo "--- Switching from component to standalone mode: auto-cleaning src/ ---"; \
		cd src; $(MAKE) clean; \
	fi
	@echo STANDALONE > $(BUILD_MODE_FILE)
	cd src; $(MAKE) EXE STANDALONE=YES FLEKS_DIR=$(CURDIR)

bin:
	mkdir bin

install: bin include/Constants.h include/UserSource.h

LIB: include/Constants.h compile_commands
	@if [ ! -d $(_SWMF_SHARE)/Library/src ]; then \
		echo "--- Building FLEKS library in standalone mode ---"; \
		if [ -f $(BUILD_MODE_FILE) ] && [ "$$(cat $(BUILD_MODE_FILE))" != "STANDALONE" ]; then \
			echo "--- Switching from component to standalone mode: auto-cleaning src/ ---"; \
			cd src; $(MAKE) clean; \
		fi; \
		echo STANDALONE > $(BUILD_MODE_FILE); \
		cd src; $(MAKE) LIB STANDALONE=YES FLEKS_DIR=$(CURDIR); \
	else \
		if [ ! -f $(SHARE_LIB) ]; then \
			echo ""; \
			echo "ERROR: $(SHARE_LIB) not found."; \
			echo "       The SWMF share library must be built before 'make LIB'."; \
			echo "       Please run the top-level SWMF build first."; \
			echo ""; \
			exit 1; \
		fi; \
		if [ ! -f $(_SWMF_SHARE)/include/con_comp_param.mod ]; then \
			echo ""; \
			echo "ERROR: CON Fortran modules not found in $(_SWMF_SHARE)/include/."; \
			echo "       The SWMF CON library must be built before 'make LIB'."; \
			echo "       Please run the top-level SWMF build first."; \
			echo ""; \
			exit 1; \
		fi; \
		if [ -f $(BUILD_MODE_FILE) ] && [ "$$(cat $(BUILD_MODE_FILE))" != "COMPONENT" ]; then \
			echo "--- Switching from standalone to component mode: auto-cleaning src/ ---"; \
			cd src; $(MAKE) clean; \
		fi; \
		echo COMPONENT > $(BUILD_MODE_FILE); \
		cd src; $(MAKE) LIB; \
		cd srcInterface; $(MAKE) LIB; \
	fi

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
	rm -f $(BUILD_MODE_FILE)

distclean:	
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf *~ ./bin lib ${TESTDIR} include/Constants.h)
	-@(rm -f test*.diff)

compile_commands:
	-rm -f compile_commands.json
	-python3 tools/generate_compile_commands.py

