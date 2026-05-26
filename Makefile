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

# Standalone builds can run from SWMF/PC/FLEKS or from a pure FLEKS checkout.
SWMF_ROOT  := $(abspath $(CURDIR)/../..)
SWMF_SHARE := $(SWMF_ROOT)/share
IN_SWMF_TREE := $(wildcard $(SWMF_SHARE)/Library/src)

ifneq ($(IN_SWMF_TREE),)
  SHARE_ROOT = $(SWMF_SHARE)
  UTIL_ROOT  = $(SWMF_ROOT)/util
  LIB_ROOT   = $(SWMF_ROOT)/lib
else
  SHARE_ROOT = $(abspath $(CURDIR)/share)
  UTIL_ROOT  = $(abspath $(CURDIR)/util)
  LIB_ROOT   = $(abspath $(CURDIR)/lib)
endif

SHARE_SRC = $(SHARE_ROOT)/Library/src
SHARE_LIB = $(LIB_ROOT)/libSHARE.a

GITCLONE_SHARE  = git clone git@github.com:SWMFsoftware/share $(CURDIR)/share
BUILD_MODE_FILE = src/.build_mode
STANDALONE_SRC_MAKE = $(MAKE) -C src STANDALONE=YES FLEKS_DIR=$(CURDIR) \
	SHARE_ROOT=$(SHARE_ROOT) UTIL_ROOT=$(UTIL_ROOT) LIB_ROOT=$(LIB_ROOT)

define prepare_standalone
	@if [ ! -d $(SHARE_SRC) ]; then \
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
endef

# Source objects cannot be reused across standalone/component preprocessor flags.
define set_build_mode
	@if [ -f $(BUILD_MODE_FILE) ] && [ "$$(cat $(BUILD_MODE_FILE))" != "$(1)" ]; then \
		echo "--- Switching to $(2) mode: auto-cleaning src/ ---"; \
		$(MAKE) -C src clean; \
	fi
	@echo $(1) > $(BUILD_MODE_FILE)
endef

EXE: GITINFO compile_commands
	$(call prepare_standalone)
	$(call set_build_mode,STANDALONE,standalone)
	$(STANDALONE_SRC_MAKE) EXE

bin:
	mkdir bin

install: bin include/Constants.h include/UserSource.h

LIB: bin include/Constants.h include/UserSource.h compile_commands
ifeq ($(IN_SWMF_TREE),)
	@echo "--- Building FLEKS library in standalone mode ---"
	$(call set_build_mode,STANDALONE,standalone)
	$(STANDALONE_SRC_MAKE) LIB
else
	@if [ ! -f $(SHARE_LIB) ]; then \
		echo ""; \
		echo "ERROR: $(SHARE_LIB) not found."; \
		echo "       The SWMF share library must be built before 'make LIB'."; \
		echo "       Please run the top-level SWMF build first."; \
		echo ""; \
		exit 1; \
	fi
	@if [ ! -f $(SWMF_SHARE)/include/con_comp_param.mod ]; then \
		echo ""; \
		echo "ERROR: CON Fortran modules not found in $(SWMF_SHARE)/include/."; \
		echo "       The SWMF CON library must be built before 'make LIB'."; \
		echo "       Please run the top-level SWMF build first."; \
		echo ""; \
		exit 1; \
	fi
	$(call set_build_mode,COMPONENT,component)
	$(MAKE) -C src LIB
	$(MAKE) -C srcInterface LIB
endif

CONVERTER: bin
	$(call prepare_standalone)
	$(call set_build_mode,STANDALONE,standalone)
	$(STANDALONE_SRC_MAKE) CONVERTER

rundir:
	mkdir -p ${RUNDIR}/${COMPONENT}
	(cd ${RUNDIR}/${COMPONENT}; \
		mkdir restartIN restartOUT plots;\
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${SCRIPTDIR}/pIDL .)
	(cd ${RUNDIR}; \
		if [ -f ${BINDIR}/FLEKS.exe ]; then ln -s ${BINDIR}/FLEKS.exe .; fi)


clean:
	$(MAKE) -C src clean
	$(MAKE) -C srcInterface clean
	rm -rf bin/*
	rm -f $(BUILD_MODE_FILE)

distclean:
	-@./Config.pl -uninstall

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@rm -rf *~ ./bin lib ${TESTDIR} include/Constants.h
	-@rm -f test*.diff

compile_commands:
	-rm -f compile_commands.json
	-python3 tools/generate_compile_commands.py
