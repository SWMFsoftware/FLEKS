include Makefile.conf
include Makefile.depend

EXE = ${INSTALL_DIR}/bin/mypic.exe

default: $(EXE)


$(EXE): $(OBJS)
	$(CXX) -o ${EXE} $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	$(CXX) -o $@ -c ${DEPFLAGS} $(FLAG_CXX) $<

.PHONY: clean distclean
clean:
	$(RM) *.o *.d
	$(RM) -rf .o .d

distclean: clean
	$(RM) ${EXE}




