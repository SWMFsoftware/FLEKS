include Makefile.conf

default: main.exe

main.exe: main.o 
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

# main.o: main.cpp MyParams.H
# 	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

# test_parameters.o: test_parameters.cpp MyParams.H
# 	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

# my_func.o: my_func.f90
# 	$(FC) -o $@ -c $(FFLAGS) $(CPPFLAGS) $<

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) main.exe
