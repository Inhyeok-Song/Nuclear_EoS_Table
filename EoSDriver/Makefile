SOURCES=driver.cpp readtable.cpp linterp_some.cpp cubinterp_some.cpp \
	findenergy.cpp NuclearEos.cpp test_nuc_eos.cpp
INCLUDES=NuclearEos.h
HDF5DIR=/cluster/software/hdf5-parallel/1.8.21/gcc--8.3.0/openmpi--3.1.4
HDF5INCS=-I$(HDF5DIR)/include
MPIINCS=-I/cluster/software/openmpi/3.1.4/gcc--8.3.0/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5
OBJECTS=$(SOURCES:.cpp=.o )

CXX=mpicxx
CFLAGS=-g -O3 -std=c++03
EXTRALIBS=-lm

driver: $(OBJECTS) $(INCLUDES)
	$(CXX) $(CFLAGS) -o driver $(OBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.cpp $(INCLUDES)
	$(CXX) $(CFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@

clean:
	rm -f *.o driver *.txt
