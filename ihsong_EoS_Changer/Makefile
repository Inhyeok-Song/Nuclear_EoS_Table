SOURCES=eos_converter.cpp readtable.cpp linterp_some.cpp \
	cubinterp_some.cpp findtemp.cpp linterp_for_temp.cpp \
	linterp_for_entr.cpp linterp_for_prss.cpp nuc_eos.cpp	set_bins.cpp \
	set_bins2.cpp write_eos_table.cpp
INCLUDES=nuc_eos.h
HDF5DIR=/cluster/software/hdf5-parallel/1.8.21/gcc--8.3.0/openmpi--3.1.4
HDF5INCS=-I$(HDF5DIR)/include
MPIINCS=-I/cluster/software/openmpi/3.1.4/gcc--8.3.0/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lhdf5_hl
OBJECTS=$(SOURCES:.cpp=.o )

CXX=mpicxx
CFLAGS=-g -O3 -std=c++03
EXTRALIBS=-lm

driver: $(OBJECTS) $(INCLUDES)
	$(CXX) $(CFLAGS) -o eos_converter $(OBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.cpp $(INCLUDES)
	$(CXX) $(CFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@

clean:
	rm -f *.o eos_converter
