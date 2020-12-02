# Makefile for driver using C and CUDA

############################################
# CUDA directory:
CUDAPATH=/cluster/software/cuda-10.0
############################################

############################################
# C Compiler and Options
MPIPATH=/cluster/software/openmpi/3.1.4/gcc--8.3.0
MPIINCS=-I$(MPIPATH)/include
CXX=$(MPIPATH)/bin/mpicxx
CFLAGS=-g -O3 -std=c++03
############################################

############################################
# NVCC Compiler Options
NVCC=$(CUDAPATH)/bin/nvcc
NVCCFLAGS=-ccbin /usr/bin/g++ -c
NVCCFLAGS+=-gencode arch=compute_75,code=\"compute_75,sm_75\"

# CUDA library directory:
CUDA_LIB_DIR=-L$(CUDAPATH)/lib64
# CUDA include directory:
CUDA_INC_DIR=-I$(CUDAPATH)/include
# CUDA linking libraries:
CUDA_LINK_LIBS=-lcudart
############################################


CPU_FILE=$(shell find . -type f -name "*.cpp")
GPU_FILE=$(shell find . -type f -name "*.cu")
INCLUDES=NuclearEos.h
INCLUDES_CUDA=NuclearEos.cuh
OBJECTS=$(CPU_FILE:.cpp=.o)
OBJECTS+=$(GPU_FILE:.cu=.o)


HDF5DIR =/cluster/software/hdf5-parallel/1.8.21/gcc--8.3.0/openmpi--3.1.4
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lhdf5_hl

EXTRALIBS=-lm

# Target executable name:
EXE=driver


# Compile .cpp files to objcet files:
%.o: %.cpp $(INCLUDES) $(INCLUDES_CUDA)
	$(CXX) $(CFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@

# Compile CUDA source files to objects files:
%.o: %.cu $(INCLUDES) $(INCLUDES_CUDA)
	$(NVCC) $(NVCCFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@ 

# Link C++ and CUDA compiled object files to targer executable:
$(EXE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXE) $(HDF5LIBS) $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)



clean: 
	rm -f *.o *.cuo driver *.txt
