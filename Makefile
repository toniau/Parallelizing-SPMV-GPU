#/*
# *********************************************
# *  314 Principles of Programming Languages  *
# *  Fall 2016                                *
# *********************************************
# */

# Location of the CUDA Toolkit
CUDA_PATH       ?= /usr/local/cuda-7.5


HOST_COMPILER ?= gcc
#NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)
NVCC          := $(CUDA_PATH)/bin/nvcc

# internal flags
#NVCCFLAGS   :=
NVCCFLAGS   := -m64 -gencode arch=compute_30,code=sm_30
CCFLAGS     :=
LDFLAGS     :=



ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(CCFLAGS)

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS) ALL_LDFLAGS += $(LDFLAGS)

# Common includes and paths for CUDA
INCLUDES  :=
LIBRARIES :=





#CCFLAGS = -Wall -fopenmp -std=c11

all: spmv 

utils.o: utils.cu utils.h
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ -c $<

mmio.o: mmio.cu mmio.h
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ -c $<

spmv.o: spmv.cu
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ -c $<


spmv: spmv.o utils.o mmio.o
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) -o $@ $^
	rm -f ./*.o

clean:
	rm -f spmv *.o
