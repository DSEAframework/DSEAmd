MAKEFLAGS += -j8

CPP_COMPIPLER=mpic++
CUDA_COMPIPLER=nvcc

FLAGS_ALL_CPP = -g -I./ -I $(CUDA_HOME)/include/ -I /fscratch/sw/source/ucx-1.17.0/src/ -I /fscratch/sw/ucx-1.17.0_dual_rail/include/ -I/zhome/academic/HLRS/hlrs/hpcmaros/sw/source/ucx-1.17.0/src
FLAGS_CUDA = -arch=native -g  -I./ --generate-line-info -I/sw/hawk-ai-rh8/hlrs/non-spack/rev-009_2022-09-01/mpi/openmpi/openmpi-4.1.5-cuda-11.8/include -I/zhome/academic/HLRS/hlrs/hpcmaros/sw/source/ucx-1.17.0/src

all: dsea dsea_ucx_send dsea_ucx_rec

dsea: obj/main.o obj/dsea_gpu.o obj/dsea_kernel.o obj/dsea_input.o obj/dsea_output.o obj/dsea_storage.o obj/dsea_main.o obj/dsea.o obj/dsea_block.o
	$(CPP_COMPIPLER) -g -fopenmp \
	obj/main.o \
	obj/dsea.o \
	obj/dsea_block.o \
	obj/dsea_gpu.o \
	obj/dsea_kernel.o \
	obj/dsea_input.o \
	obj/dsea_output.o \
	obj/dsea_storage.o \
	obj/dsea_main.o \
	-o dsea \
	-L$(CUDA_HOME)/lib64 -lcudart

dsea_ucx_send: obj/main_ucx_send.o obj/dsea_gpu.o obj/dsea_kernel.o obj/dsea_input.o obj/dsea_output_ucx.o obj/dsea_storage.o obj/dsea_main.o obj/dsea.o obj/dsea_block.o
	$(CPP_COMPIPLER) -g -fopenmp \
	obj/main_ucx_send.o \
	obj/dsea.o \
	obj/dsea_block.o \
	obj/dsea_gpu.o \
	obj/dsea_kernel.o \
	obj/dsea_input.o \
	obj/dsea_output_ucx.o \
	obj/dsea_storage.o \
	obj/dsea_main.o \
	-o dsea_ucx_send \
	-L$(CUDA_HOME)/lib64 -lcudart \
	-L/fscratch/sw/ucx-1.17.0_dual_rail/lib -luct -lucs -lucp

dsea_ucx_rec: obj/main_ucx_rec.o obj/dsea_gpu.o obj/dsea_kernel.o obj/dsea_input_ucx.o obj/dsea_output.o obj/dsea_storage.o obj/dsea_main.o obj/dsea.o obj/dsea_block.o
	$(CPP_COMPIPLER) -g -fopenmp \
	obj/main_ucx_rec.o \
	obj/dsea.o \
	obj/dsea_block.o \
	obj/dsea_gpu.o \
	obj/dsea_kernel.o \
	obj/dsea_input_ucx.o \
	obj/dsea_output.o \
	obj/dsea_storage.o \
	obj/dsea_main.o \
	-o dsea_ucx_rec \
	-L$(CUDA_HOME)/lib64 -lcudart \
	-L/fscratch/sw/ucx-1.17.0_dual_rail/lib -luct -lucs -lucp

obj/main.o: main.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/main.o -c main.cpp

obj/main_ucx_send.o: main.cpp dsea.h makefile
	$(CPP_COMPIPLER) -D MRUCX_SEND $(FLAGS_ALL_CPP) -fopenmp -o obj/main_ucx_send.o -c main.cpp

obj/main_ucx_rec.o: main.cpp dsea.h makefile
	$(CPP_COMPIPLER) -D MRUCX_REC $(FLAGS_ALL_CPP) -fopenmp -o obj/main_ucx_rec.o -c main.cpp

obj/dsea.o: dsea.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea.o -c dsea.cpp

obj/dsea_block.o: dsea_block.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea_block.o -c dsea_block.cpp

obj/dsea_input.o: dsea_input.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea_input.o -c dsea_input.cpp

obj/dsea_input_ucx.o: dsea_input_ucx.cpp dsea.h makefile
	$(CUDA_COMPIPLER) $(FLAGS_CUDA) -Xcompiler -fopenmp -o obj/dsea_input_ucx.o -c dsea_input_ucx.cpp

obj/dsea_output.o: dsea_output.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea_output.o -c dsea_output.cpp

obj/dsea_output_ucx.o: dsea_output_ucx.cpp dsea.h makefile
	$(CUDA_COMPIPLER) $(FLAGS_CUDA) -Xcompiler -fopenmp -o obj/dsea_output_ucx.o -c dsea_output_ucx.cpp

obj/dsea_storage.o: dsea_storage.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea_storage.o -c dsea_storage.cpp

obj/dsea_main.o: dsea_main.cpp dsea.h makefile
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -fopenmp -o obj/dsea_main.o -c dsea_main.cpp

obj/dsea_gpu.o: dsea_gpu.cpp dsea.h makefile
	$(CUDA_COMPIPLER) $(FLAGS_CUDA) -I./ -o obj/dsea_gpu.o -c dsea_gpu.cpp

obj/dsea_kernel.o: dsea_kernel.cu dsea.h makefile
	$(CUDA_COMPIPLER) $(FLAGS_CUDA) -I./ -o obj/dsea_kernel.o -c dsea_kernel.cu



obj/dsea_store.o: dsea_store.cpp
	$(CPP_COMPIPLER) $(FLAGS_ALL_CPP) -o obj/dsea_store.o -c dsea_store.cpp

clean:
	rm obj/*.o dsea dsea-store dsea_ucx_send dsea_ucx_rec
