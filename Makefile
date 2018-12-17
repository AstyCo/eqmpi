CCPOL=mpixlC
POLUS_FLAGS=-DPOLUS -qarch=pwr8
CFLAGS=-O3 -qstrict
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
EXECUTABLE_CUDA=t3
EXECUTABLE_FLOAT_SEQ=t2_seq

cuda: 
	nvcc -rdc=true -arch=sm_60 -ccbin mpixlC -Xcompiler -O0,-qarch=pwr8,-qstrict,-Wall cuda.cu $(SOURCES) -o $(EXECUTABLE_CUDA)
	
gpu_info:
	nvcc -rdc=true -arch=sm_60 -ccbin mpixlC -Xcompiler -O0,-qarch=pwr8,-qstrict,-Wall,-DSC_INFO cuda.cu $(SOURCES) -o $(EXECUTABLE_CUDA)
	
seqpol:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DSEQ -o $(EXECUTABLE_FLOAT_SEQ)
	
cl:
	rm -f *.o *.out *.err core.*
	
# submit polus
subm:
	bsub <bsub_args_cuda
