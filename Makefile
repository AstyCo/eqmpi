CC=mpixlcxx_r
CCPOL=mpixlC
OMP_FLAGS=-qsmp=omp
BLUEGENE_FLAGS=-DBGP -qarch=450d -qtune=450
POLUS_FLAGS=-DPOLUS -qarch=pwr8
CFLAGS=-O3 -qstrict
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t2
EXECUTABLE_FLOAT=t2_float
EXECUTABLE_OMP=t2_omp
EXECUTABLE_OMP_FLOAT=t2_omp_float
EXECUTABLE_FLOAT_SEQ=t2_seq
ARGUMENTS=


all: $(SOURCES) $(EXECUTABLE)
	
ompbgp:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(OMP_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -o $(EXECUTABLE_OMP)
ompbgpf:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(OMP_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -DFLOAT_P -o $(EXECUTABLE_OMP_FLOAT)
	
seqpol:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DSEQ -o $(EXECUTABLE_FLOAT_SEQ)
seqpolf:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE_FLOAT)
omppol:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(OMP_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -o $(EXECUTABLE_OMP)
	
cl:
	rm -f *.o *.out *.err core.*
	
# submit polus
subm:
	bsub <bsub_args_cuda
	
cuda: 
	nvcc -rdc=true -arch=sm_60 -ccbin mpixlC --x cu -Xcompiler -O0,-qarch=pwr8,-qstrict,-Wall cuda.cu $(SOURCES) -o t3
	
run:
	$(EXECUTABLE) $(ARGUMENTS)
	