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
EXECUTABLE_OMP=t2_omp
ARGUMENTS=


all: $(SOURCES) $(EXECUTABLE)
	
seqbgp:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)
ompbgp:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(OMP_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -o $(EXECUTABLE_OMP)
	
seqpol:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)
omppol:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(OMP_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -o $(EXECUTABLE_OMP)
	
clean: 
	rm -rf *.o
cl:
	rm -f *.o *.out *.err core.*
	
# submit polus
pol:
	bsub <bsub_args_1
polall:
	bsub <bsub_args_32
	bsub <bsub_args_16
	bsub <bsub_args_8
	
# submit bluegene
bgp: 
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE)

bgpall: 
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE)
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP) 
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE)
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP) 
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE)
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP) 
	
	
run:
	$(EXECUTABLE) $(ARGUMENTS)
	