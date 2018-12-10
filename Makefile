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
	
bgpcanc:
	llcancel -u edu-cmc-skmodel18-624-03
	
seqbgp:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)
seqbgpf:
	$(CC) $(CFLAGS) $(BLUEGENE_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DFLOAT_P -o $(EXECUTABLE_FLOAT)
	
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
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE)

bgpall:
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT)
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) 
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT)
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) 
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT)
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) 
bgpallomp:
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) 
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT)
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) 
	
bgpsep:
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_FLOAT) N=128
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_FLOAT) N=256
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_FLOAT) N=512
bgpsepomp:
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_OMP_FLOAT) N=128
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_OMP_FLOAT) N=256
	mpisubmit.bg -n 1 -m SMP -w 00:15:00 $(EXECUTABLE_OMP_FLOAT) N=512

	
bgp512:
	mpisubmit.bg -n 1 -m SMP -w 00:30:00 $(EXECUTABLE_FLOAT) N=512
bgp512omp:
	mpisubmit.bg -n 1 -m SMP -w 00:30:00 $(EXECUTABLE_OMP_FLOAT) N=512

	
bgp1024:
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1024
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1024
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1024
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1024
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1024
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1024
	
bgp1536:
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1536
	mpisubmit.bg -n 128 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1536
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1536
	mpisubmit.bg -n 256 -m SMP -w 00:10:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1536
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_FLOAT) N=1536
	mpisubmit.bg -n 512 -m SMP -w 00:5:00 -e  "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP_FLOAT) N=1536
	
	
run:
	$(EXECUTABLE) $(ARGUMENTS)
	