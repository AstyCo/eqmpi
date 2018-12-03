CC=mpixlcxx_r
EXTRAS_CFLAGS=-qarch=450d -qtune=450 -qsmp=omp
CFLAGS=-O3 -qstrict -qarch=450d -qtune=450 -qsmp=omp
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t2
EXECUTABLE_OMP=t2_omp
ARGUMENTS=


all: $(SOURCES) $(EXECUTABLE)
	
one:
	$(CC) $(CFLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)
omp:
	$(CC) $(CFLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -DWITH_OMP -o $(EXECUTABLE_OMP)
clean: 
	rm -rf *.o
cl:
	rm -f *.o *.out *.err core.*
	
# submit polus
pol:
	# TODO
	mpisubmit.pl -p 1 -w 00:05 -g ./$(EXECUTABLE) $(ARGUMENTS)

bsub_polus:
	# TODO
	bsub <bsub_args
	
# submit bluegene
bgp: 
	mpisubmit.bg -n 32 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE)
ompb: 
	mpisubmit.bg -n 32 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE_OMP) 
	
	
run:
	$(EXECUTABLE) $(ARGUMENTS)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(INC_PARAMS) $(CFLAGS) $< -o $@
	