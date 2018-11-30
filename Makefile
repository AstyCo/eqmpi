CC=mpixlcxx_r
EXTRAS_CFLAGS=-qarch=450d -qtune=450
CFLAGS=-O3 -qsmp=omp -qstrict $(EXTRAS_CFLAGS)
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t2
ARGUMENTS=


all: $(SOURCES) $(EXECUTABLE)
	
one:
	$(CC) $(CFLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)
clean: 
	rm -rf *.o
	
submit_polus:
	# TODO
	mpisubmit.pl -p 1 -w 00:30 -g ./$(EXECUTABLE) $(ARGUMENTS)

bsub_polus:
	# TODO
	bsub <bsub_args
	
submit_bluegene:
	mpisubmit.bg -n 32 -m SMP -w 00:15:00 -e "OMP_NUM_THREADS=4" $(EXECUTABLE)
	
run:
	$(EXECUTABLE) $(ARGUMENTS)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(INC_PARAMS) $(CFLAGS) $< -o $@
	