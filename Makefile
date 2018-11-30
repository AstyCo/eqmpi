CC=mpixlcxx_r
EXTRAS_CFLAGS=-qarch=450d -qtune=450
CFLAGS=-O3 -qsmp=omp -qstrict
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t2
ARGUMENTS=

one:
	$(CC) $(CFLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)

all: $(SOURCES) $(EXECUTABLE)
	
clean: 
	rm -rf *.o
	
submit_polus:
	# TODO
	mpisubmit.pl -p 1 -w 00:30 -g ./$(EXECUTABLE) $(ARGUMENTS)

bsub_polus:
	# TODO
	bsub <bsub_args
	
submit_bluegene:
	# TODO
run:
	$(EXECUTABLE) $(ARGUMENTS)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(INC_PARAMS) $(CFLAGS) $< -o $@
	