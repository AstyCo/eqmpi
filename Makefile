CC=mpixlcxx
CFLAGS=-O3 -qsmp=omp
INC=
INC_PARAMS=$(foreach d, $(INC), -I$d)
LDFLAGS=
SOURCES=main.cpp utils.cpp globals.cpp iterations.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t2
ARGUMENTS=

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
	