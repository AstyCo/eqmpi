CCPOL=mpixlC
POLUS_FLAGS=-DPOLUS -qarch=pwr8
CFLAGS=-O3 -qstrict
LDFLAGS=-lm
SOURCES=iterations.cpp globals.cpp utils.cpp main.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=t3_seq

seq:
	$(CCPOL) $(CFLAGS) $(POLUS_FLAGS) $(INC_PARAMS) $(LDFLAGS) $(SOURCES) -o $(EXECUTABLE)

cl:
	rm -f *.o *.out *.err core.*
	
subm:
	bsub <bsub_args_1
