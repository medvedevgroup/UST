CC=g++
CFLAGS=-c -w -std=c++11 -O2

SOURCES=ust.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=ust

BCALMFILE=example/list_reads.unitigs.fa
K=31

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o ust stitchedUnitigs.fa 

test:
	./ust -i  $(BCALMFILE) -k $(K) -a 1
