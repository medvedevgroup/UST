CC=g++
CFLAGS=-c -w -std=c++11 -O2

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main.out

BCALMFILE=/Volumes/exFAT/data2019/staphsub/31/list_reads.unitigs.fa
K=31

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o main.out stitchedUnitigs.fa 

test:
	./main.out -i  $(BCALMFILE) -k $(K) -f 1 -m 10 -a 1
