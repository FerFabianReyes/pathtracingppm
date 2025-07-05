CPP=g++
CPPFLAGS=-O3 -fopenmp 

default: rt

all: rt

rt: rt.cpp Makefile
	$(CPP) $(CPPFLAGS) -o rt rt.cpp

clean:
	-rm rt
