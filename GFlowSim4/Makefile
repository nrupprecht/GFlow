CC = icpc
MPICC = mpic++
CFLAG = -std=c++11 -Wall
OPTFLAGS = # -xHost -O3 -no-prec-div -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size
FLAGS = $(CFLAG) -g $(OPTFLAGS)

# All object files needed
ALL_OBJECTS = build/gflow.o build/simdata.o build/sectorization.o build/neighbors.o build/communicator.o build/base.o build/integrator.o build/velocityverlet.o build/force.o build/hard_sphere.o

all: bin/driver

###### Executables #####
bin/driver: build/driver.o $(ALL_OBJECTS)
	$(CC) $^ -o $@

#####   Objects   #####
# In /src/
build/%.o: src/%.cpp
	$(CC) -c -o $@ $< $(FLAGS)

#####   Clean    #####
.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
