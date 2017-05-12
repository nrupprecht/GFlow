CC = icpc
MPICC = mpic++
CFLAG = -std=c++14 -Wall
OPTFLAGS = -xHost -O3 -no-prec-div -static -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size #-qopt-report=5
FLAGS = $(CFLAG) -g # $(OPTFLAGS)

OBJECTS = build/VelocityVerletIntegrator.o build/Integrator.o build/SimData.o build/DataRecord.o

bin/gflow: build/gflow.o $(OBJECTS)
	$(MPICC) $^ -o $@

# Objects
build/%.o: src/applications/%.cpp
	$(CC) -c -o $@ $<

build/%.o: src/integrators/%.cpp
	$(CC) -c -o $@ $<

build/%.o: src/control/%.cpp
	$(CC) -c -o $@ $<

build/%.o: src/data/%.cpp
	$(CC) -c -o $@ $<

.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
