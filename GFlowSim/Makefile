CC = icpc
MPICC = mpic++
CFLAG = -std=c++14 -Wall
OPTFLAGS = -xHost -O3 -no-prec-div -static -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size #-qopt-report=5
FLAGS = $(CFLAG) -g # $(OPTFLAGS)

# All object files needed
CONTROL_OBJECTS = build/Creator.o build/Sectorization.o build/SimData.o build/ForceHandler.o
DATA_OBJECTS    = build/DataRecord.o
INT_OBJECTS     = build/VelocityVerletIntegrator.o build/Integrator.o
OBJ_OBJECTS     = build/Particle.o build/Wall.o
INC_OBJECTS     = build/vec2d.o
TEST_OBJECTS    = build/test_aligned_array.o

all: bin/gflow bin/test_aligned_array

bin/gflow: build/gflow.o $(CONTROL_OBJECTS) $(DATA_OBJECTS) $(INT_OBJECTS) $(OBJ_OBJECTS) $(INC_OBJECTS)
	$(MPICC) $^ -o $@

bin/test_aligned_array: build/test_aligned_array.o
	$(MPICC) $^ -o $@

# Objects
build/%.o: src/applications/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: src/integrators/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: src/control/%.cpp
	$(MPICC) $(FLAGS)-c -o $@ $<

build/%.o: src/data/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: src/objects/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: include/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: test/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
