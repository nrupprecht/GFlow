CC = icpc
MPICC = mpic++
CFLAG = -std=c++14 -Wall
OPTFLAGS = -xHost -O3 -no-prec-div -static -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size
#QOPT = # -qopt-report 5
#OPTFLAGS = -O3 -static -funroll-loops -mllvm -vectorize-slp -mllvm -bb-vectorize-aligned-only -ffast-math -ffinite-math-only -freciprocal-math -ftree-vectorize -fslp-vectorize-aggressive -fstrict-vtable-pointers -fwhole-program-vtables -flto
#QOPT = # -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize
FLAGS = $(CFLAG) -g $(QOPT) $(OPTFLAGS)

# All object files needed
CONTROL_OBJECTS = build/Sectorization.o build/SimData.o build/ForceHandler.o
CREATION_OBJECTS = build/Creator.o build/FileParser.o build/DataCreator.o
DATA_OBJECTS    = build/DataRecord.o 
INT_OBJECTS     = build/VelocityVerletIntegrator.o build/Integrator.o
OBJ_OBJECTS     = build/Particle.o build/Wall.o
INC_OBJECTS     = build/vec2d.o
TEST_OBJECTS    = build/test_aligned_array.o
FORCE_OBJECTS   = build/ViscousDrag.o build/ConstantAcceleration.o

ALL_OBJECTS     = $(CONTROL_OBJECTS) $(CREATION_OBJECTS) $(DATA_OBJECTS) $(INT_OBJECTS) $(OBJ_OBJECTS) $(INC_OBJECTS) $(FORCE_OBJECTS)

all: bin/gflow # bin/test_aligned_array

###### Executables #####
bin/gflow: build/gflow.o $(ALL_OBJECTS)
	$(MPICC) $^ -o $@

bin/test_aligned_array: build/test_aligned_array.o
	$(MPICC) $^ -o $@

#####   Objects   #####
# In /src/
build/%.o: src/applications/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/control/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/creation/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/data/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/forces/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/integrators/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)

build/%.o: src/objects/%.cpp
	$(MPICC) -c -o $@ $< $(FLAGS)


# Not in /src/
build/%.o: include/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

build/%.o: test/%.cpp
	$(MPICC) $(FLAGS) -c -o $@ $<

#####   Clean    #####
.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
