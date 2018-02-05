CC = icpc
MPICC = mpic++
COMP = $(CC)
CFLAG = -std=c++14 -Wall
OPTFLAGS = -xHost -O3 -no-prec-div -static -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size

#QOPT = -qopt-report 5
#OPTFLAGS = -O3 -static -funroll-loops -mllvm -vectorize-slp -mllvm -bb-vectorize-aligned-only -ffast-math -ffinite-math-only -freciprocal-math -ftree-vectorize -fslp-vectorize-aggressive -fstrict-vtable-pointers -fwhole-program-vtables -flto
#QOPT = -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize

FLAGS = $(CFLAG) -g $(OPTFLAGS) $(QOPT)

# All object files needed
CONTROL_OBJECTS = build/Sectorization.o build/SimData.o build/StatusObject.o build/TerminationCondition.o build/ForceHandler.o build/Characteristic.o
SIMULATION_OBJECTS = build/FractureSimulation.o build/StandardSimulation.o build/SimulationBase.o
CREATION_OBJECTS = build/Creator.o build/FileParser.o build/DataCreator.o
DATA_OBJECTS    = build/DataRecord.o 
INT_OBJECTS     = build/VelocityVerletIntegrator.o build/Integrator.o build/StandardSectorization.o
OBJ_OBJECTS     = build/Particle.o build/Wall.o build/ScalarField.o
INC_OBJECTS     = build/vec2d.o
FORCE_OBJECTS   = build/ViscousDrag.o build/ConstantAcceleration.o build/TemperatureForce.o

ALL_OBJECTS     = $(CONTROL_OBJECTS) $(SIMULATION_OBJECTS) $(CREATION_OBJECTS) $(DATA_OBJECTS) $(INT_OBJECTS) $(OBJ_OBJECTS) $(INC_OBJECTS) $(FORCE_OBJECTS) $(TEST_OBJECTS)

all: bin/gflow

###### Executables #####
bin/gflow: build/gflow.o $(ALL_OBJECTS) | bin
	$(COMP) $^ -o $@

# Create directories if they don't exist
bin:
	mkdir -p $@

build:
	mkdir -p $@

#####   Objects   #####
# In /src/
build/%.o: src/applications/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/control/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/simulation/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/creation/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/data/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/forces/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/integrators/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: src/objects/%.cpp | build
	$(COMP) -c -o $@ $< $(FLAGS)

# Not in /src/
build/%.o: include/%.cpp | build
	$(COMP) $(FLAGS) -c -o $@ $<

build/%.o: test/%.cpp | build
	$(COMP) $(FLAGS) -c -o $@ $<

#####   Clean    #####
.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
