CC = icpc
MPICC = mpic++
CFLAG = -std=c++14 -Wall
FLAGS = $(CFLAG) -g -xHost -O3 -no-prec-div -static -fp-model fast=2 -restrict -funroll-loops
#-qopt-report 5

targets = driver

files = GFlowBase.o Sectorization.o Particle.o

all: $(targets)

# Executables
driver : driver.o $(files)
	$(MPICC) $^ -o $@ $(LDLIBS)

# Object files
%.o : %.cpp
	$(MPICC) -c $(FLAGS) $^

.PHONY : clean
clean:
	rm -f *.o *.optrpt $(targets)
