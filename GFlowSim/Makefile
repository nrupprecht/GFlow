CC = icpc
MPICC = mpic++
CFLAG = -std=c++14 -Wall
FLAGS = $(CFLAG) -g -O3 #-g

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
	rm *.o $(targets)
