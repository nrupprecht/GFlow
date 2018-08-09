CC = icpc
MPICC = mpic++
COMP = $(CC)
CFLAG = -std=c++11 -Wall
# WARNING: DO NOT USE -xHost IT MAKES THINGS SLOWER! -fast uses -xHost, so don't use that either
OPTFLAGS = -O3 -funroll-loops $(OPTFLAGS2)
OPTFLAGS2 = -no-prec-div -restrict -no-inline-max-size #-fp-model fast=2
FLAGS = $(CFLAG) -g $(OPTFLAGS) # -qopt-report 5

# Option: -vec-report : Reports if loops were vectorized

# Get all object files
SRC_DIR := src
OBJ_DIR := build
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all: bin/driver bin/pressureanalysis bin/percolationanalysis

###### Executables #####
bin/driver: build/driver.o $(OBJ_FILES)
	$(COMP) $^ -o $@

bin/pressureanalysis: build/pressureanalysis.o $(OBJ_FILES)
	$(COMP) $^ -o $@

bin/percolationanalysis: build/percolationanalysis.o $(OBJ_FILES)
	$(COMP) $^ -o $@

#####   Objects   #####
# In /src/
build/%.o: src/%.cpp
	$(COMP) -c -o $@ $< $(FLAGS)

build/%.o: applications/%.cpp
	$(COMP) -c -o $@ $< $(FLAGS)

#####   Clean    #####
.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
