CC = icpc
MPICC = mpic++
CFLAG = -std=c++11 -Wall
OPTFLAGS = # -xHost -O3 -no-prec-div -fp-model fast=2 -restrict -funroll-loops -no-inline-max-size
FLAGS = $(CFLAG) -g $(OPTFLAGS)

# Get all object files
SRC_DIR := src
OBJ_DIR := build
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all: bin/driver

###### Executables #####
bin/driver: build/driver.o $(OBJ_FILES)
	$(CC) $^ -o $@

#####   Objects   #####
# In /src/
build/%.o: src/%.cpp
	$(CC) -c -o $@ $< $(FLAGS)

#####   Clean    #####
.PHONY : clean
clean:
	rm -f build/*.o build/*.optrpt *.o *.optrpt bin/*
