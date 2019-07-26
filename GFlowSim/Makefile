# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

COMPILER ?= 1
USE_MPI  ?= 0
DEBUG    ?= 0
NO_SIMD  ?= 0
DO_OPTIMIZE ?= 1

# Add flags to this variable as neccessary
FLAGS = 

MPICC = mpic++

# Compiler options
ifeq ($(USE_MPI), 1)
CC = $(MPICC)
FLAGS += -D USE_MPI=1
COMPILER = 2
else ifeq ($(COMPILER), 0)
CC = icpc
FLAGS += -D _INTEL_=1
else ifeq ($(COMPILER), 1)
CC = g++
FLAGS += -D _CLANG_=1
else
CC = clang++
FLAGS += -D _CLANG_=1
endif 

ifeq ($(DEBUG), 1)
FLAGS += -D DEBUG=1
endif

ifeq ($(NO_SIMD), 1)
FLAGS += -D SIMD_TYPE=0
endif

# Name of the application to compile
APP   = driver
# Directories for source, object, bin files
SRC   = src
OBJ   = obj
BIN   = bin
APPLICATIONS = $(SRC)/applications

EXCLUDE := $(shell find $(APPLICATIONS) -name '*.cpp')
EXCLUDE += src/utility/argparse2.cpp

SRCS1   := $(shell find $(SRC) -name '*.cpp')
SRCS    := $(filter-out $(EXCLUDE), $(SRCS1))
SRCDIRS := $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.cpp,$(OBJ)/%.o,$(SRCS))

CSTD     = -std=c++11
DEBUG    = -g
ifeq ($(COMPILER), 0) 
ifeq ($(DO_OPTIMIZE), 1)
OPTIMIZE := -O3 -funroll-loops -no-prec-div -restrict -no-inline-max-size -ffast-math
endif
else
ifeq ($(DO_OPTIMIZE), 1)
OPTIMIZE := -O3
endif
endif
#OPENCV   := $(shell pkg-config --cflags --libs opencv)
CFLAGS   := -c $(CSTD) $(DEBUG) $(OPTIMIZE) $(FLAGS)

all: $(BIN)/$(APP) $(BIN)/vistools
#$(BIN)/lineattraction $(BIN)/wallattraction $(BIN)/polymerattraction $(BIN)/demonrates

# Executables

$(BIN)/$(APP): obj/driver.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/driver.o $(OBJS)

$(BIN)/vistools: obj/vistools.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/vistools.o $(OBJS)

#$(BIN)/lineattraction: obj/lineattraction.o $(OBJS)
#	@mkdir -p `dirname $@`
#	@echo "Linking $@..."
#	@$(CC) -o $@ obj/lineattraction.o $(OBJS)

#$(BIN)/wallattraction: obj/wallattraction.o $(OBJS)
#	@mkdir -p `dirname $@`
#	@echo "Linking $@..."
#	@$(CC) -o $@ obj/wallattraction.o $(OBJS)

#$(BIN)/polymerattraction: obj/polymerattraction.o $(OBJS)
#	@mkdir -p `dirname $@`
#	@echo "Linking $@..."
#	@$(CC) -o $@ obj/polymerattraction.o $(OBJS)

#$(BIN)/demonrates: obj/demonrates.o $(OBJS)
#	@mkdir -p `dirname $@`
#	@echo "Linking $@..."
#	@$(CC) -o $@ obj/demonrates.o $(OBJS)

# General object files

$(OBJ)/%.o: %.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

# Excluded object files

$(OBJ)/driver.o: src/applications/driver.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

$(OBJ)/vistools.o: src/applications/vistools.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

#$(OBJ)/lineattraction.o: src/applications/lineattraction.cpp
#	@mkdir -p `dirname $@`
#	@echo "Compiling $<..."
#	@$(CC) $(CFLAGS) $< -o $@

#$(OBJ)/wallattraction.o: src/applications/wallattraction.cpp
#	@mkdir -p `dirname $@`
#	@echo "Compiling $<..."
#	@$(CC) $(CFLAGS) $< -o $@

#$(OBJ)/polymerattraction.o: src/applications/polymerattraction.cpp
#	@mkdir -p `dirname $@`
#	@echo "Compiling $<..."
#	@$(CC) $(CFLAGS) $< -o $@

#$(OBJ)/demonrates.o: src/applications/demonrates.cpp
#	@mkdir -p `dirname $@`
#	@echo "Compiling $<..."
#	@$(CC) $(CFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -r -f $(OBJ) $(BIN)/*

.PHONY : forces
forces:
	rm obj/src/interactions/*
	make

.PHONY : creators
creators:
	rm obj/src/creators/*
	make
