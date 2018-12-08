# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

COMPILER ?= 0
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
else
CC = g++
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

EXCLUDE = src/applications/percolationanalysis.cpp src/applications/pressureanalysis.cpp src/applications/driver.cpp src/applications/vistools.cpp
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

# Executables

$(BIN)/$(APP): obj/driver.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/driver.o $(OBJS) $(LDFLAGS)

$(BIN)/vistools: obj/vistools.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/vistools.o $(OBJS) $(LDFLAGS)

# General object files

$(OBJ)/%.o: %.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

# Excluded object files

obj/driver.o: src/applications/driver.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

obj/vistools.o: src/applications/vistools.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

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
