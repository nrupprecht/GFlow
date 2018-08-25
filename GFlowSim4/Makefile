# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

COMPILER ?= 0
USE_MPI  ?= 0
DEBUG    ?= 0
NO_SIMD  ?= 0

# Add flags to this variable as neccessary
FLAGS = 

# Compiler options
ifeq ($(COMPILER), 0)
CC = icpc
FLAGS += -D _INTEL_=1
else
CC = g++
FLAGS += -D _CLANG_=1
endif 
MPICC = mpic++

ifeq ($(DEBUG), 1)
FLAGS += -D DEBUG=1
endif

ifeq ($(NO_SIMD), 1)
FLAGS += -D SIMD_TYPE=0
endif

# Whether to use MPI or not
ifeq (USE_MPI, 0)
COMP = $(CC)
FLAGS += -D USE_MPI=1
else
COMP = $(MPICC)
endif

# Name of the application to compile
APP   = driver
# Directories for source, object, bin files
SRC   = src
OBJ   = obj
BIN   = bin

EXCLUDE = src/applications/percolationanalysis.cpp src/applications/pressureanalysis.cpp src/applications/driver.cpp

SRCS1   := $(shell find $(SRC) -name '*.cpp')
SRCS    := $(filter-out $(EXCLUDE), $(SRCS1))
SRCDIRS := $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.cpp,$(OBJ)/%.o,$(SRCS))

CSTD     = -std=c++11
DEBUG    = -g
ifeq ($(COMPILER), 0)
OPTIMIZE := -O3 -funroll-loops -no-prec-div -restrict -no-inline-max-size -ffast-math
else
OPTIMIZE := -O3
endif
#OPENCV   := $(shell pkg-config --cflags --libs opencv)
CFLAGS   := -c $(CSTD) $(DEBUG) $(OPTIMIZE) $(FLAGS)

all: $(BIN)/$(APP)

$(BIN)/$(APP): $(OBJS) obj/driver.o
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) -o $@ obj/driver.o $(OBJS) $(LDFLAGS) #$(OPENCV)

$(OBJ)/%.o: %.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

obj/driver.o: src/applications/driver.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -r -f $(OBJ) $(BIN)/*