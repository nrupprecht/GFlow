# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

USE_MPI ?= 0

# Compiler options
CC = icpc
MPICC = mpic++

# Add flags to this variable as neccessary
FLAGS = 

# Whether to use MPI or not
#ifneq (USE_MPI, 1)
#	COMP = $(MPICC)
#	FLAGS += -D USE_MPI=1
#else
#	COMP = $(CC)
#endif

APP   = driver

SRC   = src
OBJ   = obj
BIN   = bin

EXCLUDE = src/applications/percolationanalysis.cpp src/applications/pressureanalysis.cpp src/applications/driver.cpp

SRCS1   := $(shell find $(SRC) -name '*.cpp')
SRCS    := $(filter-out $(EXCLUDE), $(SRCS1))
SRCDIRS := $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.cpp,$(OBJ)/%.o,$(SRCS))

#OPENCV  := $(shell pkg-config --cflags --libs opencv)

CSTD     = -std=c++11
DEBUG    = -g -Wall
OPTIMIZE = -O3 -funroll-loops -no-prec-div -restrict -no-inline-max-size -ffast-math
CFLAGS   = -c $(CSTD) $(DEBUG) $(OPTIMIZE) $(FLAGS)

all: $(BIN)/$(APP)

$(BIN)/$(APP): $(OBJS) obj/driver.o
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) obj/driver.o $(OBJS) $(LDFLAGS) -o $@

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