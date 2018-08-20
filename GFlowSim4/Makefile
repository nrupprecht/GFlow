# Makefile snatched from <https://ubuntuforums.org/showthread.php?t=1204739>

CC = icpc
MPICC = mpic++
COMP = $(CC)
CFLAG = -std=c++11 -Wall

APP   = driver

SRC   = src
OBJ   = obj
BIN   = bin

SRCS    := $(shell find $(SRC) -name '*.cpp')
SRCDIRS := $(shell find . -name '*.cpp' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.cpp,$(OBJ)/%.o,$(SRCS))

CSTD     = -std=c++11
DEBUG    = -g -Wall
OPTIMIZE = -O3 -funroll-loops -no-prec-div -restrict -no-inline-max-size -ffast-math
CFLAGS   = -c $(CSTD) $(DEBUG) $(OPTIMIZE)

all: $(BIN)/$(APP)

$(BIN)/$(APP): obj/driver.o $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) obj/driver.o $(OBJS) $(LDFLAGS) -o $@

obj/driver.o: applications/driver.cpp
	@mkdir -p `dirname $@`
	@$(CC) $(CFLAGS) $< -o $@

$(OBJ)/%.o: %.cpp
	@mkdir -p `dirname $@`
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) $< -o $@


#@echo "Generating dependencies for $<..."
#@$(call make-depend,$<,$@,$(subst .o,.d,$@))


#####   Clean    #####
.PHONY : clean
clean:
	$(RM) -r $(OBJ)