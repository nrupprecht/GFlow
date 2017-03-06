CC = icpc
FLAGS = -std=c++14 -g -O3 
OPT = #-openmp -parallel
targets = driver
oldTargets = Jamming JamShape time
MKLROOT = /afs/crc.nd.edu/x86_64_linux/intel/15.0/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

files = Simulator.o Object.o Field.o Other/Tensor.o Sectorization.o

all: $(targets)

# Executables
driver: driver.o $(files)
	$(CC) $(OPT) $^ -o $@ $(LDLIBS)

# Currently unused Executables
Jamming: Jamming.o $(files)
	$(CC) $(OPT) $^ -o $@ $(LDLIBS)

JamShape: JamShape.o $(files)
	$(CC) $(OPT) $^ -o $@ $(LDLIBS)

time: time.o $(files)
	$(CC) $(OPT) $^ -o $@ $(LDLIBS)

# Object files
%.o : %.cpp
	$(CC) $(OPT) -c $(FLAGS) $^

.PHONY: clean
clean:
	rm *.o $(targets)
