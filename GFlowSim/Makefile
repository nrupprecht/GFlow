CC = icpc
MPICC = mpic++
CFLAG = -std=c++14
FLAGS = $(CFLAG) -g -O3 -g #-ip
#MKLROOT = /afs/crc.nd.edu/x86_64_linux/intel/15.0/mkl
LDLIBS = #-lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

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
