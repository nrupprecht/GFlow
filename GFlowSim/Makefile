CC = icpc
MPICC = mpic++
OPT = -O3 -g -ip
CFLAGS = -std=c++14 $(OPT)
MKLROOT = /afs/crc.nd.edu/x86_64_linux/intel/15.0/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

targets = driver
all:	$(targets)

# Executables
driver : driver.o GFlowBase.o Sectorization.o Particle.o
	$(MPICC) -o $@ $^ $(LDLIBS)

# Object files
%.o : %.cpp
	$(MPICC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) *.o *.stdout
