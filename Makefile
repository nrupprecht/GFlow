CC = icpc
FLAGS = -std=c++11 -O3 -g

targets = driver time
files = Simulator.o Object.o

all: $(targets)

# Executables

driver: driver.o $(files)
	$(CC) $^ -o $@

time: time.o $(files)

# Object files
%.o : %.cpp
	$(CC) -c $(FLAGS) $<

.PHONY: clean
clean:
	rm *.o $(targets)