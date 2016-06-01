CC = icpc
FLAGS = -std=c++11 -O3 -g

targets = driver

all: $(targets)

# Executables

driver: driver.o Simulator.o Object.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $(FLAGS) $<

.PHONY: clean
clean:
	rm *.o $(targets)