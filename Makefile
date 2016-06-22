CC = icpc
FLAGS = -std=c++11 -O3 -g

targets = driver time fdriver Jamming test
files = Simulator.o Object.o Field.o

all: $(targets)

# Executables
test: test.o Field.o
	$(CC) $^ -o $@

Jamming: Jamming.o $(files)
	$(CC) $^ -o $@

driver: driver.o $(files)
	$(CC) $^ -o $@

time: time.o $(files)

fdriver: fdriver.o Fluid.o
	$(CC) $^ -o $@

# Object files
%.o : %.cpp
	$(CC) -c $(FLAGS) $<

.PHONY: clean
clean:
	rm *.o $(targets)

.PHONY: fclean
fclean:
	rm Fluid.o fdriver.o fdriver