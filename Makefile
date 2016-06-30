CC = icpc
FLAGS = -std=c++11 -g #-O3

targets = driver time Jamming test ftest
files = Simulator.o Object.o Field.o

all: $(targets)

# Executables
ftest: ftest.o Field.o MAC.o
	$(CC) $^ -o $@

test: test.o Field.o
	$(CC) $^ -o $@

Jamming: Jamming.o $(files)
	$(CC) $^ -o $@

driver: driver.o $(files)
	$(CC) $^ -o $@

time: time.o $(files)
	$(CC) $^ -o $@

# Object files
%.o : %.cpp
	$(CC) -c $(FLAGS) $^

.PHONY: clean
clean:
	rm *.o $(targets)

.PHONY: fclean
fclean:
	rm Fluid.o MAC.o ftest.o ftest