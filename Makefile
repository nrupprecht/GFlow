CC = icpc
FLAGS = -std=c++14 -g -O3

targets = driver control Jamming JamShape time ftest tune
files = Simulator.o Object.o Field.o

all: $(targets)

# Executables
control: control.o $(files)
	$(CC) $^ -o $@

tune: tune.o $(files)
	$(CC) $^ -o $@

ftest: ftest.o GFlow.o MAC.o Object.o
	$(CC) $^ -o $@

Jamming: Jamming.o $(files)
	$(CC) $^ -o $@

JamShape: JamShape.o $(files)
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
	rm MAC.o ftest.o ftest
