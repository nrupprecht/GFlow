CC = icpc
OPT = -parallel
FLAGS = -std=c++14 -g -O3 $(OPT)
targets = driver control controlPhi Jamming JamShape time tune theory
files = Simulator.o Object.o Field.o

all: $(targets)

# Executables
control: control.o $(files)
	$(CC) $(OPT) $^ -o $@

controlPhi: controlPhi.o $(files)
	$(CC) $(OPT) $^ -o $@

tune: tune.o $(files)
	$(CC) $(OPT) $^ -o $@

Jamming: Jamming.o $(files)
	$(CC) $(OPT) $^ -o $@

JamShape: JamShape.o $(files)
	$(CC) $(OPT) $^ -o $@

driver: driver.o $(files)
	$(CC) $(OPT) $^ -o $@

time: time.o $(files)
	$(CC) $(OPT) $^ -o $@

theory: theory.o
	$(CC) $^ -o $@

# Object files
%.o : %.cpp
	$(CC) -c $(FLAGS) $^

theory.o : theory.cpp
	$(CC) -c $^

.PHONY: clean
clean:
	rm *.o $(targets)

.PHONY: fclean
fclean:
	rm MAC.o ftest.o ftest
