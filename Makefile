CC = icpc
OPT = 
FLAGS = -std=c++14 -g -O3 $(OPT)
targets = driver bacteria control controlPhi Jamming JamShape time tune solver
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

bacteria: bacteria.o $(files)
	$(CC) $(OPT) $^ -o $@

time: time.o $(files)
	$(CC) $(OPT) $^ -o $@

solver: solver.o Theory.o
	$(CC) $^ -o $@

# Object files
%.o : %.cpp
	$(CC) -c $(FLAGS) $^

.PHONY: clean
clean:
	rm *.o $(targets)
