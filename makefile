C=g++
CFLAGS=-c

default: task1 run

all: task1 run

task1: main.o
	$(C) main.o -llapack -lblas -o task1.out

main.o: 
	$(C) $(CFLAGS) main.cpp -o main.o

clean:
	rm task1.out
	rm main.o
run:
	./task1.out
