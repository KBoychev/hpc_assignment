C=g++
CFLAGS=-c -I./
FNCSD=./functions/

default: task1 clean run

all: task1 clean run

task1: main.o set_M_e.o set_K_e.o set_F_e.o get_M.o get_K.o get_F.o inv.o disp.o	
	$(C) main.o $(FNCSD)set_M_e.o $(FNCSD)set_K_e.o $(FNCSD)set_F_e.o $(FNCSD)get_M.o $(FNCSD)get_K.o $(FNCSD)get_F.o $(FNCSD)inv.o $(FNCSD)disp.o -llapack -lblas -o task1.out

main.o: 	
	$(C) $(CFLAGS) main.cpp -o main.o

set_M_e.o:
	$(C) $(CFLAGS) $(FNCSD)set_M_e.cpp -o $(FNCSD)set_M_e.o

set_K_e.o:
	$(C) $(CFLAGS) $(FNCSD)set_K_e.cpp -o $(FNCSD)set_K_e.o

set_F_e.o:

	$(C) $(CFLAGS) $(FNCSD)set_F_e.cpp -o $(FNCSD)set_F_e.o

get_M.o:
	$(C) $(CFLAGS) $(FNCSD)get_M.cpp -o $(FNCSD)get_M.o

get_K.o:
	$(C) $(CFLAGS) $(FNCSD)get_K.cpp -o $(FNCSD)get_K.o

get_F.o:
	$(C) $(CFLAGS) $(FNCSD)get_F.cpp -o $(FNCSD)get_F.o

inv.o:
	$(C) $(CFLAGS) $(FNCSD)inv.cpp -o $(FNCSD)inv.o

disp.o:
	$(C) $(CFLAGS) $(FNCSD)disp.cpp -o $(FNCSD)disp.o

clean:
	rm main.o
	rm $(FNCSD)set_M_e.o
	rm $(FNCSD)set_K_e.o
	rm $(FNCSD)set_F_e.o
	rm $(FNCSD)get_M.o
	rm $(FNCSD)get_K.o
	rm $(FNCSD)get_F.o

run:
	./task1.out
