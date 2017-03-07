C=g++ -std=c++11
CFLAGS=-c -I./
FNCSD=./functions/

default: task clean run1 run2 run3

all: task clean run1 run2 run3

task: main.o set_M_e.o set_K_e.o set_F_e.o get_M.o get_K.o get_F.o inv.o disp.o	
	$(C) main.o $(FNCSD)set_M_e.o $(FNCSD)set_K_e.o $(FNCSD)set_F_e.o $(FNCSD)get_M.o $(FNCSD)get_K.o $(FNCSD)get_F.o $(FNCSD)inv.o $(FNCSD)disp.o -llapack -lblas -o task.out

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
	rm $(FNCSD)inv.o
	rm $(FNCSD)disp.o

run1:
	./task.out -L 10.0 -N_e 4 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -N_t 10000 -rho 7850.0 -eq 0

run2:
	./task.out -L 10.0 -N_e 4 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -N_t 10000 -rho 7850.0 -eq 1 -sch 0
	
run3:
	./task.out -L 10.0 -N_e 4 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -N_t 10000 -rho 7850.0 -eq 1 -sch 1
