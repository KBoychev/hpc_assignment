C=mpicxx -std=c++11 -O3
CFLAGS=-c -I./
FNCSD=./functions/

default: app clean

all: app clean

app: main.o get_M.o get_K.o get_K_eff.o get_F.o disp.o log.o
	$(C) main.o $(FNCSD)get_M.o $(FNCSD)get_K.o $(FNCSD)get_K_eff.o $(FNCSD)get_F.o $(FNCSD)disp.o $(FNCSD)log.o -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -o app.out

main.o: 	
	$(C) $(CFLAGS) main.cpp -o main.o

get_M.o:
	$(C) $(CFLAGS) $(FNCSD)get_M.cpp -o $(FNCSD)get_M.o

get_K.o:
	$(C) $(CFLAGS) $(FNCSD)get_K.cpp -o $(FNCSD)get_K.o

get_K_eff.o:
	$(C) $(CFLAGS) $(FNCSD)get_K_eff.cpp -o $(FNCSD)get_K_eff.o

get_F.o:
	$(C) $(CFLAGS) $(FNCSD)get_F.cpp -o $(FNCSD)get_F.o

disp.o:
	$(C) $(CFLAGS) $(FNCSD)disp.cpp -o $(FNCSD)disp.o

log.o:
	$(C) $(CFLAGS) $(FNCSD)log.cpp -o $(FNCSD)log.o

clean:
	rm main.o
	rm $(FNCSD)get_M.o
	rm $(FNCSD)get_K.o
	rm $(FNCSD)get_K_eff.o
	rm $(FNCSD)get_F.o
	rm $(FNCSD)disp.o
	rm $(FNCSD)log.o

task1:
	clear
	time ./app.out -L 10.0 -N_e 24 -A 0.012 -I 0.0000144 -E 210000000000.0 -eq 0

task2:
	clear
	time ./app.out -L 10.0 -N_e 24 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -Tl 1.0 -N_t 8000.0 -rho 7850.0 -eq 1 -sch 0
	
task3:
	clear
	time ./app.out -L 10.0 -N_e 24 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -Tl 1.0 -N_t 1000.0 -rho 7850.0 -eq 1 -sch 1

task4:
	clear
	time mpiexec -n 2 app.out -L 10.0 -N_e 24 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -Tl 1.0 -N_t 8000.0 -rho 7850.0 -eq 1 -sch 0 

task5:
	clear
	time mpiexec -n 2 app.out -L 10.0 -N_e 24 -A 0.012 -I 0.0000144 -E 210000000000.0 -T 1.0 -Tl 1.0 -N_t 1000.0 -rho 7850.0 -eq 1 -sch 1
