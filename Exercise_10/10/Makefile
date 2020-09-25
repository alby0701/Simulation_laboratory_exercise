CC = mpicxx
CFLAGS = -Wall -O3 --std=c++17 -march=native
# CFLAGS = -Wall -g

main: main.cpp tsp_ga.o
	$(CC) tsp_ga.o main.cpp -o main $(CFLAGS)
    
all: main
tsp_ga.o: tsp_ga.cpp tsp_ga.h
	$(CC) -c tsp_ga.cpp -o tsp_ga.o $(CFLAGS)

clean :
	rm *.o main
	
run:
	module load mpi/openmpi-x86_64
	mpiexec -np 4 main cities config
