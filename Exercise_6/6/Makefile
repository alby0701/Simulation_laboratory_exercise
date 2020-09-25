CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
main : main.o
	$(CC)  main.o funzioni.o random.o -o main
main.o : main.cpp funzioni.o random.o
	$(CC) random.o funzioni.o -c main.cpp -o main.o $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o: funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o Monte_Carlo_ISING_1D main output*.* seed.out config.final
