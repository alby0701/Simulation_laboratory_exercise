CC = g++
#CFLAGS = -Wall -O3 --std=c++17 -march=native
CFLAGS = -Wall -g

main: main.cpp tsp_sa.o
	$(CC) tsp_sa.o main.cpp -o main $(CFLAGS)

all: main
tsp_sa.o: tsp_sa.cpp tsp_sa.h
	$(CC) -c tsp_sa.cpp -o tsp_sa.o $(CFLAGS)
	
clean :
	rm *.o main
