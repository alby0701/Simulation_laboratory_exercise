CC = g++
CFLAGS = -Wall -O3 --std=c++17 -march=native
# CFLAGS = -Wall -g

main: main.cpp tsp_ga.o
	$(CC) tsp_ga.o main.cpp -o main $(CFLAGS)
    
all: main
tsp_ga.o: tsp_ga.cpp tsp_ga.h
	$(CC) -fPIC -c tsp_ga.cpp -o tsp_ga.o $(CFLAGS)
tsp_ga.so: tsp_ga.o
	$(CC) -shared tsp_ga.o -o tsp_ga.so $(CFLAGS)


clean :
	rm *.o *.so main
