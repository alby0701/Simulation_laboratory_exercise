CC = g++
CFLAGS = -Wall -O3 --std=c++17 -march=native
# CFLAGS = -Wall -g

main: main.cpp tsp.o
	$(CC) tsp.o main.cpp -o main $(CFLAGS)
    
all: main
tsp.o: tsp.cpp tsp.h
	$(CC) -fPIC -c tsp.cpp -o tsp.o $(CFLAGS)
tsp.so: tsp.o
	$(CC) -shared tsp.o -o tsp.so $(CFLAGS)


clean :
	rm *.o *.so main
