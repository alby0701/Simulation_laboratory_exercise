CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o
	$(CC)  main.o funzioni.o random.o -o main
rannyu: rannyu.o
	$(CC)  rannyu.o random.o -o rannyu
main.o : main.cpp funzioni.o random.o
	$(CC) random.o funzioni.o -c main.cpp -o main.o $(CFLAGS)

rannyu.o : rannyu.cpp random.o
		$(CC) random.o -c rannyu.cpp -o rannyu.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

funzioni.o: funzioni.cpp funzioni.hh
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)


clean :
	rm *.o main rannyu seed.out
