CC = g++
CFLAGS = -Wall -O3 --std=c++11

integral : integral.o random.o
	$(CC)  integral.o random.o -o integral
integral.o : integral.cpp
		$(CC) -c integral.cpp -o integral.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o integral seed.out
