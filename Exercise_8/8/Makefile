CC = g++
CFLAGS = -Wall -O3 --std=c++17 -march=native

main: main.cpp    
	$(CC) main.cpp random.cpp -o main $(CFLAGS)
    
all: main
random.o: random.cpp random.h
	$(CC) -fPIC -c random.cpp -o random.o $(CFLAGS)
random.so: random.o
	$(CC) -shared random.o -o random.so $(CFLAGS)

#main: main.cpp random.o
#s	$(CC) random.o -c main.cpp -o main $(CFLAGS)

clean :
	rm *.o *.so *.dat main
