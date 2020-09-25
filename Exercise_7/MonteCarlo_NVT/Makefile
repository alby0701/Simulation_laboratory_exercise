CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT
	rm -rf *block*
	rm -rf *gdir*
	rm -rf *equilibration*
	rm -rf *autocorrelation*
