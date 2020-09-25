CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: gbm gbm_d

gbm : gbm.o random.o
	$(CC)  gbm.o random.o -o gbm
gbm_d : gbm_d.o random.o
		$(CC)  gbm_d.o random.o -o gbm_d
gbm.o : gbm.cpp
	$(CC) -c gbm.cpp -o gbm.o $(CFLAGS)
gbm_d.o : gbm_d.cpp
		$(CC) -c gbm_d.cpp -o gbm_d.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o gbm gbm_d seed.out gbm.dat gbm_d.dat
