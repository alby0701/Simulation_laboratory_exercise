CC = g++
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o
	$(CC) random.o main.o -o main
Exponential : Exponential.o random.o
		$(CC) random.o Exponential.o -o Exponential
Lorentzian : Lorentzian.o random.o
			$(CC) random.o Lorentzian.o -o Lorentzian
1_3 : 1_3.o random.o
			$(CC) random.o 1_3.o -o 1_3
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
1_3.o : 1_3.cpp
		$(CC) -c 1_3.cpp -o 1_3.o $(CFLAGS)
Lorentzian.o : Lorentzian.cpp
		$(CC) -c Lorentzian.cpp -o Lorentzian.o $(CFLAGS)
Exponential.o : Exponential.cpp
			$(CC) -c Exponential.cpp -o Exponential.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main 1_3 Exponential Lorentzian seed.out
