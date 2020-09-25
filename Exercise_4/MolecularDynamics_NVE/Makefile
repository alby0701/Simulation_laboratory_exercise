CC = g++
#CC = clang++
CFLAGS = -Wall -O3 --std=c++17 -march=native #-ffast-math

MolDyn_NVE: MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) MolDyn_NVE.cpp -o MolDyn_NVE $(CFLAGS)
clean:
	rm MolDyn_NVE
	rm -rf output.dat
	rm -rf frames/*.xyz
