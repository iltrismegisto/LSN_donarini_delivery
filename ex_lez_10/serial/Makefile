CC = g++
CFLAGS = -Wall -O3 --std=c++11

sim_ann.exe : sim_ann.o random.o
	$(CC) random.o sim_ann.o -o sim_ann.exe
sim_ann.o : sim_ann.cpp
	$(CC) -c sim_ann.cpp -o sim_ann.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o sim_ann.exe *.gene evolution.dat
