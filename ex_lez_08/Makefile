CC = g++
CFLAGS = -Wall -O3 --std=c++11

opt_par.exe : opt_par.o random.o
	$(CC) random.o opt_par.o -o opt_par.exe
opt_par.o : opt_par.cpp main.h
	$(CC) -c opt_par.cpp -o opt_par.o $(CFLAGS)
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out opt_par.exe *.dat
all : opt_par.exe main.exe
