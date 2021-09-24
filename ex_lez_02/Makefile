CC = g++
CFLAGS = -Wall -O3 --std=c++11

es21_1.exe : es21_1.o random.o
	$(CC) random.o es21_1.o -o es21_1.exe
es21_1.o : es21_1.cpp
	$(CC) -c es21_1.cpp -o es21_1.o $(CFLAGS)
es21_2.exe : es21_2.o random.o
	$(CC) random.o es21_2.o -o es21_2.exe
es21_2.o : es21_2.cpp
	$(CC) -c es21_2.cpp -o es21_2.o $(CFLAGS)
es22_1.exe : es22_1.o random.o
	$(CC) random.o es22_1.o -o es22_1.exe
es22_1.o : es22_1.cpp
	$(CC) -c es22_1.cpp -o es22_1.o $(CFLAGS)
es22_2.exe : es22_2.o random.o
	$(CC) random.o es22_2.o -o es22_2.exe
es22_2.o : es22_2.cpp
	$(CC) -c es22_2.cpp -o es22_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es22_1.exe es21_1.exe seed.out int_unif.dat int_import.dat

all : es21_1.exe es21_2.exe es22_1.exe es22_2.exe

execute : ./es21_1.exe ./es21_2.exe ./es22_1.exe ./es22_2.exe
