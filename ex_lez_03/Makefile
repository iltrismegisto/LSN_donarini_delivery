CC = g++
CFLAGS = -Wall -O3 --std=c++11

es31_1.exe : es31_1.o random.o
	$(CC) random.o es31_1.o -o es31_1.exe
es31_1.o : es31_1.cpp
	$(CC) -c es31_1.cpp -o es31_1.o $(CFLAGS)
es32_1.exe : es32_1.o random.o
	$(CC) random.o es32_1.o -o es32_1.exe
es32_1.o : es32_1.cpp
	$(CC) -c es32_1.cpp -o es32_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es31_1.exe es32_1.exe seed.out
all: es31_1.exe es32_1.exe
