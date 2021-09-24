CC = g++
CFLAGS = -Wall -O3 --std=c++11

test_step.exe : test_step.o random.o
	$(CC) random.o test_step.o -o test_step.exe
test_step.o : test_step.cpp
	$(CC) -c test_step.cpp -o test_step.o $(CFLAGS)

ex_51_1.exe : ex_51_1.o random.o
	$(CC) random.o ex_51_1.o -o ex_51_1.exe
ex_51_1.o : ex_51_1.cpp
	$(CC) -c ex_51_1.cpp -o ex_51_1.o $(CFLAGS)

ex_51_2.exe : ex_51_2.o random.o
	$(CC) random.o ex_51_2.o -o ex_51_2.exe
ex_51_2.o : ex_51_2.cpp
	$(CC) -c ex_51_2.cpp -o ex_51_2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out test_step.exe ex_51_1.exe ex_51_2.exe

all: test_step.exe ex_51_1.exe ex_51_2.exe
