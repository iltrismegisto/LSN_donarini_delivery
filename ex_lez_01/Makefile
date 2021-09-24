CC = g++
CFLAGS = -Wall -O3 --std=c++11

es11_1.exe : es11_1.o random.o
	$(CC) random.o es11_1.o -o es11_1.exe
es11_1.o : es11_1.cpp
	$(CC) -c es11_1.cpp -o es11_1.o $(CFLAGS)

es11_2.exe : es11_2.o random.o
	$(CC) random.o es11_2.o -o es11_2.exe
es11_2.o : es11_2.cpp
	$(CC) -c es11_2.cpp -o es11_2.o $(CFLAGS)

es11_3.exe : es11_3.o random.o
	$(CC) random.o es11_3.o -o es11_3.exe
es11_3.o : es11_3.cpp
	$(CC) -c es11_3.cpp -o es11_3.o $(CFLAGS)

es12_1.exe : es12_1.o random.o
	$(CC) random.o es12_1.o -o es12_1.exe
es12_1.o : es12_1.cpp
	$(CC) -c es12_1.cpp -o es12_1.o $(CFLAGS)

es13_1.exe : es13_1.o random.o
	$(CC) random.o es13_1.o -o es13_1.exe
es13_1.o : es13_1.cpp
	$(CC) -c es13_1.cpp -o es13_1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es11_1.exe es11_2.exe es11_3.exe res1_1.dat res1_2.dat res1_3.dat es12_1.exe distr_std_N_1.dat distr_exp_N_1.dat distr_lor_N_1.dat distr_std_N_2.dat distr_exp_N_2.dat distr_lor_N_2.dat distr_std_N_10.dat distr_exp_N_10.dat distr_lor_N_10.dat distr_std_N_100.dat distr_exp_N_100.dat distr_lor_N_100.dat es13_1.exe res3_1.dat seed.out

all: es11_1.exe es11_2.exe es11_3.exe es12_1.exe es13_1.exe
