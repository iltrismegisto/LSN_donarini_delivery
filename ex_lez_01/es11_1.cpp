#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(ave[n],2))/n);
};

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	int M = 1000000;
	int block_number = 100;
	int block_size = M/block_number;
	double sum = 0.;
	double *ave = new double [block_number];
	double *ave2 = new double [block_number];
	double *sum_prog = new double [block_number];
	double *sum2_prog = new double [block_number];
	double *err_prog = new double [block_number];

	for (size_t i = 0; i < block_number; i++) {
		ave[i] = 0.;
		ave2[i] = 0.;
		sum_prog[i] = 0.;
		sum2_prog[i] = 0.;
		err_prog[i] = 0.;
	}

	for (size_t i = 0; i < block_number; i++) {
		sum = 0.;
		for (size_t j = 0; j < block_size; j++) {
			sum+=rnd.Rannyu();
		}
		ave[i] = sum/block_size;
		ave2[i]= pow(ave[i], 2);
	}

	for (size_t i = 0; i < block_number; i++) {
		for (size_t j = 0; j < (i+1); j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i]+= ave2[j];
		}
		sum_prog[i] = sum_prog[i]/double(i+1);
		sum2_prog[i]= sum2_prog[i]/double(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
	}

	ofstream out;
	out.open("res1_1.dat");

	for(int i = 0; i < block_number; i++){
		out << (i+1)*block_size << ";" << (sum_prog[i]-0.5) << ";" << err_prog[i] << endl;
	}

	out.close();

	return 23;

}
