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
	rnd.SaveSeed();

	int num_sim = 100;
	int num_throws = 1000000;
	int expect = double(num_throws)/num_sim;
	int l = 0;
	double chi = 0;
	double *hits = new double [int(expect)];

	for (size_t i = 0; i < num_sim; i++) {
		hits[i]=0;
	}

	ofstream out;
	out.open("res1_3.dat");

	for (size_t i = 0; i < num_sim; i++) {
		l = 0;
		chi = 0;

		for (size_t j = 0; j < num_throws; j++){
			l = int(100*rnd.Rannyu());
			hits[l] += 1;
		}

		for (size_t j = 0; j < num_sim; j++){
			chi += pow(hits[j] - expect, 2) / expect;
			hits[j]=0;
		}
		out << i+1 << ";" << chi <<endl;
	}

	out.close();

	return 23;

}
