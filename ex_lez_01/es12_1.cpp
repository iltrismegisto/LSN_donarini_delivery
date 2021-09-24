#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

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

	int realizations = 100000;
	int N;
	double sum_std=0.;
	double sum_exp=0.;
	double sum_lor=0.;
	double lambda=1.;
	double mu=0.;
	double gamma=1.;

  N=1;
	ofstream distr_std;
	ofstream distr_exp;
	ofstream distr_lor;
	distr_std.open("distr_std_N_1.dat");
	distr_exp.open("distr_exp_N_1.dat");
	distr_lor.open("distr_lor_N_1.dat");

	for (size_t i = 0; i < realizations; i++) {
	  distr_std << rnd.Rannyu() << endl;
	  distr_exp << rnd.Exp(lambda) << endl;
	  distr_lor << rnd.Lor(gamma, mu) << endl;
	}

	distr_std.close();
	distr_exp.close();
	distr_lor.close();

	N=2;
	distr_std.open("distr_std_N_2.dat");
	distr_exp.open("distr_exp_N_2.dat");
	distr_lor.open("distr_lor_N_2.dat");


	for(size_t i=0;i<realizations;i++){
		sum_std = 0.;
		sum_exp = 0.;
		sum_lor = 0.;
		for(size_t j=0;j<N;j++){
			sum_std+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_lor+=rnd.Lor(gamma,mu);
		}
		distr_std << sum_std/N << endl;
		distr_exp << sum_exp/N << endl;
		distr_lor << sum_lor/N << endl;

	}

	distr_std.close();
	distr_exp.close();
	distr_lor.close();

	N=10;
	distr_std.open("distr_std_N_10.dat");
	distr_exp.open("distr_exp_N_10.dat");
	distr_lor.open("distr_lor_N_10.dat");

	for(size_t i=0;i<realizations;i++){
		sum_std = 0.;
		sum_exp = 0.;
		sum_lor = 0.;
		for(size_t j=0;j<N;j++){
			sum_std+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_lor+=rnd.Lor(gamma,mu);
		}
		distr_std << sum_std/N << endl;
		distr_exp << sum_exp/N << endl;
		distr_lor << sum_lor/N << endl;

	}

	distr_std.close();
	distr_exp.close();
	distr_lor.close();

	N=100;
	distr_std.open("distr_std_N_100.dat");
	distr_exp.open("distr_exp_N_100.dat");
	distr_lor.open("distr_lor_N_100.dat");

	for(size_t i=0;i<realizations;i++){
		sum_std = 0.;
		sum_exp = 0.;
		sum_lor = 0.;
		for(size_t j=0;j<N;j++){
			sum_std+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_lor+=rnd.Lor(gamma,mu);
		}
		distr_std << sum_std/N << endl;
		distr_exp << sum_exp/N << endl;
		distr_lor << sum_lor/N << endl;
	}

	distr_std.close();
	distr_exp.close();
	distr_lor.close();

	return 23;

}
