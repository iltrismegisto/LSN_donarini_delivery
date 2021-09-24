#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

const double pi = M_PI;

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0.;
	else
		return sqrt((av2[n]-pow(ave[n],2))/double(n));
}

int main(int argc, char const *argv[]) {

  Random rnd;
	int seed[23];
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

  int M = 10000000;
  int block_number = 100;
  int block_size = M/block_number;

  double *un = new double [block_number];
  double *un2 = new double [block_number];
	double *sig_un = new double [block_number];
	double *sum_pr_un = new double [block_number];
	double *sum2_pr_un = new double [block_number];

  ofstream unif_sample;
  unif_sample.open("int_unif.dat");

	for (size_t i = 0; i < block_number; i++) {
		un[i] = 0;
	  un2[i]= 0;
		sig_un[i] = 0;
		sum_pr_un[i] = 0;
		sum2_pr_un[i] = 0;

		for (size_t j = 0; j < block_size; j++) {
			un[i] += pi*0.5*cos(pi*(rnd.Rannyu()*0.5));
		}
    un[i]/=block_size;
    un2[i]=pow(un[i],2);

    for (int j = 0; j < i+1; j++) {
      sum_pr_un[i]+=un[j];
      sum2_pr_un[i]+= un2[j];
    }
    sum_pr_un[i]/=double(i+1);
    sum2_pr_un[i]/=double(i+1);
		sig_un[i]=error(sum_pr_un,sum2_pr_un,i);
		unif_sample << (i+1)*block_size << ";" << sum_pr_un[i] << ";" << sig_un[i] << endl;

  }

	unif_sample.close();

	return 23;
}
