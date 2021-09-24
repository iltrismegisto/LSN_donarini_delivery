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

  int block_number = 100;
  int shot = 1000000;
  int block_size = shot/block_number;
  int intersection = 0;
  int ruler = 100;
  double d=1/double(ruler);
  double long_stick= 0.75*d;
  double x, y;
  double *pi = new double [block_number];
  double *pi_2 = new double [block_number];
  double *prog_pi = new double [block_number];
  double *prog_pi_2 = new double [block_number];
	double *prog_sigma = new double [block_number];

 for (size_t i = 0; i < block_number; i++) {
   pi[i]=0;
   pi_2[i]=0;
   prog_pi[i]=0;
   prog_pi_2[i]=0;
   prog_sigma[i]=0;
 }

  for (size_t i = 0; i < block_number; i++) {
    intersection = 0;
    for (size_t j = 0; j < block_size; j++) {
      x=rnd.Rannyu();
      y=rnd.Rannyu();
      if (pow(x,2) + pow(y,2) < 1) {
        if(d*rnd.Rannyu()<(long_stick*x/sqrt(x*x+y*y))) intersection++;
      } else j --;
    }
    pi[i]=2.*long_stick*double(block_size)/(double(intersection)*d);
    pi_2[i]=pow(pi[i],2);
  }

	ofstream out;
	out.open("pi.dat");

  for(size_t i = 1; i < block_number; i++) {
    for (size_t j = 0; j < (i+1); j++) {
      prog_pi[i]+=pi[j];
      prog_pi_2[i]+=pi_2[j];
    }
    prog_pi[i]  /=double(i+1);
    prog_pi_2[i]/=double(i+1);
    prog_sigma[i]=error(prog_pi, prog_pi_2, i);

    out << (i+1) << ";" << prog_pi[i] << ";" << prog_sigma[i] << endl;
  }

	out.close();
  return 23;

}
