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

	double imp=0;
  double *im = new double [block_number];
  double *im2 = new double [block_number];
  double *sig_im = new double [block_number];
  double *sum_pr_im = new double [block_number];
  double *sum2_pr_im = new double [block_number];

  ofstream import_sample;
  import_sample.open("int_import.dat");

	for (size_t i = 0; i < block_number; i++) {
		im[i] = 0;
	  im2[i] = 0;
		sig_im[i] = 0;
		sum_pr_im[i] = 0;
		sum2_pr_im[i] = 0;

		for (size_t j = 1; j < block_size+1; j++) {
			imp=(1-sqrt(1-rnd.Rannyu()));
			im[i] += (pi*0.5*cos(pi*0.5*imp))/(2*(1-imp));

		}
	  im[i]/=block_size;
	  im2[i]=pow(im[i],2);

	  for (size_t j = 0; j < i+1; j++) {
			sum_pr_im[i] += im[j];
			sum2_pr_im[i]+= im2[j];
	  }
	  sum_pr_im[i] /= (i+1);
	  sum2_pr_im[i]/= (i+1);
		sig_im[i]=error(sum_pr_im,sum2_pr_im,i);
		import_sample << (i+1)*block_size << ";" << sum_pr_im[i] << ";" << sig_im[i] << endl;

	 }

   import_sample.close();

   return 23;

 }
