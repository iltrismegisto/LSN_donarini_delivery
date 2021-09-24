#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

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

	int M = 2.5*1E4;
	int block_number = 100;
	int block_size = M/block_number;
	double s_0 = 100.;
	double f_t = 1;
	double strike_price = 100;
	double inter_rate = 0.1;
	double volatility = 0.25;
	double d = 0.;

//call direct
	double *c_d = new double [block_number];
	double *c_d2 = new double [block_number];
	double *c_d_sig = new double [block_number];
	double *c_d_sum = new double [block_number];
	double *c_d_sum2 = new double [block_number];
//put direct
	double *p_d = new double [block_number];
	double *p_d2 = new double [block_number];
	double *p_d_sig = new double [block_number];
	double *p_d_sum = new double [block_number];
	double *p_d_sum2 = new double [block_number];

	ofstream direct;
	direct.open("direct.dat");

	for (size_t i = 0; i < block_number; i++) {

		c_d[i]=0.;
		c_d2[i]=0.;
		c_d_sig[i]=0.;
		c_d_sum[i]=0.;
		c_d_sum2[i]=0.;

		p_d[i]=0.;
		p_d2[i]=0.;
		p_d_sig[i]=0.;
		p_d_sum[i]=0.;
		p_d_sum2[i]=0.;

		for (size_t j = 0; j < block_size; j++) {

			d = 0.;
			d = s_0*exp((inter_rate-0.5*pow(volatility,2))+volatility*rnd.Gauss(0., 1.));

			if ((d - strike_price) > 0) {
				c_d[i]+=exp(- inter_rate*f_t)*(d - strike_price);
			} else {
				c_d[i]+=0;
			}
			if ((strike_price - d) > 0) {
				p_d[i]+=exp(- inter_rate*f_t)*(strike_price - d);
			} else {
				p_d[i]+=0;
			}

		}

		p_d[i]/=double(block_size);
		p_d2[i]=pow(p_d[i],2.);
		c_d[i]/=double(block_size);
		c_d2[i]=pow(c_d[i],2.);

		for(size_t k = 0 ; k < i+1 ; k++){

			p_d_sum[i] +=p_d[k];
			p_d_sum2[i]+=p_d2[k];
			c_d_sum[i] +=c_d[k];
			c_d_sum2[i]+=c_d2[k];

		}

		c_d_sum[i] /=double(i+1);
		c_d_sum2[i]/=double(i+1);
		p_d_sum[i] /=double(i+1);
		p_d_sum2[i]/=double(i+1);

		p_d_sig[i]=error(p_d_sum, p_d_sum2, i);
		c_d_sig[i]=error(c_d_sum, c_d_sum2, i);

		direct << i+1 << ";" << p_d_sum[i] << ";" << p_d_sig[i] << ";" << c_d_sum[i] << ";" << c_d_sig[i] << endl;

	}

	direct.close();

  return 23;

}
