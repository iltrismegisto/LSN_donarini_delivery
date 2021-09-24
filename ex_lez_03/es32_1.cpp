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
  int time_int = 100;
  double s_0 = 100.;
  double f_t = 1;
  double strike_price = 100;
  double inter_rate = 0.1;
  double volatility = 0.25;
  double t_step = f_t/double(time_int);
  double *intervals = new double [time_int+1];
  intervals[0]=s_0;

  //call time intervals
	double *c_i = new double [block_number];
	double *c_i2 = new double [block_number];
	double *c_i_sig = new double [block_number];
	double *c_i_sum = new double [block_number];
	double *c_i_sum2 = new double [block_number];
  //put time intervals
  double *p_i = new double [block_number];
  double *p_i2 = new double [block_number];
  double *p_i_sig = new double [block_number];
  double *p_i_sum = new double [block_number];
  double *p_i_sum2 = new double [block_number];

  ofstream interv;
  interv.open("intervals.dat");

  for (size_t i = 0; i < block_number; i++) {

		c_i[i]=0.;
		c_i2[i]=0.;
		c_i_sig[i]=0.;
		c_i_sum[i]=0.;
		c_i_sum2[i]=0.;

		p_i[i]=0.;
		p_i2[i]=0.;
		p_i_sig[i]=0.;
		p_i_sum[i]=0.;
		p_i_sum2[i]=0.;

		for (size_t j = 0; j < block_size; j++) {

			for (size_t k = 0; k < time_int; k++) {
				intervals[k+1] = intervals[k]*exp((inter_rate-0.5*pow(volatility,2.))*t_step + volatility*rnd.Gauss(0.,1.)*sqrt(t_step));
			}

			if ((intervals[time_int] - strike_price) > 0) {
				c_i[i]+=exp(- inter_rate*f_t)*(intervals[time_int] - strike_price);
			} else {
				c_i[i]+=0;
			}
			if ((strike_price - intervals[time_int]) > 0) {
				p_i[i]+=exp(- inter_rate*f_t)*(strike_price - intervals[time_int]);
			} else {
				p_i[i]+=0;
			}

		}

		p_i[i]/=double(block_size);
		p_i2[i]=pow(p_i[i],2.);
		c_i[i]/=double(block_size);
		c_i2[i]=pow(c_i[i],2.);

		for(size_t k = 0 ; k < i+1 ; k++){

      c_i_sum[i] += c_i[k];
			c_i_sum2[i]+= c_i2[k];
			p_i_sum[i] += p_i[k];
			p_i_sum2[i]+= p_i2[k];

		}

		p_i_sum[i] /=double(i+1);
		p_i_sum2[i]/=double(i+1);
    c_i_sum[i] /=double(i+1);
		c_i_sum2[i]/=double(i+1);

		p_i_sig[i]=error(p_i_sum, p_i_sum2, i);
		c_i_sig[i]=error(c_i_sum, c_i_sum2, i);

		interv << i+1 << ";" << p_i_sum[i] << ";" << p_i_sig[i] << ";" << c_i_sum[i] << ";" << c_i_sig[i] << endl;

	}

	interv.close();

  return 23;

}
