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

  ofstream continuum;
  continuum.open("rnd_wlk_continuum.dat");
	ofstream continuum_step;
	continuum_step.open("rnd_wlk_continuum_step.dat");

  int M=10000;
  int step=100;
  double *moves = new double [3];

  double phi;
  double theta;
  double *con = new double [M];
  double *con2 = new double [M];
  double *sum_pr_con = new double [M];
	double *sum2_pr_con = new double [M];
  double *sig_con = new double [M];

  for (size_t i = 0; i < M; i++) {

    moves[0]=0;
    moves[1]=0;
    moves[2]=0;

		for (size_t j = 0; j < step; j++) {
			phi=rnd.Rannyu(0,2*pi);
      theta=acos(1-2*rnd.Rannyu());
      moves[0]+= sin(theta)*cos(phi);
			moves[1]+= sin(theta)*sin(phi);
			moves[2]+= cos(theta);
      con[j]+=sqrt(pow(moves[0],2)+pow(moves[1],2)+pow(moves[2],2));
			if(i==0){
	      continuum_step << j+1 << ";"<<con[j] << endl;
	    }
    }
	}

  for (size_t i = 0; i < step; i++) {
    con[i]/=M;
    con2[i]=pow(con[i],2);
    for (int j = 0; j < i+1; j++) {
      sum_pr_con[i]+= con[j];
      sum2_pr_con[i]+= con2[j];
    }
    sum_pr_con[i]/=double(i+1);
    sum2_pr_con[i]/=double(i+1);
    sig_con[i]=error(sum_pr_con,sum2_pr_con,i);
    continuum << (i+1) << ";" << con[i] << ";" << sig_con[i] << endl;
  }

   continuum.close();
	 continuum_step.close();

   return 23;

  }
