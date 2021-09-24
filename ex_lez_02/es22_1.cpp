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

  ofstream discrete;
  discrete.open("rnd_wlk_disc.dat");
	ofstream discrete_step;
	discrete_step.open("rnd_wlk_disc_step.dat");

  int M=10000;
  int step=100;
  int k;

  int *moves = new int [3];
  double *dis = new double [M];
  double *dis2 = new double [M];
  double *sum_pr_dis = new double [M];
	double *sum2_pr_dis = new double [M];
  double *sig_dis = new double [M];

  for (size_t i = 0; i < M; i++) {

    moves[0]=0;
    moves[1]=0;
    moves[2]=0;

		for (size_t j = 0; j < step; j++) {
			k=rnd.Rannyu(0,6);
      if(k==0) moves[0]++;
			if(k==1) moves[0]--;
			if(k==2) moves[1]++;
			if(k==3) moves[1]--;
			if(k==4) moves[2]++;
			if(k==5) moves[2]--;
      dis[j]+=sqrt(pow(moves[0],2)+pow(moves[1],2)+pow(moves[2],2));
			if(i==0){
      	discrete_step << j+1 << ";"<<dis[j] << endl;
      }
    }
  }

  for (size_t i = 0; i < step; i++) {
    dis[i]/=M;
    dis2[i]=pow(dis[i],2);
    for (int j = 0; j < i+1; j++) {
      sum_pr_dis[i]+= dis[j];
      sum2_pr_dis[i]+= dis2[j];
    }
    sum_pr_dis[i]/=double(i+1);
    sum2_pr_dis[i]/=double(i+1);
    sig_dis[i]=error(sum_pr_dis,sum2_pr_dis,i);
    discrete << (i+1) << ";" << dis[i] << ";" << sig_dis[i] << endl;
  }

   discrete.close();
	 discrete_step.close();

   return 23;

  }
