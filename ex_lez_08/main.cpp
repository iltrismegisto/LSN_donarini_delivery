#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "main.h"

using namespace std;

double psi(double x, double mu, double sigma);
double psi_2(double x, double mu, double sigma);
double v(double x);
double error(double av, double av2, int n);
double Psi_2(double x, double mu, double sigma);

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

  unsigned int M = 100000;
  unsigned int block_number = 100;
  unsigned int block_size = M/block_number;

  double step = 2.6;
  double ratio;
  double mean;
  double mu = 0.8;
  double sigma = 0.62;
  double x = 0.;
  double x_new;
  int attempts=0;
  int accepted=0;

  double *wf = new double [block_number];
  double *wf_2 = new double [block_number];
  double *wf_sum = new double [block_number];
  double *wf_sum_2 = new double [block_number];
  double *wf_sig = new double [block_number];

  int nbin = 300;
  double x_min=-5., x_max=5.;
	double bin_size = (x_max-x_min)/double(nbin);
  int index;
  int hist[nbin][block_number];

  for (size_t i=0 ; i<block_number ; i++) {
    wf[i]=0;
    wf_2[i]=0;
    wf_sum[i]=0;
    wf_sum_2[i]=0;
    wf_sig[i]=0;
  }

  ofstream pot("pot.dat");
  ofstream histo("hist.dat");

  for (size_t j = 0 ; j < block_number ; j++) {
    mean = 0;
    x = 0.;
    for (size_t i = 0 ; i < block_size ; i++) {
      x_new = x + (rnd.Rannyu()-0.5)*2*step;
      ratio = pow(psi(x_new, mu, sigma),2)/pow(psi(x, mu, sigma),2);
      if (min(1.,ratio) == 1 || min(1.,ratio)>rnd.Rannyu()){
        x = x_new;
        accepted++;
      }
      attempts++;
      mean += (-psi_2(x,mu,sigma)*0.5 + v(x)*psi(x,mu,sigma))/(psi(x,mu,sigma));

      index = int((x+x_max)/bin_size);
      hist[index][j]++;
    }
    cout << "Acceptance : " << double(accepted)/double(attempts)*100 << " %"<< endl;
    wf[j] = mean/double(block_size);
    wf_2[j] = wf[j]*wf[j];
  }

  for (size_t i=0 ; i<block_number ; i++ ) {
    for (size_t j=0 ; j<(i+1) ; j++ ) {
      wf_sum[i] += wf[j];
      wf_sum_2[i] += wf_2[j];
    }
    wf_sum[i] /= double(i+1);
    wf_sum_2[i] /= double(i+1);
    wf_sig[i]=error(wf_sum[i], wf_sum_2[i],i);

    pot <<(i+1)<<";"<< wf_sum[i] <<";"<< wf_sig[i] << endl;
  }
  cout << "Energy : " <<  wf_sum[block_number-1]<< endl;
  pot.close();
  pot.clear();

  cout << "HISTOGRAM DATA BLOCKING" << endl << endl;

  for(size_t ibin=0 ; ibin<nbin ; ibin++) {

    for(size_t i=0 ; i < block_number ; i++) {
      wf_sum[i] = 0;
      wf_sum_2[i] = 0;
      wf_sig[i] = 0;
    }

    for(size_t i=0 ; i<block_number ; i++) {

      for(size_t j=0 ; j<i+1 ; j++) {
        wf_sum[i] += hist[ibin][j];
        wf_sum_2[i] += hist[ibin][j]*hist[ibin][j];
      }

      wf_sum[i] /= double(i+1);
      wf_sum_2[i] /= double(i+1);
      wf_sig[i] = error(wf_sum[i],wf_sum_2[i],i);

      if(i == block_number-1){
        histo << (x_min+bin_size)+(ibin*bin_size) << ";" << wf_sum[i]/(double(block_size)*bin_size) << ";" << wf_sig[i]/(double(block_size)*bin_size) << endl;
      }
    }
  }
  histo.close();
  histo.clear();

  cout << "Done" << endl;

  return 23;

}
