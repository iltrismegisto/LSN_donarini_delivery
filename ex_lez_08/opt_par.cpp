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

  int M = 100000;
  int block_number = 100;
  int block_size = M/block_number;

  double attempts=0.;
  double accepted=0.;
  double x = 0.;
  double ratio;

  double step = 2.5;
  double mu = 0;
  double sigma = 0;
  double mu_test;
  double sigma_test;
  double mean = 0.;
  double x_new = 0.;

  double min_energy = 99999.0;

  cout << "-------- Optimizing parameters mu and sigma --------" << endl;

  for ( mu_test = 0.75 ; mu_test <= 0.85 ; mu_test += 0.01) {
    for ( sigma_test = 0.55 ; sigma_test <= 0.65 ; sigma_test += 0.01) {
      attempts=0.;
      accepted=0.;

      for(int i=0 ; i<M ; i++) {
        x_new = x + (rnd.Rannyu()-0.5)*2*step;
        ratio = pow( psi(x_new, mu_test, sigma_test), 2) / pow( psi(x, mu_test, sigma_test) , 2);
        if(min(1.,ratio)>rnd.Rannyu()) {
          x = x_new;
          accepted++;
        }
        attempts++;
        mean += (-psi_2(x,mu_test,sigma_test)*0.5 + v(x)*psi(x,mu_test,sigma_test))/(psi(x,mu_test,sigma_test));
      }

      mean /= double(M);
      cout << "Energy : " << mean << endl;
      if(mean<min_energy) {
        min_energy=mean;
        mu = mu_test;
        sigma = sigma_test;
      }
      cout << "Acceptance : " << double(accepted)/attempts*100 << " %"<< endl;
    }
  }

  cout << "Best values for " << endl;
  cout << "   Energy: " << min_energy << endl;
  cout << "   Mu: " << mu << endl;
  cout << "   Sigma: " << sigma << endl;

  return 23;
}
