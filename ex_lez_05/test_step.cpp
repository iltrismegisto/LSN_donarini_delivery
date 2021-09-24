#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double S(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return pow(pow(M_PI,-0.5)*exp(-r),2);
}

double P(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return pow(pow(M_PI/2.,-0.5)/8.*r*exp(-r/2.)*z/r,2);
}

int main (int argc, char *argv[]){
  //Inizialization of the Random Number Generator
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

  int test = 10000;

  double delta_0 = 1.23;
  double delta_1 = 2.92;
  double delta_0_g = 0.76;
  double delta_1_g = 1.89;

  double *x_0 = new double [3];
  double *y_0 = new double [3];
  double *x_1 = new double [3];
  double *y_1 = new double [3];
  double *x_0_g = new double [3];
  double *y_0_g = new double [3];
  double *x_1_g = new double [3];
  double *y_1_g = new double [3];

  double ave_a0 = 0.;
  double ave_a1 = 0.;
  double ave_a0_g = 0.;
  double ave_a1_g = 0.;
  double a = 0.;

  x_0[0]=0.;  x_0_g[0]=0.;
  x_0[1]=0.;  x_0_g[1]=0.;
  x_0[2]=0.;  x_0_g[2]=0.;

  x_1[0]=0.;  x_1_g[0]=0.;
  x_1[1]=0.;  x_1_g[1]=0.;
  x_1[2]=2.;  x_1_g[2]=2.;

  for(size_t i=0;i<test;i++){
    for(size_t j=0;j<3;j++){
      y_0[j]=rnd.Rannyu(x_0[j]-delta_0,x_0[j]+delta_0);
      y_1[j]=rnd.Rannyu(x_1[j]-delta_1,x_1[j]+delta_1);
      y_0_g[j]=rnd.Gauss(x_0_g[j],delta_0_g);
      y_1_g[j]=rnd.Gauss(x_1_g[j],delta_1_g);
    }

    //Uniform Metropolis 1S
    a=min(1.,(S(y_0[0], y_0[1], y_0[2])/S(x_0[0], x_0[1], x_0[2])));
    if(rnd.Rannyu()<=a){
      for(size_t j=0;j<3;j++){
        x_0[j]=y_0[j];
      }
    }
    ave_a0+=a;

    //Gauss Metropolis 1S
    a=min(1.,(S(y_0_g[0], y_0_g[1], y_0_g[2])/S(x_0_g[0], x_0_g[1], x_0_g[2])));
    if(rnd.Rannyu()<=a){
      for(size_t j=0;j<3;j++){
        x_0_g[j]=y_0_g[j];
      }
    }
    ave_a0_g+=a;

    //Uniform Metropolis 2P
    a=min(1.,(P(y_1[0], y_1[1], y_1[2])/P(x_1[0], x_1[1], x_1[2])));
    if(rnd.Rannyu()<=a){
      for(size_t j=0;j<3;j++){
        x_1[j]=y_1[j];
      }
    }
    ave_a1+=a;

    //Gauss Metropolis 2P
    a=min(1.,(P(y_1_g[0], y_1_g[1], y_1_g[2])/P(x_1_g[0], x_1_g[1], x_1_g[2])));
    if(rnd.Rannyu()<=a){
      for(size_t j=0;j<3;j++){
        x_1_g[j]=y_1_g[j];
      }
    }
    ave_a1_g+=a;

  }

  ave_a0/=test/100.;
  ave_a0_g/=test/100.;
  ave_a1/=test/100.;
  ave_a1_g/=test/100.;

  cout << "S1 lin " << ave_a0 << endl;
  cout << "S1 gau " << ave_a0_g << endl;
  cout << "2P lin " << ave_a1 << endl;
  cout << "2P gau " << ave_a1_g << endl;

  return 23;
}
