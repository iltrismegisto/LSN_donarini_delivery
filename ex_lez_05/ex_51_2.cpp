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

double S(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return pow(pow(M_PI,-0.5)*exp(-r),2);
}

double P(double x, double y, double z){
  double r = sqrt(x*x + y*y + z*z);
  return pow(pow(M_PI/2.,-0.5)/8.*r*exp(-r/2.)*z/r,2);
}

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

  // Metropolis Uniform

  int M = 1E6;
  int block_number = 100;
  int block_size = M/block_number;

  double delta_0 = 0.76;
  double delta_1 = 1.89;

  double *x_0 = new double [3];
  double *y_0 = new double [3];
  double *x_1 = new double [3];
  double *y_1 = new double [3];

  double *ave_r0 = new double [block_number];
  double *ave_r0_2 = new double [block_number];
  double *ave_r1 = new double [block_number];
  double *ave_r1_2 = new double [block_number];
  double *prog0 = new double [block_number];
  double *prog1 = new double [block_number];
  double *prog0_2 = new double [block_number];
  double *prog1_2 = new double [block_number];
  double *s_r0 = new double [block_number];
  double *s_r1 = new double [block_number];

  double a = 0.;

  x_0[0]=0.;
  x_0[1]=0.;
  x_0[2]=0.;

  x_1[0]=0.;
  x_1[1]=0.;
  x_1[2]=2.;

  ofstream s_1, p_2, s_1_ave, p_2_ave;
  s_1.open("1s_g.dat");
	p_2.open("2p_g.dat");
  s_1_ave.open("1s_g_ave.dat");
	p_2_ave.open("2p_g_ave.dat");

  for(size_t i=0;i<block_number;i++){
  	for(size_t j=0;j<block_size;j++){
  		for(size_t j=0;j<3;j++){
        y_0[j]=rnd.Gauss(x_0[j],delta_0);
				y_1[j]=rnd.Gauss(x_1[j],delta_1);
			}
      //Gauss Metropolis 1S
      a=min(1.,(S(y_0[0], y_0[1], y_0[2])/S(x_0[0], x_0[1], x_0[2])));
      if(rnd.Rannyu()<=a){
        for(size_t j=0;j<3;j++){
          x_0[j]=y_0[j];
        }
      }
      ave_r0[i]+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);
      s_1 << x_0[0] << ";" << x_0[1] << ";" << x_0[2] << endl;

      //Gauss Metropolis 2P
      a=min(1.,(P(y_1[0], y_1[1], y_1[2])/P(x_1[0], x_1[1], x_1[2])));
      if(rnd.Rannyu()<=a){
        for(size_t j=0;j<3;j++){
          x_1[j]=y_1[j];
        }
      }
      ave_r1[i]+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);
      p_2 << x_1[0] << ";" << x_1[1] << ";" << x_1[2] << endl;
    }

   	ave_r0[i]/=double(block_size);
   	ave_r0_2[i]=pow(ave_r0[i],2.);
   	ave_r1[i]/=double(block_size);
   	ave_r1_2[i]=pow(ave_r1[i],2.);

    for(size_t k=0;k<i+1;k++){
      prog0[i]+=ave_r0[k];
      prog1[i]+=ave_r1[k];
      prog0_2[i]+=ave_r0_2[k];
      prog1_2[i]+=ave_r1_2[k];
		}

		prog0[i]/=double(i+1);
		prog1[i]/=double(i+1);
		prog0_2[i]/=double(i+1);
		prog1_2[i]/=double(i+1);

		if(i==0){
			s_r0[i]=0.;
			s_r1[i]=0.;
		}

		else{
			s_r0[i]=sqrt((prog0_2[i]-pow(prog0[i],2))/double(i));
			s_r1[i]=sqrt((prog1_2[i]-pow(prog1[i],2))/double(i));
		}

		s_1_ave << i+1 << ";" << prog0[i] << ";" << s_r0[i] << endl;
		p_2_ave << i+1 << ";" << prog1[i] << ";" << s_r1[i] << endl;
  }

  s_1.close();
 	p_2.close();
 	s_1_ave.close();
 	p_2_ave.close();

  return 23;

}
