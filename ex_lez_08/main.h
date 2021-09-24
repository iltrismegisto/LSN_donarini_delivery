#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double psi(double x, double mu, double sigma){
  double den = 2*sigma*sigma;
  double r = exp(-pow(x-mu,2)/den) + exp(-pow(x+mu,2)/den);
  return r;
}

double v(double x){
  return pow(x,4.)-2.5*pow(x,2.);
}

double error(double av, double av2, int n){
  double err =0;
  if(n==0)
    return 0.;
  else {
    err=sqrt(1./n*(av2-av*av));
  }
  return err;
}

double psi_2(double x, double mu, double sigma){
  double den = 2*sigma*sigma;
  double b = 1./pow(sigma,2);
  double a1 = pow((x-mu)/(sigma*sigma),2);
  double a2 = pow((x+mu)/(sigma*sigma),2);
  return exp(-pow(x-mu,2)/den)*(a1-b) + exp(-pow(x+mu,2)/den)*(a2-b);
}
