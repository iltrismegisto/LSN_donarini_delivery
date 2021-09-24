#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props;
int iu;
int ic;
int im;
int ix;
int ig;

double nbins;
double walker[m_props];

// averages
double blk_av[m_props];
double blk_norm;
double accepted;
double attempted;

double glob_av[m_props];
double glob_av2[m_props];

double stima_u;
double stima_c;
double stima_m;
double stima_x;
double stima_g;

double err_u;
double err_c;
double err_m;
double err_x;
double err_g;

//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
int restart;
double beta;
double temp;
double J;
double h;

// simulation
int nstep;
int nblk;
int metro;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
