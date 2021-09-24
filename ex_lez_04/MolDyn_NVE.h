#include <string>
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;
std::string phase;

// averages
double acc,att;

//restart
bool restart = false;
bool real_measure = false;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint;
double delta;

//Data Blocking
int block_number = 100.;
int block_size, h, c;
double *ave_pot = new double [block_number];
double *ave_kin = new double [block_number];
double *ave_tot = new double [block_number];
double *ave_t = new double [block_number];
double *ave_pot_2 = new double [block_number];
double *ave_kin_2 = new double [block_number];
double *ave_tot_2 = new double [block_number];
double *ave_t_2 = new double [block_number];
double *prog_pot = new double [block_number];
double *prog_kin = new double [block_number];
double *prog_tot = new double [block_number];
double *prog_t = new double [block_number];
double *prog_pot_2 = new double [block_number];
double *prog_kin_2 = new double [block_number];
double *prog_tot_2 = new double [block_number];
double *prog_t_2 = new double [block_number];
double *s_pot = new double [block_number];
double *s_kin = new double [block_number];
double *s_tot = new double [block_number];
double *s_t = new double [block_number];

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void DataBlock();
void Argon();
