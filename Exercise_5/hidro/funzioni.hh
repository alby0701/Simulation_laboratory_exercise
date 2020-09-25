#include "random.h"

using namespace std;

const double a0 = 1;//0.0529;//nm!!!!
const double exp_r100 =(3.0/2.0)*a0;
const double exp_r210 =5.0*a0;
double psi_100(double *p);
double psi_210(double *p);
double norm(double *p);
void move_pos(double *pos, double *move);
bool Metropolis(double *p, double m_d, double theta, double phi, Random rnd_s, bool type);
