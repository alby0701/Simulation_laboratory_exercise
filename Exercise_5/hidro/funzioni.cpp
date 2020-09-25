#include "funzioni.hh"
#include "random.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double norm(double *p) {
   return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

void move_pos(double *pos, double *move) {
  for (unsigned int i = 0; i < 3; i++)
    pos[i] += move[i];
}

double psi_100(double *p)
{
  return 1./M_PI*exp(-2.*(norm(p)));
}
double psi_210(double *p)
{
  return 1./(32.*M_PI)*exp(-(norm(p)))*pow(p[2],2);//cos(acos(p[3]/((norm(p)))));
}

bool Metropolis(double *p, double m_d, double theta, double phi, Random rnd_s, bool type)
{
  double p_n[3] {p[0],p[1],p[2]};
  //cout << "r = " << m_d << " theta = " << theta << " phi= " << phi <<  endl;
  double move[3] {cos(phi)*sin(theta)*m_d,sin(phi)*sin(theta)*m_d,cos(theta)*m_d};//{cos(phi)*sin(theta)*m_d,cos(phi)*cos(theta)*m_d,sin(phi)*m_d};
  move_pos(p_n, move);
  double alpha = 0;
  //cout << p[0] << " " << p[1] << " " << p[2] << endl;
  // cout << p_n[0] << " " << p_n[1] << " " << p_n[2] << endl;
  // cout << psi_100(p_n) << endl;
  if(type == false)
  {
    alpha = min(1.,psi_100(p_n)/psi_100(p));//min(1.,(abs(psi_100(p_n)*psi_100(p_n)))/(abs(psi_100(p)*psi_100(p))));
  }
  else
  {
    alpha = min(1.,psi_210(p_n)/psi_210(p));//min(1.,(abs(psi_210(p_n)*psi_210(p_n)))/(abs(psi_210(p)*psi_210(p))));
  }
  double r = rnd_s.Rannyu(0, 1);
  //cout << "alpha = " << alpha << " r = " << r <<  endl;
  if(r <= alpha)
  {
    //cout << "accept" << endl;
    p[0] = p_n[0];
    p[1] = p_n[1];
    p[2] = p_n[2];
    return true;
  }
  else
    return false;
}
void Metropolis_h_gaussian(double *p, double *p_old, double *trans, Random rnd_s, bool type)
{}
