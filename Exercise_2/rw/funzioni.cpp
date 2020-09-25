#include "funzioni.hh"
#include "random.h"
#include <fstream>

using namespace std;

double norm(double *p) {
   return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
}

void move_pos(double *pos, double *move) {
  for (unsigned int i = 0; i < 3; i++)
    pos[i] += move[i];
}
