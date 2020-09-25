#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.hh"

using namespace std;



int main (int argc, char *argv[]){
   if (argc != 4) {
     cerr << "Missing parameters!" << endl;
     cerr << "Usage: " << argv[0] << " <throws> <steps> <bool:lattice>" << endl;
     exit(1);
   }
   unsigned int M = atoi(argv[1]);
   unsigned int steps = atoi(argv[2]);
   bool lattice = atoi(argv[3]);
   Random rnd;
   int seed[4] {0,0,0,1};
   int p1 = 2892;
   int p2 = 2587;
   rnd.SetRandom(seed, p1, p2);

   double *mean = new double[steps + 1];
   double *err = new double[steps + 1];
   for (unsigned int i = 0; i < steps + 1; i++) {
     mean[i] = 0;
     err[i] = 0;
   }


   double possible_moves[6][3] = { {1,0,0}, {-1,0,0},
                                   {0,1,0}, {0,-1,0},
                                   {0,0,1}, {0,0,-1} };

   for (unsigned int j = 0; j < M; j++) {
     double pos[3] {0,0,0};
     for (unsigned int i = 1; i < steps + 1; i++) {
       if (lattice == 1) {
         move_pos(pos, possible_moves[rnd.Rannyu_i(0, 6)] );
       } else {
         double theta = rnd.Rannyu(0, M_PI);
         double phi = rnd.Rannyu(0, 2*M_PI);
         double move[3] {cos(phi)*sin(theta),cos(phi)*cos(theta),sin(phi)};
         //cout << move[0] << " " << move[1] << " " << move[2] << endl;
         move_pos(pos, move);
       }
       double r = norm(pos);
       mean[i] += sqrt(r);
       err[i] += r;
     }
   }

   ofstream out;
   out.open("reti.dat");
   for (unsigned int i = 0; i < steps + 1; i++) {
     mean[i] /= M;
     err[i] /= M;
     err[i] = sqrt((err[i] - mean[i]*mean[i]) / (M - 1));
     out << mean[i] << " " << err[i] << endl;
   }
   out.close();

   delete []mean;
   delete []err;
   return 0;
}
