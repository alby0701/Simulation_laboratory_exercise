#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.hh"

using namespace std;



int main (int argc, char *argv[]){
   if (argc != 3) {
     cerr << "Missing parameters!" << endl;
     cerr << "Usage: " << argv[0] << " <throws> <steps>" << endl;
     exit(1);
   }
   unsigned int M = atoi(argv[1]);
   unsigned int steps = atoi(argv[2]);
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


   double possible_moves[6][3] = { {a0,0,0}, {-a0,0,0},
                                   {0,a0,0}, {0,-a0,0},
                                   {0,0,a0}, {0,0,-a0} };

   for (unsigned int j = 0; j < M; j++) {
     double pos[3] {a0*10,a0*10,a0*10};
     double pos_old[3] {a0*10,a0*10,a0*10};
     for (unsigned int i = 1; i < steps + 1; i++) {
       do{
         move_pos(pos, possible_moves[rnd.Rannyu_i(0, 6)] );
         Metropolis_h_uniform(pos, pos_old, rnd);
         cout << pos[0] << " " << pos[1] << " " << pos[2] << endl;
       }while((pos_old[0] == pos[0]) & (pos_old[1] == pos[1]) & (pos_old[2] == pos[2]));
       pos_old[0] = pos[0];
       pos_old[1] = pos[1];
       pos_old[2] = pos[2];

       //cout << pos[0] << " " << pos[1] << " " << pos[2] << endl;
       //double r = norm(pos);
       //mean[i] += sqrt(r);
       //err[i] += r;
     }
   }

   /*ofstream out;
   out.open("reti.dat");
   for (unsigned int i = 0; i < steps + 1; i++) {
     mean[i] /= M;
     err[i] /= M;
     err[i] = sqrt((err[i] - mean[i]*mean[i]) / (M - 1));
     out << mean[i] << " " << err[i] << endl;
   }
   out.close();*/

   delete []mean;
   delete []err;
   return 0;
}
