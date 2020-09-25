#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "funzioni.hh"

using namespace std;



int main (int argc, char *argv[]){
   if (argc != 4) {
     cerr << "Missing parameters!" << endl;
     cerr << "Usage: " << argv[0] << " <n_blocks> <n_steps> <distribution>" << endl;
     cerr << "Starting point file <origin.dat>:" << endl;
     cerr << "<pos100[i]> <pos210[i]> (cartesian coordinates)" << endl;
     cerr << "<n_safe_steps>" << endl;
     exit(1);
   }
   cout << endl;
   cout << "New starting position, and safe-steps for -too far away start- test in <origin.dat>" << endl;
   cout << endl;
   unsigned int M = atoi(argv[1]);
   unsigned int steps = atoi(argv[2]);
   unsigned int distribution = atoi(argv[3]);
   //cout << distribution << endl;
   unsigned int step_block = steps/M;
   double pos[3] {0,0,0};
   double pos2[3] {0,0,0};
   unsigned int safe_steps = 0;
   ifstream in;
   in.open("origin.dat");
   in >> pos[0] >> pos[1] >> pos[2];
   in >> pos2[0] >> pos2[1] >> pos2[2];
   in >> safe_steps;
   in.close();
   for (unsigned int i = 0; i < 3; i++)
   {
     pos[i] = pos[i]*a0;
     pos2[i] = pos2[i]*a0;
   }
   Random rnd;
   int seed[4] {0,0,0,1};
   int p1 = 2892;
   int p2 = 2587;
   rnd.SetRandom(seed, p1, p2);
   double move_dimension100 = 2.6*a0;
   double move_dimension210 = 6.7*a0;
   double *x100 = new double[(steps + 1)];
   double *y100 = new double[(steps + 1)];
   double *z100 = new double[(steps + 1)];
   double *mean100 = new double[M];
   double *err100 = new double[M];

   double *x210 = new double[(steps + 1)];
   double *y210 = new double[(steps + 1)];
   double *z210 = new double[(steps + 1)];
   double *mean210 = new double[M];
   double *err210 = new double[M];

   ///EQUILIBRATION PROCEDURE/////////
   for (unsigned int i = 0; i < M; i++) {
     mean100[i] = 0;
     err100[i] = 0;
     mean210[i] = 0;
     err210[i] = 0;
   }
  cout << "starting simulation from: " << endl;
  cout << "x100 = " << pos[0]/a0 << "*a0" << ", y100 = " << pos[1]/a0 << "*a0" << ", z100 = " << pos[2]/a0 << "*a0" << endl;
  cout << "x210 = " << pos2[0]/a0 << "*a0" << ", y210 = " << pos2[1]/a0 << "*a0" << ", z210 = " << pos2[2]/a0 << "*a0" << endl;
  cout << endl;
  cout << "Step dimension Y100 for A = 50%:" << endl;
  cout <<  move_dimension100/a0 << "*a0;" << endl;
  cout << "Step dimension Y210 for A = 50%:" << endl;
  cout <<  move_dimension210/a0 << "*a0;" << endl;
  cout << endl;
  cout << "Equilibrate procedure, safe step chosen:" << endl;
  cout <<  safe_steps << " Steps." << endl;
  cout << endl;
  double lim_100[3] {(4/3)*a0,(4/3)*a0,(4/3)*a0};
  double lim_210[3] {10*a0,10*a0,10*a0};
  double norm_lim100 = norm(lim_100);
  double norm_lim210 = norm(lim_210);

    ////VERIFICATION POSITION//////// (TO BE VERIFIED AND CORRECT)
  // while((norm(pos) > norm_lim100)){
  //   double pos_old[3] {pos[0],pos[1],pos[2]};
  //   for (unsigned int i = 1; i < safe_steps; i++) {
  //       double theta = rnd.Rannyu(0, M_PI);
  //       double phi = rnd.Rannyu(0, 2*M_PI);
  //       Metropolis(pos_old, move_dimension100,theta, phi, rnd, 0);
  //       pos[0] = pos_old[0];
  //       pos[1] = pos_old[1];
  //       pos[2] = pos_old[2];
  //       }
  //       if(norm(pos) > 2*norm_lim100)
  //       {
  //         cout << "The starting point is too far away, it's impossible to find convergence with actual step for 100 configuration." << endl;
  //         exit(2);
  //       }
  //     }
  // while((norm(pos2) > norm_lim210)){
  //   double pos_old[3] {pos2[0],pos2[1],pos2[2]};
  //   for (unsigned int i = 1; i < safe_steps; i++)
  //   {
  //     double theta = rnd.Rannyu(0, M_PI);
  //     double phi = rnd.Rannyu(0, 2*M_PI);
  //     bool accept = false;
  //     Metropolis(pos2, move_dimension210,theta, phi, rnd, 1);
  //     pos2[0] = pos_old[0];
  //     pos2[1] = pos_old[1];
  //     pos2[2] = pos_old[2];
  //     }
  //       //cout << "x210 = " << pos2[0]/a0 << "*a0" << ", y210 = " << pos2[1]/a0 << "*a0" << ", z210 = " << pos2[2]/a0 << "*a0" << endl;
  //   if(norm(pos2) > 2*norm_lim210)  {
  //     cout << "The starting point is too far away, it's impossible to find convergence with actual step for 210 configuration." << endl;
  //     exit(3);
  //     }
  //   }
  //   if((norm(pos2) < norm_lim210) || (norm(pos) < norm_lim100))
  //   {
  //     cout << "Equilibration procedure not necessary for this run." << endl;
  //     cout << endl;
  //   }//*/
    //////SAMPLING FASE////////
    cout << "Start sampling from configuration:" << endl;
    cout << "100 : "<< "x = " << pos[0]/a0 << "*a0" << ", y = " << pos[1]/a0 << "*a0" << ", z = " << pos[2]/a0 << "*a0" << endl;
    cout << "210 : "<< "x = " << pos2[0]/a0 << "*a0" << ", y = " << pos2[1]/a0 << "*a0" << ", z = " << pos2[2]/a0 << "*a0" << endl;
    x100[0] = pos[0];
    y100[0] = pos[1];
    z100[0] = pos[2];
    x210[0] = pos2[0];
    y210[0] = pos2[1];
    z210[0] = pos2[2];
    double  block_mean100,
            err_mean100,
            mean_cum100(0),
            mean_sq_cum100(0);
    double  block_mean210,
            err_mean210,
            mean_cum210(0),
            mean_sq_cum210(0);
            int reject = 0;
    for (unsigned int j = 0; j < M; j++) {
      block_mean100 = 0;
      err_mean100 = 0;
      block_mean210 = 0;
      err_mean210 = 0;

      for (unsigned int i = 1; i < step_block+1; i++) {
        double pos_old[3] {x100[i+step_block*j-1],y100[i+step_block*j-1],z100[i+step_block*j-1]};
        double phi, theta;
        if (distribution == 0)  {
          move_dimension100 = rnd.Rannyu(0, 0.6);
          theta = rnd.Rannyu(0, 2*M_PI);
          phi = rnd.Rannyu(0, M_PI);
        }
        else  {
          move_dimension100 = rnd.Gauss(0.0, 0.6);
          theta = rnd.Gauss(0, 2*M_PI);
          phi = rnd.Gauss(0, M_PI);
        }
        Metropolis(pos_old, move_dimension100,theta, phi, rnd, 0);
        x100[i+step_block*j] = pos_old[0];
        y100[i+step_block*j] = pos_old[1];
        z100[i+step_block*j] = pos_old[2];
        double r100 = norm(pos_old);
        block_mean100 += r100;
       //err_mean100 += pow(r100,2);
      }
       for (unsigned int i = 1; i < step_block+1; i++) {
         double pos_old[3] {x210[i+step_block*j-1],y210[i+step_block*j-1],z210[i+step_block*j-1]};
         double phi, theta;
         if (distribution == 0)  {
           move_dimension210 = rnd.Rannyu(0, 6.7);
           theta = rnd.Rannyu(0, 2*M_PI);
           phi = rnd.Rannyu(0, M_PI);
         }
         else  {
           move_dimension210 = rnd.Gauss(0, 6.7);
           theta = rnd.Gauss(0, 2*M_PI);
           phi = rnd.Gauss(0, M_PI);
         }
         Metropolis(pos_old, move_dimension210,theta, phi, rnd, 1);
       x210[i+step_block*j] = pos_old[0];
       y210[i+step_block*j] = pos_old[1];
       z210[i+step_block*j] = pos_old[2];
       //cout << pos[0] << " " << pos[1] << " " << pos[2] << endl;
       double r210 = norm(pos_old);
       block_mean210 += r210;
       //err_mean210 += pow(r210,2);
     }
     block_mean100 /= step_block;
     block_mean210 /= step_block;

     mean_cum100    += block_mean100;
     mean_sq_cum100 += pow(block_mean100, 2);
     mean100[j]      = mean_cum100 / (j + 1);
     err100[j]       = sqrt((mean_sq_cum100 / (j + 1) - pow(mean100[j],2)) / j);

     mean_cum210    += block_mean210;
     mean_sq_cum210 += pow(block_mean210, 2);
     mean210[j]      = mean_cum210 / (j + 1);
     err210[j]       = sqrt((mean_sq_cum210 / (j + 1) - pow(mean210[j],2)) / j);
   }
   err100[0] = 0;
   err210[0] = 0;
   ofstream out;
   out.open("100rm.dat");
   for (unsigned int j = 1; j < M; j++) {
    out << mean100[j] << " " << err100[j]  << endl;
  }
  out.close();
  out.open("210rm.dat");
  for (unsigned int j = 1; j < M; j++) {
      out << mean210[j] << " " << err210[j]  << endl;
  }
   out.close();
   out.open("100.dat");
   for (unsigned int i = 1; i < (steps + 1); i++) {
     out << x100[i] << " " << y100[i] << " " << z100[i] << " " << sqrt(x100[i]*x100[i] + y100[i]*y100[i] + z100[i]*z100[i]) << endl;
   }
   out.close();
   out.open("210.dat");
   for (unsigned int i = 1; i < (steps + 1); i++) {
     out << x210[i] << " " << y210[i] << " " << z210[i] << " " << sqrt(x210[i]*x210[i] + y210[i]*y210[i] + z210[i]*z210[i]) << endl;
   }
   out.close();

   delete []x100;
   delete []y100;
   delete []z100;
   delete []mean100;
   delete []err100;
   delete []x210;
   delete []y210;
   delete []z210;
   delete []mean210;
   delete []err210; //*/
   return 0;
}
