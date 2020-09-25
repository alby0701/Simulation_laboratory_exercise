#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;



int main (int argc, char *argv[]){
   if (argc != 3) {
     cerr << "Missing parameters!" << endl;
     cerr << "Usage: " << argv[0] << "<n_throws> <n_blocks>" << endl;
     exit(1);
    }
     unsigned int n_throws = atoi(argv[1]);
     unsigned int n_blocks = atoi(argv[2]);
     Random rand;
     int iv[] = {0,0,0,1};
     rand.SetRandom(iv, 2892, 2587);
     double S_0 = 100;
     double T = 1;
     double K = 100;
     double r = 0.1;
     double sigma = 0.25;
     unsigned int t_interval = 100;

     double S;
     double call_cum = 0;
     double put_cum = 0;
     double call2_cum = 0;
     double put2_cum = 0;
     double dt = T/t_interval;

     double *call = new double[n_blocks];
     double *err_call = new double[n_blocks];

     double *put = new double[n_blocks];
     double *err_put = new double[n_blocks];

     unsigned int throws_per_block = n_throws/n_blocks;
     for(unsigned int j = 0; j < n_blocks; j++){
       double call_block(0), put_block(0);
       for(unsigned int i = 0; i< throws_per_block; i++){
         S = S_0;
         for(unsigned int k = 0; k< t_interval; k++)
          S *= exp((r - pow(sigma,2)/2.)*dt + sigma*rand.Gauss(0,1)*sqrt(dt));
       call_block += max(S-K,0.)*exp(-r*T);
       put_block +=  max(K-S,0.)*exp(-r*T);
      }
      call_block /= throws_per_block;
      put_block /= throws_per_block;
      call_cum += call_block;
      call2_cum += pow(call_block,2);
      call[j] = call_cum/(double(j)+1);
      err_call[j] = sqrt((call2_cum/(double(j)+1) - (call[j]*call[j]))/(double(j)));
      put_cum += put_block;
      put2_cum += pow(put_block,2);
      put[j] = put_cum/(double(j)+1);
      err_put[j] = sqrt((put2_cum/(double(j)+1) - put[j]*put[j])/(double(j)));
     }
     err_put[0] = 0;
     err_call[0] = 0;
   ofstream out;
   out.open("gbm_d.dat");
   for (unsigned int i = 0; i < n_blocks; i++) {
     out << call[i] << " " << err_call[i] << " " << put[i] << " " << err_put[i] << endl;
   }
   out.close();

   delete []call;
   delete []put;
   delete []err_put;
   delete []err_call;
   return 0;
}
