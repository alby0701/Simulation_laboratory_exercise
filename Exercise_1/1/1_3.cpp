
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
   unsigned int n_throws = 1000000;
   unsigned int n_blocks = 100;
   double d = 2.;
   double L = 1.;
if(argc > 1)
{
  	n_throws = atoi(argv[1]);
    n_blocks = atoi(argv[2]);
}
unsigned int throws_per_block = n_throws/n_blocks;
Random rnd;
int seed[4] {0,0,0,1};
int p1 = 2892;
int p2 = 2587;
rnd.SetRandom(seed, p1, p2);

double* mean = new double[n_blocks];
double* err = new double[n_blocks];

double n_hit,
       mean_cum(0),
       mean_sq_cum(0);
for(unsigned int j = 0; j < n_blocks; j++){
  n_hit = 0;
  for(unsigned int i = 0; i < throws_per_block; i++){
    double x = rnd.Rannyu()*d;
    double theta = rnd.Rannyu()*M_PI/2.;
    double x_1 = x + L*cos(theta);
    if(x_1 > d)
      n_hit += 1;
  }
  if (n_hit == 0) {
    // discard this block, I would need to divide by 0
    j--;
    continue;
  }
  n_hit = 2*L*throws_per_block/(d*n_hit);//n_hit is now our pi approx

  mean_cum    += n_hit;
  mean_sq_cum += pow(n_hit, 2);
  mean[j]      = mean_cum / (j + 1);
  err[j]       = sqrt((mean_sq_cum / (j + 1) - pow(mean[j],2)) / j);
}
err[0] = 0; // handle that bad NaN, the fast 'n' easy way :)

ofstream out;
out.open("buffon.dat");
for(unsigned int j = 0; j < n_blocks; j++){
 out << mean[j] << " " << err[j] << endl;
}
delete []mean;
delete []err;
return 0;
}
