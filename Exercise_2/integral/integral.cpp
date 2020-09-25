#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
    if (argc != 4) {
      cerr << "Missing parameters!" << endl;
      cerr << "Usage: " << argv[0] << " <throws> <blocks> <bool:i_sampling>" << endl;
      return 1;
    }
    unsigned int n_throws = atoi(argv[1]);
    unsigned int n_blocks = atoi(argv[2]);
    bool i_sampling = atoi(argv[3]);

    Random rnd;
    int seed[4] {0,0,0,1};
    int p1 = 2892;
    int p2 = 2587;
    rnd.SetRandom(seed, p1, p2);

    function<double(double)> p, cdf_inv;

    // user defined function
    auto f = [](double x) { return (M_PI/2.)*cos((M_PI/2.)*x); };
    if (i_sampling) {
      // Importance sampling case
      p = [](double x) { return 2*(1-x); };
      cdf_inv = [](double x) { return 1+sqrt(x); };
    } else {
      // Uniform case
      p = [](double x) { return 1; };
      cdf_inv = [](double x) { return x; };
    }
    double a_int = 0;
    double b_int = 1;
    /////////////////////

    auto g = [&f,&p](double x) { return f(x)/p(x); };

    double* mean = new double[n_blocks];
    double* err = new double[n_blocks];

    double block_mean,
           mean_cum(0),
           mean_sq_cum(0);
    for(unsigned int j = 0; j < n_blocks; j++){
      block_mean = 0;
      for(unsigned int i = 0; i < n_throws/n_blocks; i++){
        double u = cdf_inv(rnd.Rannyu(a_int, b_int));
        block_mean += g(u);
      }
      block_mean /= n_throws/n_blocks;

      mean_cum    += block_mean;
      mean_sq_cum += pow(block_mean, 2);
      mean[j]      = mean_cum / (j + 1);
      err[j]       = sqrt((mean_sq_cum / (j + 1) - pow(mean[j],2)) / j);
  }
  err[0] = 0; // handle that bad NaN, the fast 'n' easy way :)

  ofstream out;
  out.open("integral.dat");
  for(unsigned int j = 0; j < n_blocks; j++){
     out << mean[j] << " " << err[j] << endl;
  }
  delete []mean;
  delete []err;
  return 0;
}
