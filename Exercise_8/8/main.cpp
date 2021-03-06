#include <cmath>
#include <valarray>
#include <functional>
#include <iostream>
#include <fstream>
#include "random.h"
   
double psiT(double x, double mu, double sigma2){
    return std::exp(-std::pow(x-mu,2)/(double)2/sigma2) + std::exp(-std::pow(x+mu, 2)/(double)2/sigma2);
}

double V(double x){
    // potential
    return std::pow(x,4) - 5/(double)2*std::pow(x,2);
}

double H(double x, double mu, double sigma2){
    double exp_minus = std::exp(-std::pow(x-mu, 2)/(2. * sigma2));
    double exp_plus = std::exp(-std::pow(x+mu, 2)/(2. * sigma2));
    return -(-exp_minus/sigma2 - exp_plus/sigma2 + (std::pow(mu - x, 2)* exp_minus)/std::pow(sigma2,2) + (std::pow(mu + x, 2) * exp_plus)/std::pow(sigma2,2) )/(2.*psiT(x, mu, sigma2)) + V(x);
    //return -(-std::exp(-std::pow(x-mu, 2)/(2 * sigma2))/sigma2 - std::exp(-std::pow(mu + x, 2)/(2 * sigma2))/sigma2 + (std::pow(mu - x, 2)* std::exp(-std::pow(mu - x, 2)/(2 * sigma2)))/std::pow(sigma2,2) + (std::pow(mu + x, 2) * std::exp(-std::pow(mu + x, 2)/(2 * sigma2)))/std::pow(sigma2,2) )/psiT(x, mu, sigma2) + V(x);
}


double T_uniform(double x, double step, Random *rand){
    return x + rand->Rannyu(-step, step);
}
struct block_value{
    double properties_block;
    double properties_cum_average;
    double properties_cum_average2;
    double properties;
    double properties_err;
};
   

double metropolis(std::function<double(double, double step, Random *rand)> T, 
                double x,
                double step,
                unsigned int n_samples,
                Random *rand,
                double mu,
                double sigma2,
                block_value *block_array,
                bool verbose_pdf) {
    double Hmean = 0;
    std::ofstream x_hist_pdf("x_hist_pdf.dat", std::ios::app);

                   
    // first step
//     std::cout << p[0] << " " << p[1] <<  " " << p[2] << "\t\t" << norm(p) << " " << pdf(p) << std::endl;
    double psi2_curr = std::pow(psiT(x, mu, sigma2), 2); // squared norm of Psi_T(x)
    for (unsigned int n = 0; n < n_samples - 1; n++) {
        double x_new = T(x, step, rand);
        double psi2_new = std::pow(psiT(x_new, mu, sigma2), 2);
        
        if (rand->Rannyu() <= std::min(1., psi2_new / psi2_curr)) { // accept
            x = x_new;
            psi2_curr = psi2_new;
        }
        if (verbose_pdf) x_hist_pdf << x << std::endl;
//         std::cout << p[0] << " " << p[1] <<  " " << p[2] << "\t\t" << norm(p) << " " << pdf(p) << std::endl;
//         std::cout << norm(p) << std::endl;
        Hmean += H(x, mu, sigma2);
        block_array->properties_block += H(x, mu, sigma2);
    }
    // std::cout << Hmean/(double)n_samples << std::endl;
    return Hmean/(double)n_samples;
}

void MeasureBlock(int block, int steps, block_value *block_array) {
    block_array->properties_block /= (double)steps; 
    block_array->properties_cum_average += block_array->properties_block;
    block_array->properties_cum_average2 += std::pow(block_array->properties_block, 2);
    block_array->properties = block_array->properties_cum_average / (double)(block + 1);
    if (block != 0)
        block_array->properties_err = std::sqrt((block_array->properties_cum_average2 / (double)(block + 1) - std::pow(block_array->properties, 2)) / (double)block);

    // reset block variables
    block_array->properties_block = 0;
    // print block
    std::ofstream output("output.dat", std::ios::app);
    output << block_array->properties    << " " << block_array->properties_err << std::endl;
}

 
int main(int argc, char* argv[]){
    Random rand;
    int iv[] = {0,0,0,1};
    rand.SetRandom(iv, 2892, 2587); 
    if (argc < 6) {
        // mu, sigma, nsteps, step
        std::cout << "Usage: " << argv[0] << " mu sigma nsteps step nblocks [verbose_pdf]" << std::endl;
        exit(1);
    }
    bool verbose_pdf = false;
    if (argc > 6) {
        verbose_pdf = true;
        std::ofstream x_hist_pdf("x_hist_pdf.dat", std::ios::out); // clean file
        x_hist_pdf.close();
    }

    double mu = std::atof(argv[1]);
    double sigma2 = std::pow(std::atof(argv[2]), 2);
    int nsteps = std::atoi(argv[3]);
    int step = std::atoi(argv[4]);
    unsigned int n_blocks = std::atoi(argv[5]);

    struct block_value hamiltonian;
    double Hmean;
    // double Hmean_err = 0;

    std::ofstream output("output.dat", std::ios::out);
    output.close();

    for (unsigned int n = 0; n < n_blocks; n++)    {
        Hmean = metropolis(T_uniform, 0, step, nsteps, &rand, mu, sigma2, &hamiltonian, verbose_pdf);
        MeasureBlock(n, nsteps,&hamiltonian); // compute mean and error for the block, reset block-specific var's
    }
    std::cout << Hmean << std::endl;
    
}
