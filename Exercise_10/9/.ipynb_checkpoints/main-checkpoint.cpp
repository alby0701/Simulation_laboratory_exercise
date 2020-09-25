#include <cmath>
#include <vector>
#include <valarray>
#include <random>
// #include "tsp.h"

#include <iostream>


    
int main(){
    std::default_random_engine generator;
    std::uniform_real_distribution<double> randunif(0.0,1.0);
    std::cout << randunif(generator) << std::endl;
    
    unsigned int NUM_CITIES = 32;
    std::valarray<std::valarray<double>> X (NUM_CITIES); // (x,y) position of each city
    for (unsigned int i = 0; i < NUM_CITIES; i++) {
        X[i].resize(2);
        X[i][0] = randunif(generator);
        X[i][1] = randunif(generator);
    }

    /* hyperparameters */
    unsigned int POP_SIZE = 1000;
    unsigned int NUM_GENERATIONS = 100;
    double MUTATION_PROB = 0.1;
    double CROSSOVER_PROB = 0.5;
    double POWER_LAW_EXPONENT = 2;
}
