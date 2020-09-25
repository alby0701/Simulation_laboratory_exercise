#include <cmath>
#include <vector>
#include <valarray>
#include <iostream>
#include <fstream>
#include "tsp_sa.h"

int main(int argc, char* argv[]){
    if (argc < 3) { // We expect 3 arguments: the program name, the data_file path and the config_file path
        std::cerr << "Usage: " << argv[0] << " data_file config_file" << std::endl;
        return 1;
    }
    std::ifstream data_file(argv[1]);
    /* file storing the coordinates of each city
        x0 y0 // the city we start from
        x1 y1
        ...
    */
    std::ifstream config_file(argv[2]);
    /* file storing the hyperparameters value
        POP_SIZE
        NUM_GENERATIONS
        MUTATION_PROB
        CROSSOVER_PROB
        POWER_LAW_EXPONENT
    */

    std::vector<std::valarray<double>> X;
    double x, y;
    while (data_file >> x >> y) {
        X.push_back(std::valarray<double> {x, y});
    }

    /* hyperparameters */
    unsigned int NUM_ITERATIONS;
    unsigned int MONTE_CARLO_STEPS;
    double MAX_TEMP;

    config_file >> NUM_ITERATIONS;
    config_file >> MONTE_CARLO_STEPS;
    config_file >> MAX_TEMP;

    TSPsa t;
    t.set_hyperparameters(NUM_ITERATIONS, MONTE_CARLO_STEPS, MAX_TEMP);
    t.set_data(X);
    t.run();

    std::ofstream cost_output;
    /* file that will store the cost of the best and the mean of the first half of he population, for each generation
       L_best[0] L_mean[0]
       L_best[1] L_mean[1]
       ...
    */
    cost_output.open("cost.out", std::ios::out);
    for (unsigned int i = 0; i < NUM_ITERATIONS; i++) {
        cost_output << t.L[i] << std::endl;
        // std::cout << t.L[i] << std::endl;
    }
    cost_output.close();

    std::vector<unsigned int> solution (t.get_solution()); // the indexes of the cities, in the best order found
    std::cout << "0 ";
    // for (unsigned int i=0; i<t.NUM_CITIES-1; ++i) std::cout << t.x[i] << " ";
    for (auto p: solution) std::cout << p << " ";
    std::cout << "0" << std::endl;

    return 0;
}
