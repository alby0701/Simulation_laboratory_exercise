#include <vector>
#include <valarray>
#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>

#include "tsp_sa.h"

double TSPsa::rand(){
    return distrib(generator);
}

int TSPsa::randint(int min, int max){
    return std::floor(min+(max-min)*rand());
}

void TSPsa::set_data(std::vector<std::valarray<double>> X) {
    NUM_CITIES = X.size();
    build_cost_matrix(X);
}

void TSPsa::set_hyperparameters(unsigned int num_iterations,
                                unsigned int mc_steps,
                                double max_temp) {
    NUM_ITERATIONS = num_iterations;
    MONTE_CARLO_STEPS = mc_steps;
    MAX_TEMP = max_temp;
    DELTA_T = (double)MAX_TEMP/NUM_ITERATIONS;

    L.resize(NUM_ITERATIONS);
}

double TSPsa::norm_s(std::valarray<double> x1, std::valarray<double> x2) {
    return std::pow(x2-x1, 2).sum();
}

std::vector<unsigned int> TSPsa::new_mc_config() {
    std::vector<unsigned int> sequence(x);
    unsigned int pos = randint(1,NUM_CITIES-2); // which city is affected
    std::swap(sequence[pos], sequence[pos+1]);

    return sequence;
}

void TSPsa::build_cost_matrix(std::vector<std::valarray<double>> X) {
    // store as a lower triangular
    cost_matrix.resize(NUM_CITIES-1);
    for (unsigned int i = 1; i < NUM_CITIES; i++) {
        cost_matrix[i-1].resize(i);
        for (unsigned int j = 0; j < i; j++)
            cost_matrix[i-1][j] = norm_s(X[i], X[j]);
    }
}

double TSPsa::get_unit_cost(unsigned int i, unsigned int j) {
    /* Return the entry (i,j) of the symmetric cost matrix */
    if (i > j) return cost_matrix[i-1][j];
    else return cost_matrix[j-1][i];
}

double TSPsa::total_cost(std::vector<unsigned int> sequence){
    double L = get_unit_cost(0, sequence.front()); /* first city is 0 */
    for (unsigned int i = 0; i < NUM_CITIES - 2; i++)
        L += get_unit_cost(sequence[i], sequence[i+1]);
    L += get_unit_cost(0, sequence.back());
    return L;
}

void TSPsa::generate_initial_config() {
    // random permutation of the cities
    x.resize(NUM_CITIES-1);
    for (unsigned int i=1; i<NUM_CITIES; ++i) x[i-1] = i;
    std::random_shuffle(x.begin(), x.end());
    L_x = total_cost(x);
}

double TSPsa::transition_prob(double L_old, double L_new, double temp) {
    return (L_new > L_old) ? std::exp((L_old - L_new)/temp) : 1;
}

/****** Simulated annealing *******/
void TSPsa::simulated_annealing(double temp) {
    std::vector<unsigned int> x_new (new_mc_config());
    double L_new = total_cost(x_new);

    // Metropolis sampling from transition_prob distrib
    double p = transition_prob(L_x, L_new, temp);
    if (p > rand()) {
        // we accept the move
        // std::cout << "p " << p << std::endl;
        // x = std::move(x_new);
        x = x_new;
        L_x = L_new;
        accepted++;
        // std::cout << "L " << L_new << std::endl;
    } else {
        // std::cout << "rejected" << std::endl;
    }
    attempted++;
}

/* Cycle over N-Monte Carlo Steps */
void TSPsa::mc(double temp) {
    for (unsigned int i = 0; i < MONTE_CARLO_STEPS; i++)
        simulated_annealing(temp);
}

void TSPsa::run() {
    generate_initial_config();
    accepted = 0;
    attempted = 0;
    double temp = MAX_TEMP;
    for (unsigned int i = 0; i < NUM_ITERATIONS; i++) {
          //std::cout << "***********************************" << std::endl;
          mc(temp);
        //   std::cout << "T " << temp << ", " << L_x << "\t a_rate " << (double)accepted/attempted << std::endl;
          L[i] = L_x;
          temp -= DELTA_T;
    }
}

std::vector<unsigned int> TSPsa::get_solution() {
    return x;
}
