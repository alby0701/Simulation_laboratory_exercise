#include <vector>
#include <valarray>
// #include <algorithm>
#include <random>
#include <functional>
#include <iostream>

#ifndef __TSPsa__
#define __TSPsa__

class TSPsa {
    public:
        TSPsa(){};
        /* methods */

        void build_cost_matrix(std::vector<std::valarray<double>>);
        double get_unit_cost(unsigned int, unsigned int);
        double total_cost(std::vector<unsigned int>);

        void generate_initial_config();
        void set_data(std::vector<std::valarray<double>>);
        void set_hyperparameters(unsigned int,
                                 unsigned int,
                                 double);

        /* genetic operators */
        std::vector<unsigned int> new_mc_config();
        void simulated_annealing(double);
        double transition_prob(double, double, double);
        void mc(double);
        void run();
        std::vector<unsigned int> get_solution();

        std::vector<double> L;

        double rand();
        int randint(int, int);
        int accepted, attempted;
        std::vector<unsigned int> x; // config
        double L_x; // cost of this config

        unsigned int NUM_CITIES;

    private:
        /* hyperparameters */
        unsigned int NUM_ITERATIONS;
        unsigned int MONTE_CARLO_STEPS;
        double MAX_TEMP;
        double DELTA_T;

        std::valarray<std::valarray<double>> cost_matrix;

        /* methods */
        double norm_s(std::valarray<double> x1, std::valarray<double> x2);

        /*PRNG*/
        // std::default_random_engine generator;
        std::random_device generator;
        std::uniform_real_distribution<double> distrib;
};


#endif // __Random__
