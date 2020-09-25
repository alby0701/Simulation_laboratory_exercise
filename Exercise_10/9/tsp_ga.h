#include <vector>
#include <valarray>
// #include <algorithm>
#include <random>
#include <functional>
#include <iostream>

#ifndef __TSPga__
#define __TSPga__

// class TSPga;
struct Individual {
  std::vector<unsigned int> sequence;
  double cost;
  bool operator>(const Individual &p) const { return cost > p.cost; };
  bool operator<(const Individual &p) const { return cost < p.cost; };

  Individual(){};
  // Individual(std::vector<unsigned int> i_sequence)
  // {
  //   sequence = i_sequence;
  //   cost = TSPga::total_cost(sequence);
  // };
  Individual(std::vector<unsigned int> i_sequence, double t_cost)
  {
    sequence = i_sequence;
    cost = t_cost;
    // std::cout << "cost " << cost << std::endl;
    // std::cout << "seq " << sequence[1] << std::endl;
  };
};

class TSPga {
    public:
        TSPga(){
          distrib = std::uniform_real_distribution<double>(0.,1.);

          /* register mutations */
          mutations.push_back(std::bind(&TSPga::mutation_perm, this, std::placeholders::_1));
          mutations.push_back(std::bind(&TSPga::mutation_shift, this, std::placeholders::_1));
          mutations.push_back(std::bind(&TSPga::mutation_swap, this, std::placeholders::_1));
          mutations.push_back(std::bind(&TSPga::mutation_inversion, this, std::placeholders::_1));
          // mutations.push_back(mutation_perm);
          // mutations.push_back(mutation_shift);
          // mutations.push_back(mutation_swap);
          // mutations.push_back(mutation_inversion);

          };
        /* methods */
        void build_cost_matrix(std::vector<std::valarray<double>>);
        double get_unit_cost(unsigned int, unsigned int);
        double total_cost(std::vector<unsigned int>);

        // std::vector<Individual> generate_initial_population();
        void generate_initial_population();
        void set_data(std::vector<std::valarray<double>>);
        void set_hyperparameters(unsigned int, 
                              unsigned int, 
                              double, 
                              double, 
                              double);
        /* genetic operators */
        std::vector<unsigned int> selection();
        std::vector<std::function<std::vector<unsigned int>(std::vector<unsigned int>)>> mutations;
        std::vector<unsigned int> mutation_perm(std::vector<unsigned int>);
        std::vector<unsigned int> mutation_shift(std::vector<unsigned int>);
        std::vector<unsigned int> mutation_swap(std::vector<unsigned int>);
        std::vector<unsigned int> mutation_inversion(std::vector<unsigned int>);
        void crossover(std::vector<unsigned int>&, std::vector<unsigned int>&);


        void evolve();

        void run();
        std::vector<unsigned int> get_solution();

        std::vector<double> L_best;
        std::vector<double> L_mean;

        double rand();
        int randint(int, int);

        std::vector<Individual> population;


    private:
        unsigned int NUM_CITIES;
        // std::valarray<std::valarray<double>> X; // (x,y) position of each city        
        
        /* hyperparameters */
        unsigned int POP_SIZE;
        unsigned int NUM_GENERATIONS;
        double MUTATION_PROB;
        double CROSSOVER_PROB;
        double POWER_LAW_EXPONENT;
      
        std::valarray<std::valarray<double>> cost_matrix;
    
        /* methods */
        double norm_s(std::valarray<double> x1, std::valarray<double> x2);

        /*PRNG*/
        // std::default_random_engine generator;
        std::random_device generator;
        std::uniform_real_distribution<double> distrib;
};


#endif // __Random__
