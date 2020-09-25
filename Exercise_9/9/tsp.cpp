#include <vector>
#include <valarray>
#include <algorithm>
#include <random>
// #include <list>
#include <iostream>

#include "tsp.h"

double TSP::rand(){
    return distrib(generator);
}

int TSP::randint(int min, int max){
    return std::floor(min+(max-min)*rand());
}

void TSP::set_data(std::vector<std::valarray<double>> X) {
    NUM_CITIES = X.size();
    build_cost_matrix(X);
}

void TSP::set_hyperparameters(unsigned int pop_size, 
                              unsigned int num_generations, 
                              double mutation_prob, 
                              double crossover_prob, 
                              double power_law_exponent) {
    POP_SIZE = pop_size;
    NUM_GENERATIONS = num_generations;
    MUTATION_PROB = mutation_prob;
    CROSSOVER_PROB = crossover_prob;
    POWER_LAW_EXPONENT = power_law_exponent;                      

    L_best.resize(NUM_GENERATIONS);
    L_mean.resize(NUM_GENERATIONS);    
}

double TSP::norm_s(std::valarray<double> x1, std::valarray<double> x2) {
    return std::pow(x2-x1, 2).sum();
}

void TSP::build_cost_matrix(std::vector<std::valarray<double>> X) {
    // store as a lower triangular
    cost_matrix.resize(NUM_CITIES-1);
    for (unsigned int i = 1; i < NUM_CITIES; i++) {
        cost_matrix[i-1].resize(i);
        for (unsigned int j = 0; j < i; j++)
            cost_matrix[i-1][j] = norm_s(X[i], X[j]);
    }
}

double TSP::get_unit_cost(unsigned int i, unsigned int j) {
    /* Return the entry (i,j) of the symmetric cost matrix */
    if (i > j) return cost_matrix[i-1][j];
    else return cost_matrix[j-1][i];
}
        
double TSP::total_cost(std::vector<unsigned int> sequence){
    double L = get_unit_cost(0, sequence.front()); /* first city is 0 */
    for (unsigned int i = 0; i < NUM_CITIES - 2; i++)
        L += get_unit_cost(sequence[i], sequence[i+1]);
    L += get_unit_cost(0, sequence.back());
    return L;
}

void TSP::generate_initial_population() {
    population.reserve(POP_SIZE);

    std::vector<unsigned int> cities_range;
    for (unsigned int i=1; i<NUM_CITIES; ++i) cities_range.push_back(i);
   
    for (unsigned int i = 0; i < POP_SIZE; i++) {
        std::random_shuffle(cities_range.begin(), cities_range.end());        
        population.push_back(Individual(cities_range, total_cost(cities_range)));
        std::push_heap(population.begin(), population.end());
    }

    std::sort_heap(population.begin(), population.end());
    // for (auto p: population) std::cout << p.cost << " ";
    // std::cout << std::endl;
}


/****** Genetic operators *******/
/* Selection */
std::vector<unsigned int> TSP::selection() {
    /* Select an individual from the ordered (by fitness) population */
    unsigned int j = std::round((POP_SIZE-1) * std::pow(rand(),POWER_LAW_EXPONENT)); // WRONG IN THE TEXT!!
    return population[j].sequence;
}

/* Mutations */
std::vector<unsigned int> TSP::mutation_perm(std::vector<unsigned int> sequence) {
    /* pair permutation of cities (except for the first city)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 2, 4, 3, 5 \right]$ */
    unsigned int pos = randint(1,NUM_CITIES-2); // start position of the mutation
    std::swap(sequence[pos], sequence[pos+1]);
    return sequence;
}

std::vector<unsigned int> TSP::mutation_shift(std::vector<unsigned int> sequence) {
    /* shift of $+n$ positions for $m$ contiguous cities (except for the first city and $m \lt N-1$)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a $+2$ shift of the second and third cities. */
    unsigned int aff = randint(1,NUM_CITIES-2); // cities affected by this mutation
    unsigned int start = randint(1,NUM_CITIES-2-aff);  // start position of the mutation
    unsigned int shift = randint(1,aff);

    std::rotate(sequence.begin()+start, sequence.begin()+start+shift, sequence.begin()+start+aff);
    return sequence;
}

std::vector<unsigned int> TSP::mutation_swap(std::vector<unsigned int> sequence) {
    /* permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2. */
    unsigned int aff = randint(1,(NUM_CITIES-2)/2); // cities affected by this mutation
    unsigned int start1 = randint(1,(NUM_CITIES-2)/2-aff);  // start position of the mutation
    unsigned int start2 = randint(start1+aff,(NUM_CITIES-2)/2); // start position of the second part of the mutation

    for (unsigned int i = start1; i<start1+aff; i++)
        std::swap(sequence[i], sequence[start2+i]);
    return sequence;
}

std::vector<unsigned int> TSP::mutation_inversion(std::vector<unsigned int> sequence) {
    /* permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2. */
    unsigned int m = randint(2,(NUM_CITIES-2)/2); // cities affected by this mutation
    unsigned int pos = randint(1,(NUM_CITIES-2)/2-m);  // start position of the mutation

    std::reverse(sequence.begin()+pos, sequence.begin()+pos+m);
    return sequence;
}


void TSP::crossover(std::vector<unsigned int>& parent1, std::vector<unsigned int>& parent2) {
    std::vector<unsigned int> child1(NUM_CITIES-1);
    std::vector<unsigned int> child2(NUM_CITIES-1);

    unsigned int cut_pos = randint(1,(NUM_CITIES-2));
    std::copy(parent1.begin(), parent1.begin()+cut_pos, child1.begin());
    std::copy(parent2.begin(), parent2.begin()+cut_pos, child2.begin());

    // std::cout << cut_pos << std::endl;
    // std::cout << "p1 ";
    // for (auto p: parent1) std::cout << p << " ";
    // std::cout << std::endl;
    // std::cout << "p2 ";
    // for (auto p: parent2) std::cout << p << " ";
    // std::cout << std::endl;

    // std::list<unsigned int> missing1 (parent1.begin()+cut_pos+1, parent1.end());
    // std::list<unsigned int> missing2 (parent2.begin()+cut_pos+1, parent2.end());

    // std::cout << "m1 ";
    // for (auto p: missing1) std::cout << p << " ";
    // std::cout << std::endl;

    // unsigned int j = 0;
    // for (auto p: parent2) {
    //     for (std::list<unsigned int>::iterator i = missing1.begin(); i != missing1.end(); ++i) {
    //         if (p == *i) {
    //             std::cout << *i << " " << j << std::endl;
    //             child1[cut_pos + j] = *i;
    //             i++;
    //             j++;
    //             i = missing1.erase(i);
    //             break;
    //         }
    //     }
    // }

    std::vector<unsigned int> missing1 (parent1.begin()+cut_pos, parent1.end());
    std::vector<unsigned int> missing2 (parent2.begin()+cut_pos, parent2.end());

    unsigned int i = 0;
    for (auto p: parent2) {
        for (auto l: missing1) {
            if (p == l) {
                child1[cut_pos + i] = l;
                i++;
                break;
            }
        }
    }
    i = 0;
    for (auto p: parent1) {
        for (auto l: missing2) {
            if (p == l) {
                child2[cut_pos + i] = l;
                i++;
                break;
            }
        }
    }

    // for (auto p: child1) std::cout << p << " ";
    // std::cout << std::endl << std::endl;

    parent1 = std::move(child1);
    parent2 = std::move(child2);
}                          
                            
/* Evolution */
void TSP::evolve() {
    std::vector<Individual> new_population;
    new_population.reserve(POP_SIZE);
    for (unsigned int i = 0; i < POP_SIZE / 2; i++) {
        /*  select two parents (which will become the offspring) */
        std::vector<unsigned int> parent1 = selection();
        std::vector<unsigned int> parent2 = selection();

        if (rand() < CROSSOVER_PROB) {crossover(parent1, parent2);}

        /* mutate children */
        for (auto mutation: mutations) {
            if (rand() < MUTATION_PROB) parent1 = mutation(parent1);
            if (rand() < MUTATION_PROB) parent2 = mutation(parent2);
        }

        new_population.push_back(Individual(parent1, total_cost(parent1)));
        new_population.push_back(Individual(parent2, total_cost(parent2)));
    }
    
    std::sort_heap(new_population.begin(), new_population.end());
    population = std::move(new_population); // save new generation
    // population = new_population;
}
    
void TSP::run() {
    generate_initial_population();
    for (unsigned int i = 0; i < NUM_GENERATIONS; i++) {
        evolve();

        L_best[i] = population[0].cost;
        for (unsigned int j = 0; j < POP_SIZE / 2; j++) // mean cost of the best half of the population
            L_mean[i] += population[j].cost;
        L_mean[i] /= POP_SIZE / 2;
    }
}

std::vector<unsigned int> TSP::get_solution() {
    return population[0].sequence; // TODO: add (0, ***, 0)
}