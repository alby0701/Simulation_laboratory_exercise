#include <vector>
#include <valarray>
#include <algorithm>
#include <random>
#include <list>
#include <iostream>
#include "mpi.h"

#include "tsp_ga.h"

double TSPga::rand(){
    return distrib(generator);
}

int TSPga::randint(int max){
    return std::floor(max*rand());
}

int TSPga::randint(int min, int max){
    return std::floor(min+(max-min)*rand());
}

void TSPga::set_data(std::vector<std::valarray<double>> X) {
    NUM_CITIES = X.size();
    build_cost_matrix(X);
}

void TSPga::set_hyperparameters(unsigned int pop_size, 
                              unsigned int num_generations, 
                              double mutation_prob, 
                              double crossover_prob, 
                              double power_law_exponent,
                              unsigned int num_mpi_gen_migr) {
    POP_SIZE = pop_size;
    NUM_GENERATIONS = num_generations;
    MUTATION_PROB = mutation_prob;
    CROSSOVER_PROB = crossover_prob;
    POWER_LAW_EXPONENT = power_law_exponent;
    NUM_MPI_GEN_MIGR = num_mpi_gen_migr;
}

void TSPga::mpi_setup_worker(int rank, int size) {
    mpi_rank = rank;
    mpi_size = size;
}

void TSPga::mpi_sync_initial_data() {
    MPI_Bcast(&NUM_CITIES, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&POP_SIZE, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NUM_GENERATIONS, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&MUTATION_PROB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&CROSSOVER_PROB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&POWER_LAW_EXPONENT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NUM_MPI_GEN_MIGR, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mpi_rank > 0) cost_matrix.resize(NUM_CITIES-1);
    for (unsigned int i = 1; i < NUM_CITIES; i++) {
        if (mpi_rank > 0) cost_matrix[i-1].reserve(i);
        MPI_Bcast(&cost_matrix[i-1].front(), i, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

double TSPga::norm_s(std::valarray<double> x1, std::valarray<double> x2) {
    return std::pow(x2-x1, 2).sum();
}

void TSPga::build_cost_matrix(std::vector<std::valarray<double>> X) {
    // store as a lower triangular
    cost_matrix.resize(NUM_CITIES-1);
    for (unsigned int i = 1; i < NUM_CITIES; i++) {
        cost_matrix[i-1].resize(i);
        for (unsigned int j = 0; j < i; j++)
            cost_matrix[i-1][j] = norm_s(X[i], X[j]);
    }
}

double TSPga::get_unit_cost(unsigned int i, unsigned int j) {
    /* Return the entry (i,j) of the symmetric cost matrix */
    if (i > j) return cost_matrix[i-1][j];
    else return cost_matrix[j-1][i];
}
        
double TSPga::total_cost(std::vector<unsigned int> sequence){
    double L = get_unit_cost(0, sequence.front()); /* first city is 0 */
    for (unsigned int i = 0; i < NUM_CITIES - 2; i++)
        L += get_unit_cost(sequence[i], sequence[i+1]);
    L += get_unit_cost(0, sequence.back());
    return L;
}

void TSPga::generate_initial_population() {
    population.reserve(POP_SIZE);

    // std::vector<unsigned int> cities_range;
    // cities_range.reserve(NUM_CITIES-1);
    // for (unsigned int i=1; i<NUM_CITIES; ++i) cities_range.push_back(i);
    std::vector<unsigned int> cities_range(NUM_CITIES-1);
    for (unsigned int i=1; i<NUM_CITIES; ++i) cities_range[i-1] = i;
   
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
std::vector<unsigned int> TSPga::selection() {
    /* Select an individual from the ordered (by fitness) population */
    unsigned int j = std::round((POP_SIZE-1) * std::pow(rand(),POWER_LAW_EXPONENT)); // WRONG IN THE TEXT!!
    return population[j].sequence;
}

/* Mutations */
std::vector<unsigned int> TSPga::mutation_perm(std::vector<unsigned int> sequence) {
    /* pair permutation of cities (except for the first city)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 2, 4, 3, 5 \right]$ */
    unsigned int pos = randint(1,NUM_CITIES-2); // start position of the mutation
    std::swap(sequence[pos], sequence[pos+1]);
    return sequence;
}

std::vector<unsigned int> TSPga::mutation_shift(std::vector<unsigned int> sequence) {
    /* shift of $+n$ positions for $m$ contiguous cities (except for the first city and $m \lt N-1$)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a $+2$ shift of the second and third cities. */
    unsigned int aff = randint(1,NUM_CITIES-2); // cities affected by this mutation
    unsigned int start = randint(1,NUM_CITIES-2-aff);  // start position of the mutation
    unsigned int shift = randint(1,aff);

    std::rotate(sequence.begin()+start, sequence.begin()+start+shift, sequence.begin()+start+aff);
    return sequence;
}

std::vector<unsigned int> TSPga::mutation_swap(std::vector<unsigned int> sequence) {
    /* permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2. */
    unsigned int aff = randint(1,(NUM_CITIES-2)/2); // cities affected by this mutation
    unsigned int start1 = randint(1,(NUM_CITIES-2)/2-aff);  // start position of the mutation
    unsigned int start2 = randint(start1+aff,(NUM_CITIES-2)/2); // start position of the second part of the mutation

    for (unsigned int i = start1; i<start1+aff; i++)
        std::swap(sequence[i], sequence[start2+i]);
    return sequence;
}

std::vector<unsigned int> TSPga::mutation_inversion(std::vector<unsigned int> sequence) {
    /* permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2. */
    unsigned int m = randint(2,(NUM_CITIES-2)/2); // cities affected by this mutation
    unsigned int pos = randint(1,(NUM_CITIES-2)/2-m);  // start position of the mutation

    std::reverse(sequence.begin()+pos, sequence.begin()+pos+m);
    return sequence;
}


void TSPga::crossover(std::vector<unsigned int>& parent1, std::vector<unsigned int>& parent2) {
    std::vector<unsigned int> child1(NUM_CITIES-1);
    std::vector<unsigned int> child2(NUM_CITIES-1);

    unsigned int cut_pos = randint(1,(NUM_CITIES-2));
    std::copy(parent1.begin(), parent1.begin()+cut_pos, child1.begin());
    std::copy(parent2.begin(), parent2.begin()+cut_pos, child2.begin());

    std::list<unsigned int> missing1 (parent1.begin()+cut_pos, parent1.end());
    std::list<unsigned int> missing2 (parent2.begin()+cut_pos, parent2.end());

    unsigned int i = 0;
    for (auto p: parent2) {
        for (std::list<unsigned int>::iterator l = missing1.begin(); l != missing1.end(); l++) { 
            if (p == *l) {
                child1[cut_pos + i] = *l;
                l = missing1.erase(l);
                i++;
                break;
            }
        }
    }
    i = 0;
    for (auto p: parent1) {
        for (std::list<unsigned int>::iterator l = missing2.begin(); l != missing2.end(); l++) { 
            if (p == *l) {
                child2[cut_pos + i] = *l;
                l = missing2.erase(l);
                i++;
                break;
            }
        }
    }

    parent1 = std::move(child1);
    parent2 = std::move(child2);
}                          
                            
/* Evolution */
void TSPga::evolve() {
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

void TSPga::mpi_best_exchange() {
    // randomly swap the best individuals (without re-sorting) between workers
    std::vector<int> workers(mpi_size);
    // MPI_Gather(&full_best_sequence.front(), NUM_CITIES-1, MPI_INT, &dest, 1, MPI_INT, 0, MPI_COMM_WORLD);

/*
    std::cout << "W" <<  mpi_rank << ": sending to " << dest << std::endl;
    // MPI_Isend(&imesg[0], n, MPI_INTEGER,0, itag, MPI_COMM_WORLD, &req);
    MPI_Send(&mpi_rank, 1, MPI_INT, dest, dest, MPI_COMM_WORLD);//, &req);
    MPI_Recv(&dest2, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_rank, MPI_COMM_WORLD, &status);
    // std::cout << "W" <<  mpi_rank << ": received " << dest2 << std::endl;
    // MPI_recv(&, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    // MPI_Recv(&dest2, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat1);
    std::cout << "W" <<  mpi_rank << ": received from " << dest2 << std::endl;
*/

    for (unsigned int i = 0; i<POP_SIZE/20; i++)  {
        if (mpi_rank == 0) {
            for (int i=0; i<mpi_size; ++i) workers[i] = i;
            std::random_shuffle(workers.begin(), workers.end());
        }
        MPI_Status status;
        int dest; // the worker to whom we'll send our baby
        MPI_Scatter(&workers.front(), 1, MPI_INT, &dest, 1, MPI_INT, 0, MPI_COMM_WORLD); // each worker gets its destination
        std::vector<unsigned int> old_sequence (population[i].sequence);
        // randomly re-assigning best individual to another worker
        /*                            , DEST, TAG ,                            */
        MPI_Send(&old_sequence.front(), NUM_CITIES-1, MPI_UNSIGNED, dest, dest, MPI_COMM_WORLD);//, &req);
        MPI_Recv(&population[i].sequence.front(), NUM_CITIES-1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_rank, MPI_COMM_WORLD, &status);
    }

    std::sort(population.begin(), population.end());
}

    
void TSPga::run() {
    generate_initial_population();
    for (unsigned int i = 0; i < NUM_GENERATIONS; i++) {
        evolve();

        if (i % NUM_MPI_GEN_MIGR) {
            // MPI_Barrier(MPI_COMM_WORLD);
            mpi_best_exchange();
        }
    }
}

void TSPga::get_solution() {
    /* get only the best cost
    double finalres;
    MPI_Reduce(&population[0].cost, &finalres, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    std::cout << "W" << mpi_rank << ": "  << population[0].cost << std::endl;
    if (mpi_rank==0) std::cout << finalres << std::endl;
    */

    struct { 
        double cost; 
        int   rank; 
    } worker_cost, res_cost; // we need this to access the actual sequence (not just the cost) w/o actually passing all the sequences

    worker_cost.cost = population[0].cost;
    worker_cost.rank = mpi_rank;

    // MPI_Reduce(&worker_cost, &res_cost, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&res_cost.rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&worker_cost, &res_cost, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if (mpi_rank == res_cost.rank) {
        // std::cout << "cost: " << population[0].cost  << std::endl;
        std::cout << population[0].cost  << std::endl;
        std::cout << "0 ";
        for (auto p: population[0].sequence) std::cout << p << " ";
        std::cout << "0" << std::endl;
    }
}