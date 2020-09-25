#include <cmath>
#include <vector>
#include <valarray>
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "tsp_ga.h"
    
int main(int argc, char* argv[]){
    if (argc < 3) { // We expect 3 arguments: the program name, the data_file path and the config_file path
        std::cerr << "Usage: " << argv[0] << " data_file config_file" << std::endl;
        return 1;
    }

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* hyperparameters */
    unsigned int POP_SIZE;
    unsigned int NUM_GENERATIONS;
    double MUTATION_PROB;
    double CROSSOVER_PROB;
    double POWER_LAW_EXPONENT;

    unsigned int NUM_MPI_GEN_MIGR;

    TSPga t;
    t.mpi_setup_worker(rank, size);

    if  (rank==0) {
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
        while (data_file >> x >> y)
            X.push_back(std::valarray<double> {x, y});
        t.set_data(X);

        config_file >> POP_SIZE;
        config_file >> NUM_GENERATIONS;
        config_file >> MUTATION_PROB;
        config_file >> CROSSOVER_PROB;
        config_file >> POWER_LAW_EXPONENT;
        config_file >> NUM_MPI_GEN_MIGR;
        
        t.set_hyperparameters(POP_SIZE, NUM_GENERATIONS, MUTATION_PROB, CROSSOVER_PROB, POWER_LAW_EXPONENT, NUM_MPI_GEN_MIGR);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t.mpi_sync_initial_data();

    t.run();
    t.get_solution();

    MPI_Finalize();
    return 0;
}
