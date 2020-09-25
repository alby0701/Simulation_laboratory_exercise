#include <vector>
#include <valarray>
#include <algorithm>
#include <random>

#include <functional>

/* Individual of the population */
struct Individual {
    std::vector<unsigned int> sequence;
    double cost;
    bool operator> (const int &i, const Individual &p) {return i >  p.cost;}
    bool operator< (const int &i, const Individual &p) {return i <  p.cost;}
    
    Individual(){};
    Individual(std::vector<unsigned int> i_sequence): 
        {
            sequence = i_sequence;
            cost = TSP::total_cost(sequence);
        }
}

class TSP {
    public:
        /* methods */
        void build_cost_matrix(std::valarray<std::valarray<double>>);
        double get_unit_cost(unsigned int, unsigned int);
        double total_cost(std::vector<unsigned int>);
    
        std::vector<Individual> generate_initial_population();


    private:
        unsigned int NUM_CITIES;
        std::valarray<std::valarray<double>> X; // (x,y) position of each city
    
        std::vector<Individual> population;
        
        /* hyperparameters */
        unsigned int POP_SIZE;
        unsigned int NUM_GENERATIONS;
        double MUTATION_PROB;
        double CROSSOVER_PROB;
        double POWER_LAW_EXPONENT;
    
        std::vector<double> L_best;
        std::vector<double> L_mean;
    
        std::valarray<std::valarray<double>> cost_matrix;
    
        /* methods */
        double norm_s(std::valarray<double> x1, std::valarray<double> x2);
};

void TSP::set_data(std::valarray<std::valarray<double>> X) {
    NUM_CITIES = X.size();
}

void TSP::set_hyperparameters(std::valarray<std::valarray<double>> X) {
    /*********************
//     set population (POP_SIZE)
    L_best (NUM_GENERATIONS)
    L_mean (NUM_GENERATIONS)
    *********************/
    
}

double TSP::norm_s(std::valarray<double> x1, std::valarray<double> x2) {
    return std::pow(x2-x1, 2).sum();
}

void TSP::build_cost_matrix(std::valarray<std::valarray<double>> X) {
    for (unsigned int i = 0; i < NUM_CITIES - 1; i++) {
        cost_matrix[i] = std::valarray<double> (NUM_CITIES-1-i);
        for (unsigned int j = i + 1; j < NUM_CITIES - 1; j++) {
            cost_matrix[i][j] = norm_s(X[i], X[j]);
        }
    }
}

double TSP::get_unit_cost(unsigned int i, unsigned int j) {
    /* Return the entry (i,j) of the symmetric cost matrix */
    if (i < j) return cost_matrix[i][j];
    else return cost_matrix[j][i];
}
        
double TSP::total_cost(std::vector<unsigned int> sequence){
    double L = get_unit_cost(0, sequence.front()); /* first city is 0 */
    for (unsigned int i = 0; i < NUM_CITIES - 1; i++)
        L += get_unit_cost(sequence[i-1], sequence[i]);
    L += get_unit_cost(0, sequence.back());    
    return L;
}

std::vector<Individual> TSP::generate_initial_population() {
    std::vector<unsigned int> cities_range;
    for (int i=1; i<NUM_CITIES; ++i) cities_range.push_back(i);
   
    std::vector<Individual> population (POP_SIZE);
    for (unsigned int i = 0; i < POP_SIZE; i++) {
        population.push_back(Individual(std::random_shuffle(cities_range.begin(), cities_range.end()))); 
        std::push_heap(population.begin(), population.end());
    }
    std::sort_heap(population.begin(), population.end());
    return population;
}


/****** Genetic algorithms *******/
/* Selection */
/*
def selection(pop):
    """Select an individual from the ordered (by fitness) population `pop`"""
    j = round((POP_SIZE-1) * np.random.rand()**POWER_LAW_EXPONENT) # WRONG IN THE TEXT!!
    choices.append(j)
    return pop[j][1] # just the sequence

# mutation operators
def mutation_swap(sequence):
    """ pair permutation of cities (except for the first city)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 2, 4, 3, 5 \right]$
    """
    pos = np.random.randint(1,NUM_CITIES-2) # start position of the mutation
    
    mutated = list(sequence)
    mutated[pos] = sequence[pos+1]
    mutated[pos+1] = sequence[pos]
#     return (0, *np.random.permutation(sequence[1:-1]), 0) # keep first and last fixed
    return tuple(mutated)

def mutation_shift(sequence):
    """ shift of $+n$ positions for $m$ contiguous cities (except for the first city and $m \lt N-1$)
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a $+2$ shift of the second and third cities.
    """
    m = np.random.randint(1,NUM_CITIES-1) # cities affected by this mutation
    pos = np.random.randint(NUM_CITIES-1-m) # start position of the mutation
    
    shift = np.random.randint(NUM_CITIES-1)
    mutated = np.roll(sequence[1+pos:1+pos+m], shift)
    return (*sequence[:1+pos], *mutated, *sequence[1+pos+m:])

def mutation_permswap(sequence):
    """ permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2.
    """
    m = np.random.randint(1,NUM_CITIES//2-1) # cities affected by this mutation
    pos = np.random.randint(NUM_CITIES//2-1-m) # position of the mutation
    
    mutated1 = np.random.permutation(sequence[1+pos:1+pos+m])
    mutated2 = np.random.permutation(sequence[1+pos+m:1+pos+2*m])
    return (*sequence[:1+pos], *mutated1, *mutated2, *sequence[1+pos+2*m:])

def mutation_inversion(sequence):
    """ inversion of the order in which they appear in the path of $m$ cities (except for the first city and $m \le N$),
    e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 3, 2, 5 \right]$ for the inversion of the cities from 2 to 4.
    """
    m = np.random.randint(1,NUM_CITIES-1) # cities affected by this mutation
    pos = np.random.randint(NUM_CITIES-1-m) # start position of the mutation
    
    mutated = sequence[1+pos:1+pos+m][::-1]
    return (*sequence[:1+pos], *mutated, *sequence[1+pos+m:])

# crossover operator
def crossover(parent1, parent2):
    cut_pos = np.random.randint(1, NUM_CITIES-1)
    
    child1 = np.empty_like(parent1)
    child1[:cut_pos] = parent1[:cut_pos]
    child2 = np.empty_like(parent2)
    child2[:cut_pos] = parent2[:cut_pos]
    
    missing1 = set(parent1[cut_pos:-1])
    pos = cut_pos
    for i in parent2[1:-1]:
        if i in missing1:
            child1[pos] = i
            missing1.remove(i)
            pos += 1
    child1[-1] = 0
            
    missing2 = set(parent2[cut_pos:-1])
    pos = cut_pos
    for i in parent1[1:-1]:
        if i in missing2:
            child2[pos] = i
            missing2.remove(i)
            pos += 1
    child2[-1] = 0
    
    return tuple(child1), tuple(child2)
*/
/**************************************************/
                            
                            
/*                            
# evolve population
def evolve(pop):
    new_pop = []
    for i in np.arange(POP_SIZE//2): # each couple of parents gives two children
        # select two parents (which will become the offspring)
        parent1 = selection(pop)
        parent2 = selection(pop)
                
        if np.random.rand() < CROSSOVER_PROB:
            parent1, parent2 = crossover(parent1, parent2)
        
        # mutate children
        for mutation in [mutation_swap, mutation_shift, mutation_permswap, mutation_inversion]:
            if np.random.rand() < MUTATION_PROB:
                parent1 = mutation(parent1)
            if np.random.rand() < MUTATION_PROB:
                parent2 = mutation(parent2)
        
        heapq.heappush(new_pop, (cost(parent1), parent1))
        heapq.heappush(new_pop, (cost(parent2), parent2))
        
    new_pop.sort()
    return new_pop
    */
    
void TSP::run() {
    population = TSP::generate_initial_population();
    for (unsigned int i = 0; i < NUM_GENERATIONS; i++) {
        population = evolve(population);
        L_best[i] = population[0].cost;
//         L_mean[i] = np.mean([a[0] for a in pop[:POP_SIZE//2]])
    }
}

std::vector<unsigned int> TSP::get_solution() {
    return population[0].sequence; // TODO: add (0, ***, 0)
}