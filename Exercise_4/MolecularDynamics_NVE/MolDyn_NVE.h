/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//#include <string>
#include <array>
#include <vector>

const unsigned int ndim = 3; // dimension of our system (3D)
unsigned int n_particles;

std::vector< std::array<double,ndim> > q, v, q_old; // q stores the configuration for the n_particles particles, v their velocities

const std::string separator = "    ";

// measurable properties
const unsigned int n_properties = 5;
double epot, ekin, etot, temp, pres;
// E_pot, E_kin, E_tot, T, P
std::array<double, n_properties> properties,          // block-averaged values
                                 properties_err,      // their error
                                 properties_block{0},
                                 properties_cum_average{0},
                                 properties_cum_average2{0};

// input (and input-derived) variables
double vol, rho, box, rcut;
unsigned int nstep, nblocks, iprint, seed;
double delta;
bool initial_velocity_from_file, rescale_velocity;

// functions
void Input(const std::string inputfile);
void ConfFinal(void);
void VelocityFinal(void);
void ConfXYZ(int);
void SaveMeasure(void);
void MeasureBlock(int);
void measure(void);
void move(void);
std::array<double, ndim> force(unsigned int i_p);
std::array<double, ndim> pbc(std::array<double, ndim> r);
double norm(const std::array<double, ndim> array);
double norm_sq(const std::array<double, ndim> array);

// operators overloading
std::ifstream& operator>>(std::ifstream& input, std::array<double, ndim>& array);
std::ofstream& operator<<(std::ofstream& output, std::array<double, ndim>& array);
std::ostream& operator<<(std::ostream& output, std::array<double, ndim>& array);
std::array<double, ndim> operator*(const std::array<double, ndim>& A, double c);
std::array<double, ndim>& operator*=(std::array<double, ndim>& A, double c);
std::array<double, ndim>& operator+=(std::array<double, ndim>& A, const std::array<double, ndim> B);
std::array<double, ndim>& operator/=(std::array<double, ndim>& A, double c);
std::array<double, ndim> operator/(const std::array<double, ndim> A, double c);
std::array<double, ndim>& operator-=(std::array<double, ndim>& A, const std::array<double, ndim> B);
std::array<double, ndim> operator-(const std::array<double, ndim> A, const std::array<double, ndim> B);
std::array<double, ndim> operator+(const std::array<double, ndim> A, const std::array<double, ndim> B);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
