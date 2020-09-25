/****************************************************************
 *****************************************************************
     _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
    _/_/  _/ _/       _/       Physics Department
   _/  _/_/    _/    _/       Universita' degli Studi di Milano
  _/    _/       _/ _/       Prof. D.E. Galli
 _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
#include "MolDyn_NVE.h"
#include <cmath>    // rint, pow
#include <fstream>  // Stream class to both read and write from/to files.
#include <iostream> // cin, cout: Standard Input/Output Streams Library
#include <stdlib.h> // srand, rand: to generate random number
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>
#include <numeric>
//#include <execution>

using namespace std;

int main(int argc, char *argv[]) {
    if (argc > 1)
        Input(argv[1]); // Inizialization
    else
        Input("input.dat");

    unsigned int steps_per_block = nstep / nblocks;
    for (unsigned int n = 0; n < nblocks; n++) {
        for (unsigned int s = 0; s < steps_per_block; s++) {
            move(); // Move particles with Verlet algorithm
            measure();
        }
        //ConfXYZ(n);
        if (n % 10 == 0)
            cout << "Block " << n + 10  << "/" << nblocks << endl;
        MeasureBlock(n); // compute mean and error for the block, reset block-specific var's
        SaveMeasure();
    }

    ConfFinal();     // Write final configuration to restart
    VelocityFinal(); // Write final velocities to restart
    return 0;
}

void Input(std::string inputfile) { // Prepare all stuff for the simulation
    ifstream ReadInput, ReadConf, ReadVelocity;
    // double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    // Set seed for random numbers
    srand(seed); // Initialize random number generator

    ReadInput.open(inputfile); // Read input

    ReadInput >> temp;
    ReadInput >> n_particles;
    cout << "Number of particles = " << n_particles << endl;
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)n_particles / rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol, 1.0 / 3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> nblocks;
    ReadInput >> iprint;
    ReadInput >> initial_velocity_from_file;
    ReadInput >> rescale_velocity;

    cout << "The program integrates Newton equations with the Verlet method" << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();

    // Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");

    // read particle initial configuration
    q.reserve(n_particles);
    std::ifstream config0("config.0");
    for (unsigned int i = 0; i < n_particles; ++i) {
        std::array<double,ndim> qi;
        config0 >> qi;
        q.push_back(qi * box);
    }

    // initial velocities
    v.reserve(n_particles);
    std::array<double,ndim> sum_v {0};
    if (initial_velocity_from_file) {
      cout << "Read initial velocity from file velocity.0 " << endl << endl;
      std::ifstream vel0("velocity.final");
      //std::ifstream vel0("velocity.0");
      for (unsigned int i = 0; i < n_particles; ++i) {
          std::array<double,ndim> vi;
          vel0 >> vi;
          sum_v += vi;
          v.push_back(vi);
      }
    } else {
      cout << "Prepare random velocities with center of mass velocity equal to zero" << endl << endl;
      for (unsigned int i = 0; i < n_particles; ++i) {
          std::array<double,ndim> vi;
          for (unsigned int d = 0; d < ndim; ++d)
              vi[d] = rand() - 0.5;
          sum_v += vi;
          v.push_back(vi);
      }
    }

    // rescale velocities to match given temperature (if instructed)
    if (rescale_velocity || !initial_velocity_from_file) { // do it anyway if we randomly generated the velocities
        sum_v /= (double)n_particles;
        double sum_v2(0);
        std::for_each(v.begin(), v.end(),
                      [sum_v, &sum_v2](std::array<double,ndim> &vi){
                        vi -= sum_v;
                        sum_v2 += norm_sq(vi);
                      });
        sum_v2 /= (double)n_particles;

        double scaling_factor = std::sqrt(3 * temp / sum_v2);
        std::for_each(v.begin(), v.end(),
                      [scaling_factor](std::array<double,ndim> &vi){
                        vi *= scaling_factor;
                      });
    }

    // compute q_old
    q_old = q;
    for (unsigned int i = 0; i < n_particles; ++i) {
        q_old[i] -= v[i] * delta;
    }

    // clear output files
    ofstream output;
    output.open("output.dat", ios::trunc);
    output.close();
}

void move(void) {
  // move particles using Verlet integration
    std::array<double,ndim> q_new;
    // compute forces here: q changes at every step of the next for loop
    std::vector< std::array<double,ndim> > forces;
    forces.reserve(n_particles);
    for (unsigned int i = 0; i < n_particles; ++i) {
        forces.push_back(force(i));
    }
    for (unsigned int i = 0; i < n_particles; ++i) {
        //q_new = pbc(q[i] * 2.0 - q_old[i] + force(i) * std::pow(delta, 2));
        q_new = pbc(q[i] * 2.0 - q_old[i] + forces[i] * std::pow(delta, 2));
        v[i] = pbc(q_new - q_old[i]) / (2.0 * delta);

        q_old[i] = q[i];
        q[i] = q_new;
    }
}

std::array<double, ndim> force(unsigned int i_p) {
    // compute forces acting on particle i_p as -Grad_ip V(r)
    std::array<double,ndim> F {0};
    for (unsigned int i = 0; i < n_particles; ++i) {
        if (i != i_p) {
            std::array<double,ndim> dvec = pbc(q[i_p] - q[i]);
            double dr = norm(dvec);

            if (dr < rcut)
              F += dvec * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
        }
    }
    return F;
}

void MeasureBlock(int block) {
    unsigned int steps_per_block = nstep / nblocks;

    for (unsigned int s = 0; s < properties.size(); ++s) {
        properties_block[s] /= (double)steps_per_block;
        properties_cum_average[s] += properties_block[s];
        properties_cum_average2[s] += std::pow(properties_block[s], 2);
        properties[s] = properties_cum_average[s] / (double)(block + 1);
        if (block != 0)
            properties_err[s] = std::sqrt((properties_cum_average2[s] / (double)(block + 1) - std::pow(properties[s], 2)) / (double)block);
    }

    // reset block variables
    for (unsigned int s = 0; s < properties.size(); ++s)
        properties_block[s] = 0;
}

void measure(){
    double V(0), V1(0);
    // cycle over pairs of particles
    for (unsigned int i = 0; i < n_particles - 1; ++i) {
        for (unsigned int j = i + 1; j < n_particles; ++j) {
            double dr = norm(pbc(q[i] - q[j]));
            if (dr < rcut) {
                V += 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);   // Potential energy
                V1 += 16.0 / pow(dr, 12) - 8.0 / pow(dr, 6); // modified potential energy (for the pressure)
            }
        }
    }

    // kinetic energy
    double t = 0.5 * std::accumulate(v.begin(), v.end(), 0.0, [](double s, std::array<double,ndim> x) {return s + norm_sq(x);});
    epot = V / (double)n_particles;                     // Potential energy
    ekin = t / (double)n_particles;                     // Kinetic energy
    temp = (2.0 / 3.0) * t / (double)n_particles;       // Temperature
    etot = (t + V) / (double)n_particles;               // Total enery
    pres = rho * temp + rho * V1 / (double)n_particles; // pressure

    // E_pot, E_kin, E_tot, T, P
    properties[0] = epot; // potential energy
    properties[1] = ekin; // kinetic energy
    properties[2] = etot; // total energy
    properties[3] = temp; // temperature
    properties[4] = pres; // pressure

    //std::for_each(properties.begin(), properties.end(),)
    for (unsigned int s = 0; s < properties.size(); ++s)
        properties_block[s] += properties[s];
}

void SaveMeasure() { // Save properties measurement
    std::ofstream output("output.dat", ios::app);
    for (unsigned int s = 0; s < properties.size(); ++s) {
        output << properties[s]    << separator
               << properties_err[s] << separator;
    }
    output << endl;
}

void ConfFinal(void) { // Write final configuration
    cout << endl << "Print final configuration to file config.final" << endl << endl;
    std::ofstream conffinal("config.final");
    for (auto qi : q) {
        std::array<double, ndim> qi1 = qi / box;
        conffinal << qi1 << endl;
    }
}

void VelocityFinal(void) { // Write final velocities
    cout << endl << "Print final velocities to file velocity.final" << endl << endl;
    std::ofstream velfinal("velocity.final");
    for (auto vi : v)
        velfinal << vi << endl;
}

void ConfXYZ(int nconf) { // Write configuration in .xyz format
    std::ofstream configxyz("frames/config_" + to_string(nconf) + "A.xyz");
    configxyz << n_particles << endl;
    configxyz << "Lennard-Jones potential" << endl;
    for (auto qi : q) {
        std::array<double, ndim> quiz = pbc(qi);
        configxyz << "LJ  ";
        //configxyz << pbc(qi) << endl;
        configxyz << quiz << endl;
    }
}

std::array<double, ndim> pbc(std::array<double, ndim> r) {
  // apply periodic boundary conditions with side L=box
    std::array<double, ndim> r_pbc;
    for (unsigned int d = 0; d < ndim; d++)
        r_pbc[d] = r[d] - box * rint(r[d] / box);
    return r_pbc;
}

double norm(const std::array<double, ndim> array) {
  //return sqrt(std::reduce(array.begin(), array.end(), 0, [](double s, double x) {return s + x*x;}));
    return std::sqrt(std::accumulate(array.begin(), array.end(), 0., [](double s, double x) {return s + x*x;}));
}

double norm_sq(const std::array<double, ndim> array) {
  //return sqrt(std::reduce(array.begin(), array.end(), 0, [](double s, double x) {return s + x*x;}));
    return std::accumulate(array.begin(), array.end(), 0., [](double s, double x) {return s + x*x;});
}

// overload some operators for convenience
std::ifstream& operator>>(std::ifstream& input, std::array<double, ndim>& array) {
    for (unsigned int d = 0; d < ndim; d++)
        input >> array[d];
    return input;
}
std::ofstream& operator<<(std::ofstream& output, std::array<double, ndim>& array) {
    for (unsigned int d = 0; d < ndim; d++)
        output << array[d] << separator;
    return output;
}
std::ostream& operator<<(std::ostream& output, std::array<double, ndim>& array) {
    for (unsigned int d = 0; d < ndim; d++)
        output << array[d] << separator;
    return output;
}
std::array<double, ndim> operator*(const std::array<double, ndim>& A, double c) {
    std::array<double, ndim> B;
    for (unsigned int d = 0; d < ndim; d++)
        B[d] = A[d] * c;
    return B;
}
std::array<double, ndim> operator-(const std::array<double, ndim> A, const std::array<double, ndim> B) {
    std::array<double, ndim> C;
    for (unsigned int d = 0; d < ndim; d++)
        C[d] = A[d] - B[d];
    return C;
}
std::array<double, ndim> operator+(const std::array<double, ndim> A, const std::array<double, ndim> B) {
    std::array<double, ndim> C;
    for (unsigned int d = 0; d < ndim; d++)
        C[d] = A[d] + B[d];
    return C;
}
std::array<double, ndim>& operator*=(std::array<double, ndim>& A, double c) {
    for (unsigned int d = 0; d < ndim; d++)
        A[d] *= c;
    return A;
}
std::array<double, ndim>& operator/=(std::array<double, ndim>& A, double c) {
    for (unsigned int d = 0; d < ndim; d++)
        A[d] /= c;
    return A;
}
std::array<double, ndim> operator/(const std::array<double, ndim> A, double c) {
  std::array<double, ndim> B;
    for (unsigned int d = 0; d < ndim; d++)
        B[d] = A[d] / c;
    return B;
}
std::array<double, ndim>& operator+=(std::array<double, ndim>& A, const std::array<double, ndim> B) {
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::plus<double>());
    /*for (unsigned int d = 0; d < ndim; d++)
      A[d] += B[d];*/
    return A;
}
std::array<double, ndim>& operator-=(std::array<double, ndim>& A, const std::array<double, ndim> B) {
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::minus<double>());
    /*for (unsigned int d = 0; d < ndim; d++)
      A[d] -= B[d];*/
    return A;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
