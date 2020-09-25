#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input();                                 //Inizialization
  Warmup(metro);
  for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk); //Reset block averages
    for (int istep = 1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk); //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}

////FUNCTIONS////
void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl
       << endl;
  cout << "Nearest neighbour interaction      " << endl
       << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl
       << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl
       << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> old_conf;

  ReadInput >> nblk;

  ReadInput >> nstep;
  if (old_conf == 1)
    cout << "The program utilizes the previous configuration" << endl;
  if (metro == 1)
    cout << "The program performs Metropolis moves" << endl;
  else
    cout << "The program performs Gibbs moves" << endl;
  ReadInput >> equilibration_step;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Number of equilibration steps = " << equilibration_step << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables
               //initial configuration chosen to be the previous one
  if (old_conf == 0)
  {
    ReadInput.open("config.final");
    for (int i = 0; i < nspin; ++i)
      ReadInput >> s[i];
    ReadInput.close();
  }
  //initial configuration if not chosen to be the previous one
  else
  {
    for (int i = 0; i < nspin; ++i)
    {
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
    }
  }

  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu] / (double)nspin << endl;
}

void Warmup(int metro)
{
  // Wait till the system is warmed up for MC
  //Read 2 seeds for random numbers
  int p1, p2, p3, p4;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 >> p3 >> p4;
  Primes.close();
  int seed2[4];
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();
  for (int i = 0; i < nspin; ++i)
  {
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
  }
  ofstream Ene;
  for (int istep = 1; istep <= equilibration_step; ++istep) {
    Move(metro);
    Measure();
    Ene.open("equi.mag.0", ios::app);
    Ene << walker[ix]/(double)nspin  << endl;
    Ene.close();
  }
  input.open("seed_2.in");
  input >> seed2[0] >> seed2[1] >> seed2[2] >> seed2[3];
  rnd.SetRandom(seed2, p3, p4);
  for (int i = 0; i < nspin; ++i)
  {
    if (rnd.Rannyu() >= 0.5)
      s[i] = 1;
    else
      s[i] = -1;
    }
  for (int istep = 1; istep <= equilibration_step; ++istep) {
  Move(metro);
  Measure();
  Ene.open("equi.mag.1", ios::app);
  Ene << walker[ix]/(double)nspin  << endl;
  Ene.close();
  }
}

void Move(int metro)
{
  for (int i = 0; i < nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    if (metro == 1) //Metropolis
    {
      int o = (int)(rnd.Rannyu() * nspin);
      attempted++;
      double alpha = min(1., exp(2 * beta * Boltzmann(s[o], o)));
      if (rnd.Rannyu() < alpha)
        s[o] = -s[o];
      else
        accepted++;
    }
    else //Gibbs sampling
    {
      double p = 1. / (1. + exp(-2 * beta * J * Boltzmann(1, i)));
      s[i] = rnd.Rannyu() > p ? 1 : -1;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
  return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0, m = 0.0, chiA = 0.0, chiB = 0.0, uA = 0.0, uB = 0.0;

  //cycle over spins
  for (int i = 0; i < nspin; ++i)
  {
    if (h != 0.0)
      m += s[i];
    if (h == 0.0)
    {
      chiA += s[i];
      chiB += s[i] * s[i];
      uB += (-J * s[i] * s[Pbc(i + 1)]) * (-J * s[i] * s[Pbc(i + 1)]);
      uA += -J * s[i] * s[Pbc(i + 1)];
    }
    u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
  }
  walker[iu] = u;
  walker[ic] = (double)nspin * beta * beta * (uB / (double)nspin - (uA / (double)nspin) * (uA / (double)nspin)); //beta*beta*(uA*uA - uB);
  walker[im] = m;
  walker[ix] = (beta * (chiA * chiA - chiB));

}

void Reset(int iblk) //Reset block averages
{

  if (iblk == 1)
  {
    for (int i = 0; i < n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) //Update block averages
{

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;

  //For block dimension verification
  ofstream Ene("block.ene.dat", ios::app);
  Ene << walker[iu]/ (double)nspin << endl;
  Ene.close();
}

void Averages(int iblk) //Print results for current block
{

/*
  ofstream Ene, Heat, Mag, Chi, Fin;
  const int wd = 12;

  cout << "Block number " << iblk << endl;
  if (metro == 1)
    cout << "Acceptance rate " << accepted / attempted << endl
         << endl;
  else
    cout << "Gibbs sampling selected " << endl
         << endl;
  Ene.open("output.ene.0", ios::app);
  */
  stima_u = blk_av[iu] / blk_norm / (double)nspin; //Energy
  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u * stima_u;
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);
  /*
  Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
  cout << "Energy status " << setw(wd) << glob_av[iu] / (double)iblk << " +- " << err_u << endl;
  Ene.close();*/
  if (h != 0.0)
  {

    //Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im] / blk_norm / (double)nspin; //Magnetization
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);
    /*Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
    cout << "Magnetization status " << setw(wd) << glob_av[im] / (double)iblk << " +- " << err_m << endl;
    Mag.close();*/
  }
  if (h == 0.0)
  {
    //Heat.open("output.heat.0", ios::app);
    stima_c = blk_av[ic] / blk_norm / (double)nspin; //Heat
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c * stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);
    /*Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
    cout << "Heat status " << setw(wd) << glob_av[ic] / (double)iblk << " +- " << err_c << endl;
    Heat.close();*/

    //Chi.open("output.chi.0", ios::app);
    stima_x = blk_av[ix] / blk_norm / (double)nspin; //Susceptibility
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x * stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
    /*Chi << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
    cout << "Heat status " << setw(wd) << glob_av[ic] / (double)iblk << " +- " << err_c << endl;
    Chi.close();*/

  }
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl
       << endl;
  WriteConf.open("config.final");
  for (int i = 0; i < nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
  ofstream Fin;
  const int wd = 12;
  Fin.open("outputFINAL.dat");
  Fin << glob_av[iu] / (double)nblk << setw(wd) << err_u << endl;
  Fin << glob_av[im] / (double)nblk << setw(wd) << err_m << endl;
  Fin << glob_av[ic] / (double)nblk << setw(wd) << err_c << endl;
  Fin << glob_av[ix] / (double)nblk << setw(wd) << err_x << endl;
  Fin.close();

  rnd.SaveSeed();
}

int Pbc(int i) //Algorithm for periodic boundary conditions
{
  if (i >= nspin)
    i = i - nspin;
  else if (i < 0)
    i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk == 1)
    return 0.0;
  else
    return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
}
