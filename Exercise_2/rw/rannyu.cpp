
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
   int M = 1000000;
if(argc > 1)
  	M = atoi(argv[1]);

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	ofstream out;
	out.open("random.dat");
   for(int i=0; i<M; i++){
      out << rnd.Rannyu_i(0, 6) << endl;
   }

   rnd.SaveSeed();
   return 0;
}
