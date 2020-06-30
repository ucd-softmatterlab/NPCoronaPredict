#include <cstdio>
#include <random>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iostream>
#include <vector>
//conversion utility for mapping the PMF from a tube onto a plane.
#include "TubePotential.h"
 
 using namespace std;
int main(int argc, char *argv[]){

 
            std::string filename = "ILE.dat";

            std::ifstream handle(filename.c_str());

            if (!handle.is_open()) {
                std::cerr << "Error: Failed to find PMF file '" << filename << "'\n";
                std::exit(1);
            }

            std::string line;
            std::vector<float> energy, distance;

            try {
                while (std::getline(handle, line)) {
                    if (line.empty() || line[0] == '#') {
                        continue;
                    }
 double d =  std::stod(line.substr(0, 8));
                        std::cout <<  d << " " << std::stod(line.substr(8, 11) )  << " " <<    ( HamakerAtomTubeUnit(1000,d,0.00001) -      HamakerAtomTubeUnit(1000,d,1.0)  )/( HamakerAtomTubeUnit(0.75,d,0.00001) -      HamakerAtomTubeUnit(0.75,d,1.0)  ) <<  "\n";
 
                     //   energy.emplace_back(std::stod(line.substr(8, 11)) * 0.40092); // kJ/mol -> kBT
 
                }
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Error: Failed to parse pmf file '" << filename << "'\n";
                std::cerr << "Error: Failed on line: " << line << "\n";
            }

 

            handle.close();


 
}

