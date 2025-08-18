#ifndef SURFACE__H__
#define SURFACE__H__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "CylinderPotential.h"
#include "CubePotential.h"
#include "TubePotential.h"

class SurfaceData {
public:
    std::string         m_aminoAcid;
    std::vector<double> m_distance;
    std::vector<double> m_energy;

public:
    SurfaceData(const std::string& aminoAcid, const std::vector<double> distance, const std::vector<double>& energy)
        : m_aminoAcid(aminoAcid), m_distance(distance), m_energy(energy)
    {}
};

class SurfacePotential : public std::vector<double> {
private:
    double m_start, m_stop, m_dr;

public:
    SurfacePotential(const std::vector<double>& energy, const std::vector<double>& distance)
        : std::vector<double>(energy), m_start(distance[0]), m_stop(distance[distance.size() - 1]), m_dr((m_stop - m_start) / (energy.size() - 1.0))
    {}


    double Value(const double distance, const double radius, const double cutoff, int correctionType=1) {
  		if (distance > m_stop) {
       		return 0.0;
        }
       	else if (distance < m_start) {
                double distanceUnderflow = m_start - distance;
                
                double gradEst = (   (*this)[2] - (*this)[0] )/(2*m_dr);                
                
                //linear blocking
                //gradEst = std::min( gradEst, -100.0 );
                //gradEst = 0; 
                //double distUFPower = std::pow( distanceUnderflow, 4);  
       		//return (*this)[0] - distanceUnderflow*gradEst ;
                gradEst = std::min( gradEst, -100.0 );
                //repulsive r^-6 (roughly r^-12 integrated over 6 times)
                double aCoeff = (*this)[0] + m_start*gradEst/6.0 ;
                double bCoeff = -1 * gradEst * std::pow(m_start,7)/6.0 ; 
                double r = std::max( distance, 0.01 );
                //std::cout << distance << " " << m_start << " " << gradEst << " " << (*this)[0] << " coeffs:" << aCoeff << " " << bCoeff << ": "<<  aCoeff + bCoeff/(std::pow(r,6) ) << "\n";
                return aCoeff + bCoeff/(std::pow(r,6) );
                //double distUFPower = std::pow( distanceUnderflow, 6);

        }
        
        // Linear interpolation of the energy
        const int      index  = static_cast<int>((distance - m_start) * (this->size() - 1.0) / (m_stop - m_start));
        const double   factor = (distance - m_start) / m_dr - index;
        const double 	energy = (*this)[index + 1] * factor - (*this)[index] * (factor - 1.0);




        // Sphere correction factor
 


 double correction  = 1.0;
      if(radius > cutoff ){ //if the np bead itself is smaller than the cutoff distance, no correction is applied
            if(correctionType == 1){
            double top        = cutoff * cutoff * (distance - 2.0 * radius) + 2.0 * cutoff * distance * (distance - 2.0 * radius) - 3.0 * distance * distance * (distance + 2.0 * radius);
            double bottom     = 2.0 * (cutoff * cutoff + 2.0 * cutoff * distance + 3.0 * distance * distance) * (distance + radius);
            correction = -1.0 * top / bottom;
            }
            //The remaining potentials we use typically have at least one numerical integration involved and therefore we don't have simple analytical expressions. We also usually don't have limits for R->Infty and the expressions are not always numerically stable.
//The general method we use is as follows. An inclusive segment at the target radius is found by calculating the potential with no exclusion and subtracting the potential with the known exclusion. This is then divided by the inclusive segment obtained for a flat plane (for titania PMFs) or for a CNT at a radius of 0.75nm. At large distances this fails and so a ceiling is set on d to prevent the expressions diverging. 
            else if(correctionType == 2){ //cylinder.  
            double distCeilVal = cutoff;
              if(distance > distCeilVal){
                correction =  (  HamakerAtomCylinderUnit(radius,distCeilVal,0.000000011) - HamakerAtomCylinderUnit(radius,distCeilVal,cutoff)  )/ ( HamakerAtomCylinderUnit(1000,distCeilVal,0.000000011) - HamakerAtomCylinderUnit(1000,distCeilVal,cutoff) );
               }
              else{
                correction = (  HamakerAtomCylinderUnit(radius,distance,0.000000011) - HamakerAtomCylinderUnit(radius,distance,cutoff)  )/ (  HamakerAtomCylinderUnit(1000,distance,0.000000011) - HamakerAtomCylinderUnit(1000,distance,cutoff) );
              }
            }
            else if(correctionType == 3){ //cube. we assume that the exclusion distance is less than the half-length of the cube and because it's a flat surface there's no correction needed.
              correction = 1;
            }
            else if(correctionType == 4 || correctionType==5){//tube. here the PMFs are calculated for a radius of 0.75 and so we adapt the potentials based on this. As with cylinder, this is ratio of the inclusive segment for the actual radius divided by inclusive segment for the calculated tube.
              double distCeilVal = 0.9;
              if(distance > distCeilVal){
                correction = ( HamakerAtomTubeUnit(radius,distCeilVal,0.00001) -      HamakerAtomTubeUnit(radius,distCeilVal,cutoff)  )/( HamakerAtomTubeUnit(0.75,distCeilVal,0.00001) -      HamakerAtomTubeUnit(0.75,distCeilVal,cutoff)  );
              }
              else{
               correction = ( HamakerAtomTubeUnit(radius,distance,0.00001) -      HamakerAtomTubeUnit(radius,distance,cutoff)  )/( HamakerAtomTubeUnit(0.75,distance,0.00001) -      HamakerAtomTubeUnit(0.75,distance,cutoff)  );
              }
            }
      }
      
 
        // Multiply and return
        return energy * correction;
    }
};

class SurfacePMFs : public std::vector<SurfaceData> {
public:
    SurfacePMFs(const std::string& directory, const std::string& prefix, const std::vector<std::string>& aminoAcids) {
        ParseAllFilesInPMFDirectory(directory, prefix, aminoAcids);
    }

    void ParseAllFilesInPMFDirectory(const std::string& directory, const std::string& prefix,
        const std::vector<std::string>& aminoAcids) {
        if (aminoAcids.empty()) {
            std::cerr << "Error: No amino acid tags were found in the config file\n";
            std::exit(1);
        }
        
        for (const auto& aminoAcid : aminoAcids) {
            std::string filename = directory + "/" + prefix + aminoAcid + ".dat";

            std::ifstream handle(filename.c_str());

            if (!handle.is_open()) {
                std::cerr << "Error: Failed to find PMF file '" << filename << "'\n";
                std::exit(1);
            }

            std::string line;
            std::vector<double> energy, distance;

            try {

                while (std::getline(handle, line)) {
                    if (line.empty() || line[0] == '#') {
                        continue;
                    }
                    
                    std::size_t pos = line.find(',');

                    if (pos == std::string::npos) {
                        distance.emplace_back(std::stod(line.substr(0, 8)));
                        energy.emplace_back(std::stod(line.substr(8, 11)) * 0.40092); // kJ/mol -> kBT
                    }
                    else {
                        distance.emplace_back(std::stod(line.substr(0, pos)));
                        energy.emplace_back(std::stod(line.substr(pos + 1, line.size() - pos - 1)) * 0.40092); // kJ/mol -> kBT
                    }
                }
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Error: Failed to parse pmf file '" << filename << "'\n";
                std::cerr << "Error: Failed on line: " << line << "\n";
            }

            //energy = Average(energy);

            this->emplace_back(aminoAcid, distance, energy);

            handle.close();
        }
    }
};

#endif
