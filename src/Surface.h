#ifndef SURFACE__H__
#define SURFACE__H__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

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


    double Value(const double distance, const double radius, const double cutoff) {
  		if (distance > m_stop) {
       		return 0.0;
        }
       	else if (distance < m_start) {
       		return (*this)[0];
        }
        
        // Linear interpolation of the energy
        const int      index  = static_cast<int>((distance - m_start) * (this->size() - 1.0) / (m_stop - m_start));
        const double   factor = (distance - m_start) / m_dr - index;
        const double 	energy = (*this)[index + 1] * factor - (*this)[index] * (factor - 1.0);

        // Sphere correction factor
        const double top        = cutoff * cutoff * (distance - 2.0 * radius) + 2.0 * cutoff * distance * (distance - 2.0 * radius)
            - 3.0 * distance * distance * (distance + 2.0 * radius);
        const double bottom     = 2.0 * (cutoff * cutoff + 2.0 * cutoff * distance + 3.0 * distance * distance) * (distance + radius);
        const double correction = -1.0 * top / bottom;

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
                        distance.emplace_back(std::stod(line.substr(0, 8)));
                        energy.emplace_back(std::stod(line.substr(8, 11)) * 0.40092); // kJ/mol -> kBT
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
