#ifndef SURFACE__H__
#define SURFACE__H__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

std::vector<double> Average(const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for(int i = 0; i < (int) vec.size(); ++i) {
        int size = 0;
        double avg = 0.0;
        for(int j = 0; j < 21; ++j) {
            int k = i + j - 10;
            if (k < 0 || k >= (int) vec.size()) {
                continue;
            }
            ++size;
            avg += vec[k];
        }
        result[i] = avg / size;
    }
    return result;
}


class SurfaceData {
public:
    std::string         m_aminoAcid;
    std::vector<double> m_distance;
    std::vector<double> m_energy;

public:
    SurfaceData(const std::string& aminoAcid, const std::vector<double> distance, const std::vector<double>& energy)
        : m_aminoAcid(aminoAcid), m_distance(distance), m_energy(energy)
    {}

    std::vector<double> ShapeCorrectedEnergy(const double r, const double c) const {
        std::vector<double> correctedEnergy(m_energy.size());
        for (std::size_t i = 0; i < m_energy.size(); ++i) {
            const double d          = m_distance[i];
            const double top        = c * c * (d - 2.0 * r) + 2.0 * c * d * (d - 2.0 * r) - 3.0 * d * d * (d + 2.0 * r);
            const double bottom     = 2.0 * (c * c + 2.0 * c * d + 3.0 * d * d) * (d + r);
            const double correction = -1.0 * top / bottom;
            correctedEnergy[i]      = m_energy[i] * correction;
        }
        return correctedEnergy;
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
                    if (pos == std::string::npos) // .dat file
                    {
                        distance.emplace_back(std::stod(line.substr(0, 8)));
                        energy.emplace_back(std::stod(line.substr(8, 11)) * 0.40092); // kJ/mol -> kBT
                    }
                    else    // .csv file
                    {
                        distance.emplace_back(std::stod(line.substr(0, pos++)));
                        energy.emplace_back(std::stod(line.substr(pos, line.size() - pos)) * 0.40092); // kJ/mol -> kBT
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
