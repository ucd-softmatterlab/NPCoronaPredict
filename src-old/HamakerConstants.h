#ifndef HAMAKER_CONSTANTS__H__
#define HAMAKER_CONSTANTS__H__

#include "StringFormat.h"

#include <fstream>
#include <iostream>
#include <unordered_map>

class HamakerConstants {
private:
    std::unordered_map<std::string, double> m_aminoAcidToHamakerConstant;

public:
    HamakerConstants(const std::string& filename) {
        ReadHamakerFile(filename);
    }

    void ReadHamakerFile(const std::string& filename) {
        std::ifstream handle(filename.c_str());
        
        if (!handle.is_open()) {
            std::cerr << "Error: Could not find hamaker file '" << filename << "'\n";
            std::exit(1);
        }
        
        std::string line;
        
        while (std::getline(handle, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::string aminoAcid  = line.substr(0, 7);
            std::string kTValueStr = line.substr(7, 10);

            StringFormat::Strip(aminoAcid);
            StringFormat::Strip(kTValueStr);

            double kTValue;

            try {
                kTValue = std::stod(kTValueStr);
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Error: Failed to read hamaker constant from file '" << filename << "'\n";
                std::exit(1);
            }

            m_aminoAcidToHamakerConstant.insert(std::pair<std::string, double>(aminoAcid, kTValue));
        }
        
        handle.close();
    }

    double operator[] (const std::string& aminoAcid) const {
        if (!m_aminoAcidToHamakerConstant.count(aminoAcid)) {
            std::cerr << "Error: Failed to find a hamaker constant for the amino acid '" << aminoAcid << "'\n";
            std::exit(1);
        }
        return m_aminoAcidToHamakerConstant.at(aminoAcid);
    }
};

#endif
