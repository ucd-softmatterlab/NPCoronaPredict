#ifndef NP_FILE__H__
#define NP_FILE__H__

#include "StringFormat.h"
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <unordered_map>

class NP {
public:
    const std::vector<double>   m_x;
    const std::vector<double>   m_y;
    const std::vector<double>   m_z;
    const std::vector<int>      m_id;
    const double                m_length;
    const std::string           m_name;
    const std::vector<double>   m_radius;
    const std::vector<double>   m_zeta;
    const std::vector<double>   m_coreFactor;
    const std::vector<double>   m_surfFactor;
    const std::vector<int>      m_shape;
    const std::vector<string>   m_hamakerFile;
    const std::vector<string>   m_pmfFile;
    const double                m_boundRadius;
    const std::vector<double>   m_pmfCutoff;
    const std::vector<int>      m_correctionType;
 
    NP(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
    const std::vector<int>& id, const double length, const std::string& name, const std::vector<double>& radius,  
const std::vector<double>& zeta,const std::vector<double>& coreFactor,const std::vector<double>& surfFactor   ,const std::vector<int>& shape,   
const std::vector<string>& hamakerFile , const std::vector<string>& pmfFile , const double boundRadius, const std::vector<double>& pmfCutoff, std::vector<int>& correctionType )
        : m_x(x), m_y(y), m_z(z), m_id(id), m_length(length), m_name(name), m_radius(radius), m_zeta(zeta), 
m_coreFactor(coreFactor), m_surfFactor(surfFactor), m_shape(shape), m_hamakerFile(hamakerFile), m_pmfFile(pmfFile), m_boundRadius(boundRadius), m_pmfCutoff(pmfCutoff), m_correctionType(correctionType)
    {}
};

NP ReadNPFile(const std::string&, const std::unordered_map<std::string, std::size_t>&);

class NPs : public std::vector<NP> {
public:
    NPs(const std::vector<std::string>& filenames,
        const std::unordered_map<std::string, std::size_t>& aminoAcidIdMap) {
        this->reserve(filenames.size());
        for (const auto& filename : filenames) {
            this->push_back(ReadNPFile(filename, aminoAcidIdMap));
        }
    }
};

NP ReadNPFile(const std::string& filename, const std::unordered_map<std::string, std::size_t>& aminoAcidIdMap) {
    std::ifstream handle(filename.c_str());
    if (!handle.is_open()) {
        std::cerr << "Error: Could not find NP file '" << filename << "'\n";
        std::exit(1);
    }

    std::string line;
    std::string tag;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<int>    id;
    std::vector<double> radius;
    std::vector<double> zeta;
    std::vector<double> coreFactor;
    std::vector<double>  surfFactor;
    std::vector<int>    shape;
    std::vector<string>  hamakerFile;
    std::vector<string>  pmfFile;
    std::vector<double>  pmfCutoff;
    std::vector<int>  correctionType;

//NP definition: 1 line per NP, comma separated parameters in order:
//x,y,z,radius,zeta,coreFactor,surfFactor,shape,hamakerFile,pmfFile
double length = 0;
double boundRadius = 0;


    while (std::getline(handle, line)) {
        if(line.size() > 3 && line.substr(0, 1) != "#") {
            try {
               std::vector<std::string> results;
                boost::split(results, line, [](char c){return c == ',';});

                //x.emplace_back(0.1 * std::stod(line.substr(30, 8)));
                //y.emplace_back(0.1 * std::stod(line.substr(38, 8)));
                //z.emplace_back(0.1 * std::stod(line.substr(46, 8)));
                //std::cout << line << "\n";
                //std::cout << line.substr(56,4) << "\n";
                //occupancy.emplace_back(  std::stod(line.substr(56,4)));
                //tag = line.substr(16, 4); // When using lipids
                //tag = line.substr(17, 3);
                //StringFormat::Strip(tag);

                //id.emplace_back(aminoAcidIdMap.at(tag));

            double xval = std::stod(results[0]);
            double yval = std::stod(results[1]);
            double zval = std::stod(results[2]);
            double radiusval = std::stod(results[3]);


            x.emplace_back(xval);
            y.emplace_back(yval);
            z.emplace_back(zval);
            radius.emplace_back(radiusval);
            zeta.emplace_back(std::stod(results[4]));
            coreFactor.emplace_back(std::stod(results[5]));
            surfFactor.emplace_back(std::stod(results[6]));
            shape.emplace_back(std::stoi(results[7]));
            hamakerFile.emplace_back(results[8]);
            pmfFile.emplace_back(results[9]);
            pmfCutoff.emplace_back(std::stod(results[10]));
            correctionType.emplace_back(std::stoi(results[11]));
            
            double newBoundRadius = 0;
            newBoundRadius = sqrt( xval*xval + yval*yval + zval*zval   ) + radiusval;

            if(newBoundRadius > boundRadius){
            boundRadius = newBoundRadius;
            }

            length +=1;
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Error: Failed to parse NP file '" << filename << "'\n";
                std::cerr << "Error: Failed on line: " << line << "\n";
                std::exit(1);
            }
            catch (const std::out_of_range& oor) {
                std::cerr << "Error: Encountered a molecule in the NP file which was not in the config file\n";
                std::cerr << "Error: NP filename = '" << filename << "', molecule name = '" << tag << "'\n";
                std::exit(1);
            }
        }
        else if (line.size() > 5 && line.substr(0, 6) == "ENDMDL") {
            break;
        }
    }

    handle.close();
     /*
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        std::cout << x[i] << " " <<  y[i] << " " << z[i] << "\n";
    }
   */
    // Center protein

/*
    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
    double totalOcc = 0;
    double numRes = 0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        mean_x += x[i];
        mean_y += y[i];
        mean_z += z[i];
        totalOcc += occupancy[i];
        numRes +=1;
        //std::cout << numRes << " " << occupancy[i] << "\n";

    }
   
    mean_x /= x.size();
    mean_y /= y.size();
    mean_z /= z.size();

    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        x[i] -= mean_x;
        y[i] -= mean_y;
        z[i] -= mean_z;
    }
   */ 
    // Lenght

   /*
    double max_x = 0.0, max_y = 0.0, max_z = 0.0;

    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (std::fabs(x[i]) > max_x) {
            max_x = std::fabs(x[i]);
        }
        if (std::fabs(y[i]) > max_y) {
            max_y = std::fabs(y[i]);
        }
        if (std::fabs(z[i]) > max_z) {
            max_z = std::fabs(z[i]);
        }
    }
   */
    //const double length = std::sqrt(max_x * max_x + max_y * max_y + max_z * max_z);
  
    std::string name = NPTargetList::Filename(filename);
    std::cout << "NP filename: " << name << " generated with bounding radius " << boundRadius <<"\n";
    return NP(x, y, z, id, length, name,radius,zeta,coreFactor,surfFactor,shape,hamakerFile,pmfFile,boundRadius,pmfCutoff,correctionType);
}

#endif
