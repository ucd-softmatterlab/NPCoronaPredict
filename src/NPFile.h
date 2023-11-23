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
    const std::vector<int>      m_npBeadType;
    const std::vector<double>   m_x;
    const std::vector<double>   m_y;
    const std::vector<double>   m_z;
    const std::vector<int>      m_id;
    const int                m_numBeadTypes;
    const std::string           m_name;
    const std::vector<double>   m_radius;
    const std::vector<double>   m_zeta;
    const std::vector<double>   m_coreFactor;
    const std::vector<double>   m_surfFactor;
    const std::vector<int>      m_shape;
    const std::vector<string>   m_hamakerFile;
    const std::vector<string>   m_pmfFile;
    const double                m_boundRadius;
    const double                m_outerBoundRadius;
    const double                m_zetaName;
    const std::vector<double>   m_pmfCutoff;
    const std::vector<int>      m_correctionType;
 
    NP(const std::vector<int>& npBeadType, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
    const std::vector<int>& id, const int numBeadTypes, const std::string& name, const std::vector<double>& radius,  
const std::vector<double>& zeta,const std::vector<double>& coreFactor,const std::vector<double>& surfFactor   ,const std::vector<int>& shape,   
const std::vector<string>& hamakerFile , const std::vector<string>& pmfFile , const double boundRadius, const double outerBoundRadius, const double zetaName, const std::vector<double>& pmfCutoff, std::vector<int>& correctionType )
        : m_npBeadType(npBeadType) ,  m_x(x), m_y(y), m_z(z), m_id(id), m_numBeadTypes(numBeadTypes), m_name(name), m_radius(radius), m_zeta(zeta), 
m_coreFactor(coreFactor), m_surfFactor(surfFactor), m_shape(shape), m_hamakerFile(hamakerFile), m_pmfFile(pmfFile), m_boundRadius(boundRadius), m_outerBoundRadius(outerBoundRadius),m_zetaName(zetaName),  m_pmfCutoff(pmfCutoff), m_correctionType(correctionType)
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
    
    //prepare the list of components, each of which is one of the bead types stored below
    std::vector<int> npBeadType;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    //list of NP bead types
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


//Old:
//NP definition: 1 line per NP, comma separated parameters in order:
//x,y,z,radius,zeta,coreFactor,surfFactor,shape,hamakerFile,pmfFile

//Updated:
//First part of the NP file defines bead types, second defines individual beads
//TYPE,radius,zeta,coreFactor,surfFactor,shape,hamaker,pmf,pmfcutoff,correctiontype
//BEAD,typeIndex,x,y,z

int numBeads = 0;
int numBeadTypes = 0;

double boundRadius = 0;
double outerBoundRadius = 0;
double zetaName = 0;

double fileInnerBound = -5;
double fileOuterBound = -5;

    while (std::getline(handle, line)) {
        if(line.size() > 3 && line.substr(0, 1) != "#") {
            try {
                

            
               std::vector<std::string> results;
                boost::split(results, line, [](char c){return c == ',';});

               if( results[0].substr(0,4) == "TYPE"){
               //register a new type
                double radiusval = std::stod(results[1]);
                double zetaval = std::stod(results[2]);
                   radius.emplace_back(radiusval);
          	   zeta.emplace_back(zetaval);
         	   coreFactor.emplace_back(std::stod(results[3]));
       		   surfFactor.emplace_back(std::stod(results[4]));
         	   shape.emplace_back(std::stoi(results[5]));
         	   hamakerFile.emplace_back(results[6]);
         	   pmfFile.emplace_back(results[7]);
       		   pmfCutoff.emplace_back(std::stod(results[8]));
		   correctionType.emplace_back(std::stoi(results[9]));  
                   numBeadTypes +=1;

                if( radiusval > boundRadius){
            boundRadius =  radiusval;
            }
                 
               if( pow(zetaval,2) > pow(zetaName,2) ){
            zetaName = zetaval;
            }                
                   
               
               }
               else if( results[0] == "INNERBOUND"){
               std::cout << "Inner bound manually set, will overwrite automatic and config-file options \n";
               fileInnerBound = std::stod(results[1]);
               }
               else if( results[0] == "OUTERBOUND"){
               std::cout << "Outer bound manually set, will overwrite automatic and config-file options \n";
               fileOuterBound = std::stod(results[1]);
               }
                              
               
               else if( results[0] == "BEAD"){
               //place an NP bead 
               int    npBeadTypeVal = std::stoi(results[1]);
               double xval = std::stod(results[2]);
               double yval = std::stod(results[3]);
               double zval = std::stod(results[4]);
               npBeadType.emplace_back(npBeadTypeVal);
               x.emplace_back(xval);
               y.emplace_back(yval);
               z.emplace_back(zval);
               numBeads +=1;
               double beadRadius  = 0.0; //assign a default value
               if( radius.size() > npBeadTypeVal    ){
               beadRadius = radius[npBeadTypeVal];
               }
               else{
               std::cout<<"NP bead of type " <<  npBeadTypeVal  << " not yet defined, ignoring for bounding radius. \n";
               }
               
               double newOuterBoundRadius   = sqrt( xval*xval + yval*yval + zval*zval   ) + beadRadius;
               
                if(newOuterBoundRadius > outerBoundRadius){
                outerBoundRadius = newOuterBoundRadius;
                }
            
            
               }
               else{
               std::cout << "NP input line " << line << " \n not recognised, will attempt to parse using legacy system. WARNING: Mixing this with new-style will lead to misassigned beads. \n";
               
               double xval = std::stod(results[0]);
               double yval = std::stod(results[1]);
               double zval = std::stod(results[2]);
               
               double radiusval = std::stod(results[3]);
                double zetaval = std::stod(results[4]);
                double coreFactorVal = std::stod(results[5]);
                double surfFactorVal = std::stod(results[6]);
                int shapeVal = std::stoi(results[7]) ;
                string hamakerFileS = results[8] ;
                std::string pmfFileS = results[9];
                double pmfCutoffVal = std::stod(results[10]) ;
                int correctionTypeVal = std::stoi(results[11]) ;
                int foundBeadType = 0;
                int npBeadTypeVal = 0;
                int j = 0;
                //test to see if the bead matches any beads which have already been seen
                for(j=0 ; j < radius.size(); ++j){
                
                 if( radius[j] == radiusval && zeta[j] == zetaval && coreFactor[j] == coreFactorVal && surfFactor[j] == surfFactorVal && shape[j] == shapeVal && hamakerFile[j] == hamakerFileS && pmfFile[j] == pmfFileS && pmfCutoff[j] == pmfCutoffVal && correctionType[j] == correctionTypeVal){
                 foundBeadType = 1;
                 npBeadTypeVal = j;
                 } 
                
                }
                
                if(foundBeadType == 0){
                

                //register a new bead type

                npBeadTypeVal = radius.size();
                   radius.emplace_back(radiusval);
          	   zeta.emplace_back(zetaval);
         	   coreFactor.emplace_back(coreFactorVal);
       		   surfFactor.emplace_back(surfFactorVal);
         	   shape.emplace_back(shapeVal);
         	   hamakerFile.emplace_back(results[8]);
         	   pmfFile.emplace_back(results[9]);
       		   pmfCutoff.emplace_back(pmfCutoffVal);
		   correctionType.emplace_back(correctionTypeVal);  
		   
		   std::cout << "Legacy reader: Defining bead type " << npBeadTypeVal << " radius " << radius[npBeadTypeVal] << "\n";
		   
		   numBeadTypes += 1;
		   
		                   if( radiusval > boundRadius){
            boundRadius =  radiusval;
            }
		   
		   
		   
               }
               
               //then add the NP bead
               std::cout<<"Legacy reader: Adding bead of type " << npBeadTypeVal << "\n";
               npBeadType.emplace_back(npBeadTypeVal);
               x.emplace_back(xval);
               y.emplace_back(yval);
               z.emplace_back(zval);
               numBeads +=1;
               double beadRadius  = 0.0; //assign a default value
               if( radius.size() > npBeadTypeVal    ){
               beadRadius = radius[npBeadTypeVal];
               }
               else{
               std::cout<<"NP bead of type " <<  npBeadTypeVal  << " not yet defined, ignoring for bounding radius. \n";
               }
               
               double newOuterBoundRadius   = sqrt( xval*xval + yval*yval + zval*zval   ) + beadRadius;
               
                if(newOuterBoundRadius > outerBoundRadius){
                outerBoundRadius = newOuterBoundRadius;
                }
               
               
               
               
               }


            

            
            
            //previous expression: sqrt( xval*xval + yval*yval + zval*zval   ) + radiusval;
            



            
       
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

    if(numBeadTypes > 64){
    std::cout << "Warning: more than 64 NP bead types have been registered. This will likely lead to a segfault. \n";
    }
    std::string name = NPTargetList::Filename(filename);
    std::cout << "NP filename: " << name << " generated with bounding radius " << boundRadius <<"\n";
    
    if(fileInnerBound > -1){
    boundRadius = std::max(0.0, fileInnerBound);
    std::cout << "Inner bound set to " << boundRadius << "\n";
    }
    if(fileOuterBound > -1){
    outerBoundRadius = std::max(0.0,fileOuterBound);
    std::cout << "Outer bound set to " << outerBoundRadius << "\n";
    }
    outerBoundRadius = std::max( boundRadius, outerBoundRadius);
    
    return NP(npBeadType, x, y, z, id, numBeadTypes, name,radius,zeta,coreFactor,surfFactor,shape,hamakerFile,pmfFile,boundRadius,outerBoundRadius,zetaName,pmfCutoff,correctionType);
}

#endif
