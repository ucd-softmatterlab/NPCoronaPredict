#ifndef LIGAND_FILE__H__
#define LIGAND_FILE__H__

#include "StringFormat.h"
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <unordered_map>

//get a ligand ID to AA map 

class LigandMap {
private:
    std::unordered_map<std::string, int> m_ligandCount;

public:
    std::unordered_map<std::string, std::string> m_ligandAALookup;


    LigandMap(const std::string& filename) {
        ReadLigandFile(filename);
    }

    void ReadLigandFile(const std::string& filename) {
        std::ifstream handle(filename.c_str());
        
        if ( filename != "") {
        std::cout << "Beginning ligand parsing \n";
        if (!handle.is_open()) {
            std::cerr << "Error: Could not find ligand file '" << filename << " . To avoid using ligands, remove the config line or supply an empty file \n";
            std::exit(1);
        }
        
        std::string line;
        std::string lastLigand = "";
        bool isNewLigand = false;
        bool addLine = false;
        int lineNum = 1;
        while (std::getline(handle, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            //std::string aminoAcid  = line.substr(0, 7);
            //std::string kTValueStr = line.substr(7, 10);
            addLine = false;

            std::vector<std::string> results;
            boost::split(results, line, [](char c){return c == ',';});

            std::vector<std::string> results2;
            boost::split(results2, results[0], [](char c){return c == '-';});

            std::string ligandName = results2[0];
            std::string ligandID = results[0];
            std::string aaTag = results[1];
            StringFormat::Strip(ligandName);
            StringFormat::Strip(ligandID);
            StringFormat::Strip(aaTag);

            if( ligandName != lastLigand){
                //if the ligand name changes we assume the previous entry is complete, no further beads will be accepted for this ligand
                m_ligandCount.insert( std::pair<std::string, int>(lastLigand,1) );
                
                //std::cout << "ligand " << ligandName << " different to last  " <<  lastLigand << "\n";
                lastLigand = ligandName;

            }

            if( m_ligandCount.count(ligandName) == 0 ){
            //std::cout<< "Mapping " << ligandID<< " to CG bead " << aaTag << "\n";
            if( m_ligandAALookup.count(ligandID) == 0 ){
                       //std::cout<< "Mapping " << ligandID<< " to CG bead " << aaTag << "\n";

            m_ligandAALookup.insert(std::pair<std::string, std::string>(ligandID,aaTag));
            }
             else{
                
                std::cout  << "Warning on line "<< lineNum <<": " << "Mapping for " << ligandID << " already assigned to " << m_ligandAALookup.at( ligandID)  << ", skipping \n";
            }

          }
          else{
             std::cout << "Warning on line "<< lineNum <<": " <<  "An entry for " << ligandName << " has already been registered, skipping \n";
          }


        lineNum++;
        }
        
        handle.close();


         /* 
         //debug - print out all registered ligand to AA keys 
         for (const auto& [key, value] : m_ligandAALookup) {
           cout << key << ", " << value << '\n';
             }
         */

        }
        else{
        std::cout << "No ligand file specified, only ATOM CA records will be read \n";
        }
    }
 /*
    std::string operator[] (const std::string& ligandID) const {
        if (!m_ligandAALookup.count(ligandID)) {
            //std::cerr << "Error: Failed to find a hamaker constant for the amino acid '" << aminoAcid << "'\n";
            //std::exit(1);
            return ""; //if not found return an empty string, which is interpreted as "place no bead"
        }
        return m_ligandAALookup.at(ligandID);
}
*/


    
};

#endif
