#ifndef PDB_FILE__H__
#define PDB_FILE__H__

#include "StringFormat.h"

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "RelaxClasses.h"

class PDB {
public:
    const std::vector<double>   m_x;
    const std::vector<double>   m_y;
    const std::vector<double>   m_z;
    const std::vector<int>      m_id;
    const double                m_length;
    const std::string           m_name;
    const std::vector<double>   m_occupancy;
    const std::vector<double>   m_rmsd;
    const std::vector<std::string>   m_resTag; 
    const std::vector<bond> m_bondSet;
    const std::vector<std::vector<int>> m_nbExclusions;
    const std::vector<int>      m_beadType;


    PDB(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
    const std::vector<int>& id, const double length, const std::string& name, const std::vector<double>& occupancy, const std::vector<double>& rmsd, 
        const std::vector<std::string>& resTag, const std::vector<bond> bondSet, const std::vector<std::vector<int>> nbExclusions, const std::vector<int> beadType)
        : m_x(x), m_y(y), m_z(z), m_id(id), m_length(length), m_name(name), m_occupancy(occupancy), m_rmsd(rmsd), m_resTag(resTag)  , m_bondSet(bondSet), m_nbExclusions(nbExclusions) , m_beadType(beadType) 
    {}
};

PDB ReadPDBFile(const std::string&, const std::unordered_map<std::string, std::size_t>&,     const std::unordered_map<std::string,std::string>&,  const std::vector<std::string>&,   const int , const double, const double ,const bool  , const double, const bool, const std::string,const double  );

class PDBs : public std::vector<PDB> {
public:
    PDBs(const std::vector<std::string>& filenames,
        const std::unordered_map<std::string, std::size_t>& aminoAcidIdMap, 
        const std::unordered_map<std::string, std::string>& ligandMap,const std::vector<std::string>& aaTags,

const int inputDisorderStrat = 0, const double disorderMin = -5.0, const double disorderMax = 50.0, const bool readLigands = false, const double bondCutoff=0.551, const bool backboneOn = false, const std::string backboneTag="GLY", const double backboneScaleFactor = 0.6666) {
        this->reserve(filenames.size());
        for (const auto& filename : filenames) {
            this->push_back(ReadPDBFile(filename, aminoAcidIdMap, ligandMap,aaTags, inputDisorderStrat, disorderMin, disorderMax, readLigands, bondCutoff, backboneOn, backboneTag, backboneScaleFactor));
        }
    }
};

bool isAA(std::string& tag,const  std::vector<std::string>& aaTags){

    //const std::vector<std::string> aaTags{"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",  "HIS", "ILE", "LEU", "LYS", "MET", "PHE",  "SER", "THR", "TRP", "TYR", "VAL", "HIE","HID","HIP","GAN"};
    if( std::find( aaTags.begin(), aaTags.end(), tag) == aaTags.end() ){
        return false;
    }
    else{
        return true;
    }
}

PDB ReadPDBFile(const std::string& filename, const std::unordered_map<std::string, std::size_t>& aminoAcidIdMap,   const std::unordered_map<std::string,std::string>& ligandMap,  const std::vector<std::string>& aaTags,    const int inputDisorderStrat = 0, const double disorderMin = -5.0, const double disorderMax = 50.0, const bool readLigands=false, const double bondCutoff=0.551, const bool backboneOn = false, const std::string backboneTag="GLY", const double backboneScaleFactor = 0.6666 ) {
    std::ifstream handle(filename.c_str());
    if (!handle.is_open()) {
        std::cerr << "Error: Could not find pdb file '" << filename << "'\n";
        std::exit(1);
    }

    std::string line;
    std::string tag;
    int disorderStrat = inputDisorderStrat;

    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<int>    id;
    std::vector<double> occupancy;
    std::vector<double> bfactor;
    std::vector<double> rmsd;
    std::vector<std::string> resTag; 
    std::vector<int>    beadType; //0 = CA, 1 =CB, 2 = Ligand 

    std::vector<bond> bondSet;
    std::vector< std::vector<int>  > nbExclusions;
    //double bondCutoff = config.m_bondCutoffNM;

    //bool backboneOn = true;
    //double backboneScaleFactor = 2.0/3.0; //scale backbone bead occupancies by this value to account  for replicated atoms
    bool bIsAlphaFold = false;
    while (std::getline(handle, line)) {
        bool foundAtom = false;
        bool atomIsLigand = false;
        bool atomIsBackboneGly = false;
        int currentBeadType = 0;

        if( line.find("ALPHAFOLD") != std::string::npos){
        //std::cout <<"ALPHAFOLD detected, assuming this is an AlphaFold structure \n";
        bIsAlphaFold = true;
        }

        //expected behaviour: a CA bead places that residue BUT if backbone is enabled AND the residue returns as a splittable AA, actually put a protein-backbone bead here

        if( line.substr(0, 4) == "ATOM" && line.substr(13, 2) == "CA" ){
        foundAtom = true;
        tag = line.substr(17, 3);
        StringFormat::Strip(tag);
        if( isAA(tag,aaTags) && backboneOn == true  ) { // note that by default isAA(GLY, aaTags) returns false to make sure we don't override actual GLY, as does isAA(PRO,aaTags) 
        //std::cout <<"placing backbone " << line << "\n";
        atomIsBackboneGly = true; //set that we're overriding this atom to be a backbone bead (defaults to glycine)
        tag = backboneTag ; //override to whatever out backbone atom type is - normally glycine - as the protein backbone if we have this mode enabled
        }

        //StringFormat::Strip(tag);

        }

        //if we're in backbone mode, put the remainder of the AA here
        if( line.substr(0, 4) == "ATOM" && line.substr(13, 2) == "CB" && backboneOn == true ){
        //foundAtom = true;
        tag = line.substr(17, 3);
        StringFormat::Strip(tag);
            if( isAA(tag,aaTags)  ){ 
             //std::cout << "placing residue " << line << "\n";
             foundAtom = true;
             currentBeadType = 1;
             //atomIsBackboneGly = true;
             }

        }



        if( line.substr(0,6) == "HETATM" && readLigands == true){
        foundAtom = true;
        atomIsLigand = true;
        
        //attempt to parse line as a potential ligand
                std::string ligandID = line.substr(17,3);
                std::string ligandAtomID = line.substr(12,4) ;
                StringFormat::Strip(ligandID);
                StringFormat::Strip(ligandAtomID);
               std::string trialTag = "";
               if( ligandMap.count( ligandID+"-"+ligandAtomID) == 1 ){
               trialTag = ligandMap.at(ligandID+"-"+ligandAtomID)  ;
              }
              if(trialTag != ""){ 
              tag = trialTag;
              //std::cout << "Parsing " << ligandID << "-" << ligandAtomID << " as " << tag << "\n";
              currentBeadType = 2;
              }
              else{
              foundAtom = false;
              }
  
              //std::cout << "Not Parsing " << ligandID << "-" << ligandAtomID << " as no CG substitute found found \n";
              
        }


        if(line.size() > 3 && foundAtom == true) {
            try {
                x.emplace_back(0.1 * std::stod(line.substr(30, 8)));
                y.emplace_back(0.1 * std::stod(line.substr(38, 8)));
                z.emplace_back(0.1 * std::stod(line.substr(46, 8)));
                //std::cout << 0.1 * std::stod(line.substr(30, 8)) << "," << 0.1 * std::stod(line.substr(38, 8)) << "," << 0.1 * std::stod(line.substr(46, 8))<<"\n";
                //std::cout << line << "\n";
                //std::cout << line.substr(56,4) << "\n";
                double occupancyVal =  std::stod(line.substr(54,6));
                if( atomIsBackboneGly == true){ 
                    occupancyVal = occupancyVal * backboneScaleFactor;
                }
                occupancy.emplace_back(  occupancyVal   );
                bfactor.emplace_back( std::stod(line.substr(60,6)));
                //tag = line.substr(16, 4); // When using lipids
                double rmsdVal = 0;
                if( bIsAlphaFold == false){
                rmsdVal = 0.1*sqrt(  std::stod( line.substr(60,6))/( 8.0 * 3.1415 * 3.1415) ) ;
                }
                else{
                rmsdVal = 0.05;
                }
                rmsd.emplace_back(rmsdVal);
                resTag.emplace_back( tag) ;
                beadType.emplace_back( currentBeadType);
                //std::cout << "read AA with bfactor" <<  std::stod(line.substr(60,6)) << " to  rmsd " << rmsdVal << "\n"; 
                //
                    

                id.emplace_back(aminoAcidIdMap.at(tag));
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Error: Failed to parse pdb file '" << filename << "'\n";
                std::cerr << "Error: Failed on line: " << line << "\n";
                std::exit(1);
            }
            catch (const std::out_of_range& oor) {
                std::cerr << "Error: Encountered a molecule in the pdb file which was not in the config file\n";
                std::cerr << "Error: pdb filename = '" << filename << "', molecule name = '" << tag << "'\n";
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
   
   //Strategies for dealing with disordered residues:
   //0   = default = do nothing
   //1 = shift to COM defined by ordered residues, leave active
   // 2 = shift to COM defined by ordered residues, set occupancy to 0
   //note that disordered is defined using AlphaFold convention for bfactor, pIDDT < 50 = disordered. If you're using non-AF proteins then this will give bad results so set disorderStrat = 0 for those.
    // Center protein
    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
    double totalOcc = 0;
    double numRes = 0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        int isDisordered = 0;
        if(bfactor[i]>disorderMin && bfactor[i] < disorderMax){
        isDisordered = 1;
        }
        if(disorderStrat == 2 && isDisordered == 1){
        occupancy[i] = 0.01;
        }
    
    
        if(disorderStrat == 0 ||  isDisordered == 0){ //if the residue is ordered or we don't care about disorder then it contributes to the COM
        mean_x += x[i]*occupancy[i];
        mean_y += y[i]*occupancy[i];
        mean_z += z[i]*occupancy[i];
        totalOcc += occupancy[i]; 
        }
        
        numRes +=1;

    }
   if(totalOcc > 0.5){
    mean_x /= totalOcc;
    mean_y /= totalOcc;
    mean_z /= totalOcc;
    }
    else{
    mean_x = 0;
    mean_y = 0;
    mean_z = 0;
    }
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
    
    
        int isDisordered = 0;
        if(bfactor[i]>disorderMin && bfactor[i] < disorderMax){
        isDisordered = 1;
        }
    
    
    if(disorderStrat == 0 || isDisordered == 0){ //apathetic to disorder or ordered: shift relative to COM
        x[i] -= mean_x;
        y[i] -= mean_y;
        z[i] -= mean_z;
    }
    else{ //disordered when we care about that: set to the COM
        x[i] = 0;
        y[i] = 0;
        z[i] = 0;
    }
      
    }
    
    // Lenght
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

    const double length = std::sqrt(max_x * max_x + max_y * max_y + max_z * max_z);

    std::string name = TargetList::Filename(filename);
    if(totalOcc < 0.5 * numRes){
        std::cout << numRes <<  " total CA, total occupancy: " << totalOcc << "\n";
    }


    int numBonds = 0;
    int size = static_cast<int>(x.size()) ; 
   //generate the list of bonds and exclusions for relaxation - we do it here so it's done once rather than 10000 times 
    for(int i = 0; i< size; ++i){
        std::vector<int> iExclusions;

        for(int j = i+1; j < size; ++j){

            double resDist = sqrt( pow(x[i]-x[j],2 )+pow(y[i]-y[j],2 )+pow(z[i]-z[j],2 ) );
            if(   resDist < bondCutoff){
              
                //
               //std::cout << i << "-" << j << "adding bond with length " << resDist << "\n" ;
                double bondRMS = 0.02; //default values
                double bondK = 1.0/pow(bondRMS,2);  //100.0/pow(resDist,2)

                if(resDist < 0.4){
                bondRMS = 0.005;  //extremely short bonds are tagged as backbone bonds and given a lower RMS
                bondK = 1.0/pow(bondRMS,2);
                }
                else{
                bondRMS = 0.5 * resDist;  //other bonds are assigned as being more flexible
                bondK =  100.0/pow( bondRMS, 2);
                }
                //double bondK = 1.0/pow(bondRMS,2);  //100.0/pow(resDist,2)
                bondSet.emplace_back(bond(i,j,resDist,  bondK  ));
                iExclusions.emplace_back( j );
                numBonds += 1;
                //std::cout << " Bonding residues: " << i << " and " << j << " distance: " << resDist << " bond const: " << pow(0.5*resDist,2) << "\n";
            }
        }

        nbExclusions.emplace_back( iExclusions );
    }





    return PDB(x, y, z, id, length, name,occupancy,rmsd,resTag, bondSet, nbExclusions, beadType);
}


//basic PDB file writer for getting optimised structures
void writePDB(std::string name, std::string directory, const int size, const PDB& pdb, double *x, double *y, double *z, double zOffset=0.0){
    boost::filesystem::create_directories(directory);
    std::string fileLoc = directory+"/"+name+".pdb";
    std::cout << "Saving optimised PDB to " << fileLoc << "\n";
    std::ofstream handle(fileLoc.c_str() ); 
    handle << "HEADER UA-OPTIMISED " << name << "\n";
    handle << "TITLE " << name << "\n";
    for(int i =0; i<size; ++i){
     handle << "ATOM  "; 
     handle <<  std::right << std::setw(5) <<  i ; 
     if( pdb.m_beadType[i] == 1 ){
     handle << std::left << "  CB" ; //place SCAs at the CB location for consistency 
     }
     else{
     handle << std::left << "  CA" ;
     }
     handle << "  ";
     //handle << "X";
     handle << pdb.m_resTag[i];
     handle << " A" ;
     handle << std::right << std::setw(4) << i ;
     handle << " "; //insertion ID
     handle << "   ";
     handle << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 10*x[i] ;
     handle << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 10*y[i] ;
     handle << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 10*z[i] ;
     handle << std::right << std::setw(6) << std::fixed << std::setprecision(2) << pdb.m_occupancy[i] ;
     double bfactor = 8.0*3.1415*3.1415*pow(10*pdb.m_rmsd[i],2);

     handle << std::right << std::setw(6) << std::fixed << std::setprecision(2) << bfactor ;
     handle << "      "; //blank space columns 67 - 72 
     handle << "    "; //segment ID
     handle << " C";
     handle << "\n";
    }
    handle << "END";
     handle.close();
}

#endif

