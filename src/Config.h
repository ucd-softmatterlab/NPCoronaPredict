#ifndef CONFIG__H__
#define CONFIG__H__

#include "StringFormat.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

class Config {
public: // Switches
    bool        m_enableSurface         = false;
    bool        m_enableCore            = false;
    bool        m_enableElectrostatic   = false;
    bool	m_sumNPPotentials    = true;
    bool        m_enableFullScan     = false;
    bool        m_zshiftToPlane = false;
    bool        m_enableLocalBoltz = false;
    bool        m_confirmOverrideAngle = false;
    bool        m_readLigands = false;
    bool        m_relaxPDB = false;
    bool        m_saveOptimumPDB = false;
public: // Key - vaules
    std::vector<std::string>    m_pdbTargets        = {};
    std::vector<std::string>    m_npTargets         = {};
    std::vector<std::string>    m_aminoAcids        = {};
    std::vector<double>         m_aminoAcidRadii    = {};
    std::vector<double>         m_nanoparticleRadii = {};
    std::vector<double>         m_zetaPotential     = {};
    std::vector<double>         m_aminoAcidCharges  = {};
    std::vector<double>         m_testAngle         = {};
    std::vector<double>         m_omegaAngles       = {};
    std::string m_configFile            = "ua.config";
    std::string m_pmfDirectory          = ".";
    std::string m_pmfPrefix             = "";
    std::string m_hamakerFile           = "hamaker.dat";
    std::string m_outputDirectory       = "";
    std::string m_ligandFile            = "";
    int         m_simulationSteps       = 10000;
    int         m_potentialSize         = 2000;
    int         m_multiNP               = 0;
    int         m_numSamples            = 128;
    int         m_flexMethod            = 0;
    int         m_flexNumSDev           = 3;
    double      m_flexResolution        = 0.02;

    double      m_potentialCutoff       = 10.0;
    double      m_angleDelta            = 5.0;
    double      m_bejerumLength         = 1.0;
    double      m_debyeLength           = 0.7;
    double      m_imaginary_radius      = -1.0;
    double      m_PMFCutoff      = 1.0;
    int		m_npType		= 1; //this defines the type of the nanoparticle. 1 = sphere, 2 = solid cylinder, 3 = cube, 4 = tube (hollow cylinder)
    int         m_recalcZP               = 0; //if this is non-zero then the input zeta potentials are treated as reference values for an NP of radius 1 in solution with Debye and Bjerrum lengths equal to one, and used to calculate shape- and size-dependent ZPs.
    int         m_calculateMFPT = 0 ; //if this is non-zero then UA also calculates and outputs the mean first passage time for each orientation. This involves a bunch of integration and so is slow. 
   int         m_savePotentials = 0; //if this is set to 1 then UA saves a potential file for each potential and orientation
    double     m_temperature = 300.0; //set the temperature in Kelvin. this is used as a scale for the PMFs
    double m_potentialStart = 0.001; //closest approach to bounding surface
    double m_boundingRadius = -1; //if  > 0: bounding radius (effective NP radius), if less than 0 automatically estimates this parameter using a safest value
   double m_overlapPenalty = 0.0; //Penalty for overlapping NP - AA beads
   double m_overlapRadiusFactor = 1.0;
    int         m_disorderStrat = 0; //0 = ignore disorder, 1 = shift disordered atoms to COM of remaining atoms, 2 = shift and disable
  //Disordered residues are those falling between disorderMinBound and disorderMaxBound given below. 
  //For AlphaFold, setting disorderMinBound less than 0 and disorderMaxBound to 50 gets the suggested limits.
  //For physical proteins see e.g. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2083-8 for a discussion of what a>
  double m_disorderMinBound = -5; //b-factors less than this value are not recognised as disordered
  double m_disorderMaxBound = 50; //b-factors greater than this value are not recognised as disordered.
  double m_pdbJitterMag = 0.0; //set the magnitude for PDB jittering in nanometers [recommended < 0.1 nm] 
  double m_potNoiseMag = 0.0;
  double m_bondCutoffNM = 0.8; // 0.551; //default value for the bond cutoff length for relaxation
   int m_relaxSteps = 400;
   int m_rigidSteps = 20000;
   double m_relaxZPull = 0.0; //gradient of the potential used to pull during the initial relaxation, positive values = pulling to - z
    //parameters enabling backbone-mode
    bool        m_enableBackbone        = false;
    std::string m_backboneTag           = "PBB";
    double      m_backboneScale         = 1.0; //scaling factor for backbone potentials 
    //default residues to split into backbone and SCA components 
    std::vector<std::string>    m_backboneReplaceSet        = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",  "HIS", "ILE", "LEU", "LYS", "MET", "PHE",  "SER", "THR", "TRP", "TYR", "VAL", "HIE","HID","HIP","GAN"};

public:
    void UpdateSwitches(const std::vector<std::string>& switches) {
        for (std::size_t i = 0; i < switches.size(); ++i) {
            if (switches[i] == "enable-surface") {
                m_enableSurface = true;
            }
            else if (switches[i] == "enable-core") {
                m_enableCore = true;
            }
            else if (switches[i] == "enable-electrostatic") {
                m_enableElectrostatic = true;
            }
            else if(switches[i] == "disable-np-summation"){
            m_sumNPPotentials = false;
            }
            else if(switches[i] == "enable-fullscan"){
            m_enableFullScan = true;
            }
            else if(switches[i] == "enable-local-boltz"){
            m_enableLocalBoltz = true;
            }           
            else if(switches[i] == "enable-pdb-relax"){
             m_relaxPDB = true;
            }
            else if(switches[i] == "save-optimum-pdb"){
            m_saveOptimumPDB = true;
           }

            else if(switches[i] == "enable-backbone"){
            m_enableBackbone = true; 
            std::cout << "Enabling backbone \n" ; 
            }
            else if(switches[i] == "confirm-override-angle"){
            std::cout << "Warning: Using angular resolutions other than 5 degrees is not supported by the NPCoronaPredict pipeline. Higher resolutions do not correspond to realistic conformations and the resulting energies must be interpreted with extreme caution.  \n";
            m_confirmOverrideAngle = true;
            }

            else if(switches[i] == "disable-protein-smartzshift"){
            std::cout << "Warning: disable-protein-smartzshift should not be used for production runs and is provided only for testing purposes.\n" ;
            m_zshiftToPlane = true;
            }
            else {
                std::cerr << "Error: Unknown switch statment '" << switches[i] << "\n'";
                std::exit(1);
            }
        }
    }

    void UpdateKeyValues(const std::vector<std::string>& keys, const std::vector<std::string>& values) {
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (keys[i] == "config-file") {
                m_configFile = values[i];
            }
            else if (keys[i] == "nanoparticle-radius") {
                m_nanoparticleRadii = AsDoubleList(values[i]);
            }
            
            
           else if (keys[i] == "omega-angles") {
                m_omegaAngles = AsDoubleList(values[i]);
            }
            
            else if (keys[i] == "ligand-file"){
                 m_readLigands = true;
                 m_ligandFile =  values[i] ;
            }
            
            else if (keys[i] == "amino-acids") {
                m_aminoAcids = AsStringList(values[i]);
            }
            else if (keys[i] == "amino-acid-radii") {
                m_aminoAcidRadii = AsDoubleList(values[i]);
            }
            else if (keys[i] == "pdb-target") {
                m_pdbTargets = AsStringList(values[i]);
            }
            else if(keys[i] == "np-target"){
                m_npTargets = AsStringList(values[i]); 
                m_multiNP = 1;
            }

            else if(keys[i] == "replace-sca-set"){
               m_backboneReplaceSet = AsStringList(values[i]) ;
            }

            else if (keys[i] == "pmf-directory") {
                m_pmfDirectory = values[i];
            }
            else if (keys[i] == "output-directory") {
                m_outputDirectory = values[i];
            }
            else if (keys[i] == "pmf-prefix") {
                m_pmfPrefix = values[i];
            }
            else if (keys[i] == "hamaker-file") {
                m_hamakerFile = values[i];
            }
            else if (keys[i] == "simulation-steps") {
                m_simulationSteps = AsInt(values[i]);
            }
            else if (keys[i] == "potential-size") {
                m_potentialSize = AsInt(values[i]);
            }
            else if (keys[i] == "potential-cutoff") {
                m_potentialCutoff = AsDouble(values[i]);
            }
            else if (keys[i] == "angle-delta") {
                m_angleDelta = AsDouble(values[i]);
            }
            else if (keys[i] == "num-random-samples"){
                int numSamplesIn = AsInt(values[i]);
                if( numSamplesIn < 1){
                std::cout << "Samples must be greater than 0, setting to 1 \n";
                numSamplesIn = 1;
                }
                m_numSamples = numSamplesIn;
            }
            else if (keys[i] == "pmf-cutoff") {
                m_PMFCutoff = AsDouble(values[i]);
            }
            else if (keys[i] == "bjerum-length") {
                m_bejerumLength = AsDouble(values[i]);
            }
            else if (keys[i] == "debye-length") {
                m_debyeLength = AsDouble(values[i]);
            }
            else if (keys[i] == "zeta-potential") {
                m_zetaPotential = AsDoubleList(values[i]);
            }
            else if (keys[i] == "temperature"){
                m_temperature = AsDouble(values[i]);
           }
            else if (keys[i] == "amino-acid-charges") {
                m_aminoAcidCharges = AsDoubleList(values[i]);
            }
            else if (keys[i] == "imaginary-radius") {
                m_imaginary_radius = AsDouble(values[i]);
            }
            else if (keys[i] == "bounding-radius"){
               m_boundingRadius = AsDouble(values[i]);
           }
           else if (keys[i] == "overlap-penalty"){
                 m_overlapPenalty = AsDouble(values[i]);
            }   
           else if (keys[i] == "overlap-radiusfactor"){
                 m_overlapRadiusFactor = AsDouble(values[i]);
                 if(  m_overlapRadiusFactor < 0){
                   std::cout << "Warning: overlap radius factors < 0 will lead to unphysical behaviour. If you're trying to disable overlap set overlap-penalty less than zero. \n";
                  }
            }  
            else if (keys[i] == "test-angle") {
                m_testAngle = AsDoubleList(values[i]);
                if (m_testAngle.size() != 2) {
                    std::cerr << "Error: Must provide two test angles [ phi, theta ]: " << values[i] << "\n";
                    std::exit(1);
                }
            }
            else if(keys[i] == "recalculate-zp"){
            m_recalcZP = AsInt(values[i]);
               if(m_recalcZP!=0){
              std::cout << "Treating input zeta potential as reference \n";
              }
              else{
              std::cout << "Using input ZP as actual value \n";
             }
            }
            else if(keys[i] == "calculate-mfpt"){
            m_calculateMFPT = AsInt(values[i]);
              std::cout << "Calculating MFPT \n";
            }
           else if( keys[i] == "save-potentials"){
            m_savePotentials =AsInt(values[i]) ;
             std::cout <<"Enabled saving UA potentials \n";
            }

            else if(keys[i] == "bead-flexibility-method"){
                int trialFlexMethod = AsInt(values[i]) ;

               if(trialFlexMethod == 0){
                 std::cout << "Flex method 0: Biobeads are fixed in place \n ";
                 m_flexMethod = 0;
               }
               else if(trialFlexMethod == 1){
                 std::cout << "Flex method 1: Gaussian smoothing\n";
                 m_flexMethod = 1;
               }
               else if(trialFlexMethod == 2){
                 std::cout << "Flex method 2: Gaussian re-sampling of energy\n";
                 m_flexMethod = 2;
               }
               else if(trialFlexMethod == 3){
                 std::cout << "Flex method 3: Free energy (gaussian re-sampling of probability, per-potential) \n";
                 m_flexMethod = 3;
               }
               else if(trialFlexMethod == 4){
                 std::cout << "Flex method 4: Free energy (gaussian re-sampling of probability, whole molecule) \n";
                 m_flexMethod = 4;
               }


               else{
                  m_flexMethod = 0;
                  std::cout << "Using default: fixed beads \n";
               }


           }
           else if( keys[i] == "flex-sdev"){
            m_flexNumSDev =std::max(0, AsInt(values[i]) );
            
             std::cout <<"Using "<< m_flexNumSDev << "standard deviations for flexibility \n";
            }
           else if( keys[i] == "flex-resolution"){
            m_flexResolution = std::max(0.01, AsDouble(values[i])) ;
             std::cout <<" Flex resolution set to " << m_flexResolution << "  \n";
            }
            else if( keys[i] == "bond-cutoff"){
             m_bondCutoffNM = std::max(0.0, AsDouble(values[i]) );
              if(  m_relaxPDB == true){
                std::cout << "Bond cutoff for relaxation set to " << m_bondCutoffNM << "\n";
                }
            }
            else if( keys[i] == "relax-steps"){
            m_relaxSteps = std::max(0, AsInt(values[i]) );
             if( m_relaxPDB == true){
              std::cout << "Using " << m_relaxSteps << " steps to relax structure \n";
              }
           }
           else if (keys[i] == "relax-gradient"){
            m_relaxZPull = AsDouble(values[i]) ;
              if( m_relaxPDB == true ){
              std::cout << "Potential gradient (-force) for relaxation set to " << m_relaxZPull << "\n";
              }
           }
          else if (keys[i] == "rigid-dynamics-steps"){
        
         m_rigidSteps  = AsDouble(values[i]) ;
          if( m_relaxPDB == true){
                 std::cout << "Using " << m_rigidSteps << " for rigid dynamics \n";
             } 
          }

           else if( keys[i] == "disorder-strategy"){
            m_disorderStrat =AsInt(values[i]) ;
            if(m_disorderStrat == 1){
             std::cout <<"Disorder strategy: shift to COM \n";
              }
              else if(m_disorderStrat == 2){
            std::cout <<"Disorder strategy: shift to COM and ignore \n";
              }
            else{
                    std::cout <<"Disorder strategy: proceed as normal \n";    
            m_disorderStrat = 0;
            }
            }
            else if(keys[i] == "disorder-minbound"){
            m_disorderMinBound = AsDouble(values[i]);
            }
            else if(keys[i] == "disorder-maxbound"){
            m_disorderMaxBound = AsDouble(values[i]);
            }
            else if(keys[i] == "pdb-jitter-magnitude"){
            m_pdbJitterMag = AsDouble(values[i]);
            }
            else if(keys[i] == "potential-noise-magnitude"){
            m_potNoiseMag = AsDouble(values[i]);
            }
            else if(keys[i] == "backbone-scale"){
            m_backboneScale = AsDouble(values[i]);
            std::cout << "Backbone beads will be scaled by " << m_backboneScale << "\n";
            }
            else if(keys[i] == "backbone-tag"){
            std::string bbTag = values[i];
            StringFormat::Strip(bbTag);
            std::cout << "Backbone will be modelled using " << bbTag << "\n";
            m_backboneTag = bbTag;
            }


            else if (keys[i] == "np-type") {
                std::cout << "npType " << values[i] << "\n";
                int trialVal =  AsInt(values[i]);
                if(trialVal == 1 || trialVal == 2 || trialVal == 3 || trialVal==4 || trialVal==5){
                m_npType = trialVal;
                                std::cout << "found NP type " << trialVal << "\n";
                }
                else{
                std::cout << "Unknown nanoparticle type, defaulting to sphere \n";
                m_npType = 1;
                }
            }
            else {
                std::cerr << "Error: Unknown parameter '" << keys[i] << "'\n";
                std::exit(1);
            }
        }
    }

    std::unordered_map<std::string, std::size_t> AminoAcidIdMap() const {
        std::unordered_map<std::string, std::size_t> aminoAcidIdMap;
        for (std::size_t id = 0; id < m_aminoAcids.size(); ++id) {
            aminoAcidIdMap.insert(std::pair<std::string, std::size_t>(m_aminoAcids[id], id));
        }
        return aminoAcidIdMap;
    }

    double AminoAcidRadius(const std::string& aminoAcid) const {
        if (m_aminoAcids.size() != m_aminoAcidRadii.size()) {
            std::clog << "Error: N != NR (" << m_aminoAcids.size() << ", " << m_aminoAcidRadii.size() << ")\n";
            std::exit(1);
        }
        for (std::size_t i = 0; i < m_aminoAcids.size(); ++i) {
            if (aminoAcid == m_aminoAcids[i]) {
                return m_aminoAcidRadii[i];
            }
        }
        std::cerr << "Error: Could not find a radius for the amino acid '" << aminoAcid << "'\n";
        std::exit(1);
    }

    double AminoAcidCharge(const std::string& aminoAcid) const {
        for (std::size_t i = 0; i < m_aminoAcids.size(); ++i) {
            if (aminoAcid == m_aminoAcids[i]) {
                return m_aminoAcidCharges[i];
            }
        }
        std::cerr << "Error: Could not find a charge for the amino acid '" << aminoAcid << "'\n";
        std::exit(1);
    }

private:
    static int AsInt(const std::string& str) {
        try {
            return std::stoi(str);
        }
        catch (const std::invalid_argument& ia) {
	        std::cerr << "Error: Invalid argument: " << ia.what() << " (" << str << ")\n";
            std::exit(1);
        }
    }

    static double AsDouble(const std::string& str) {
        try {
            return std::stod(str);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Error: Invalid argument: " << ia.what() << " (" << str << ")\n";
            std::exit(1);
        }
    }

    static std::vector<std::string> AsStringList(const std::string& str) {
        std::vector<std::string> stringList;
        if (str.front() == '[' && str.back() == ']') {
            StringFormat::Split(stringList, str.substr(1, str.size() - 2), ",");
        }
        else {
            stringList.emplace_back(str);
        }
        return stringList;
    }

    static std::vector<double> AsDoubleList(const std::string& str) {
        std::vector<std::string> stringList = AsStringList(str);
        std::vector<double> doubleList(stringList.size());
        for (std::size_t i = 0; i < stringList.size(); ++i) {
            doubleList[i] = AsDouble(stringList[i]);
        }
        return doubleList;
    }
};

#endif
