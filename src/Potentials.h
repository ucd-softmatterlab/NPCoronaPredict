#ifndef POTENTIALS__H__
#define POTENTIALS__H__

#include "Surface.h"
#include "HamakerConstants.h"
#include "Config.h"
#include "CylinderPotential.h"
#include "CubePotential.h"
#include "TubePotential.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <boost/math/special_functions/bessel.hpp> 
#include "CubeESPotential.h"

constexpr int potential_size = 1024;
constexpr int max_np_types= 64;

class Potential {
public:
    double m_energy[max_np_types][potential_size];
    const double m_start;
    const double m_cutoff;
    const double m_dr;

public:
    Potential(const std::vector<std::vector<double>>& energy, const double start, const double cutoff)
        : m_start(start), m_cutoff(cutoff), m_dr((cutoff - start)/ (potential_size - 1.0))
    {
        for (int i = 0; i < potential_size; ++i) {
           for (int j = 0; j < max_np_types; ++j) {
            m_energy[j][i] = energy[j][i];
        }
        }
    }

    inline double Value(const double r, const int npType) const noexcept {
        if (r > m_cutoff) {
            return 0.0;
        }
        else if (r < m_start) {
            return 200.0 * (m_start - r)*(m_start-r) + m_energy[npType][0];
        }
        const int       index  = static_cast<int>((r - m_start) * (potential_size - 1.0) / (m_cutoff - m_start)); 
        const double    factor = (r - m_start) / m_dr - index;
        return m_energy[npType][index + 1] * factor - m_energy[npType][index] * (factor - 1.0);
    }

    inline double VariableSizeValue(const double r, const int size, const int npType) const noexcept {
        if (r > m_cutoff) {
            return 0.0;
        }
        const int       index  = static_cast<int>((r - m_start) * (size - 1.0) / (m_cutoff - m_start)); 
        const double    factor = r / m_dr - static_cast<int>((r - m_start)  / m_dr);
        return m_energy[npType][index + 1] * factor - m_energy[npType][index] * (factor - 1.0);
    }

};

Potential GeneratePotential(const SurfaceData&, const HamakerConstants&, const double, const double, const Config&, const NP&);
double HamakerPotential(const double, const double, const double, const double);
double ElectrostaticPotential(const double, const double, const double, const double, const double);
double ElectrostaticCylinderPotential(const double, const double, const double, const double, const double);
double ElectrostaticCubePotential(const double, const double, const double, const double, const double);
double HamakerPotentialV2(const double, const double, const double, const double, const double);
double SmallHamaker(const double, const double, const double, const double, const double);

class Potentials : public std::vector<Potential> {
public:
    Potentials(const SurfacePMFs& surfacePMFs, const HamakerConstants& hamakerConstants, const double zetaPotential, const double nanoparticleBoundingRadius, const Config& config, const NP& npComponents) {
        this->reserve(surfacePMFs.size());
        for (const auto& surfaceData : surfacePMFs) {
            this->emplace_back(GeneratePotential(surfaceData, hamakerConstants, zetaPotential, nanoparticleBoundingRadius, config, npComponents));
        }
    }
};

//Generate the set of potential as a vector of vectors such that potential[npIndex][distanceIndex] returns the potential associated with this amino acid for each NP component at a 
Potential GeneratePotential(const SurfaceData& surfaceData, const HamakerConstants& hamakerConstants, const double zetaPotentialAtBound, const double nanoparticleBoundingRadius, const Config& config, const NP& npComponents) {

//    const double        pmfStart            = surfaceData.m_distance.front();
//    const double        pmfEnd              = surfaceData.m_distance.back();  //final distance recorded in the PMF
//    const double        pmfCutoff           = config.m_PMFCutoff;  //vdW cutoff
    const double        cutoff              = config.m_potentialCutoff; //final distance of the total potential
//    const double        hamaker             = hamakerConstants[surfaceData.m_aminoAcid];
    const double        bjerumLength        = config.m_bejerumLength;
    const double        debyeLength         = config.m_debyeLength;
//    const auto          aminoAcidRadius     = config.AminoAcidRadius(surfaceData.m_aminoAcid);
//    const double        Z                   = config.AminoAcidCharge(surfaceData.m_aminoAcid);
    const int           recalculateZetaPotential =  config.m_recalcZP;
    
    std::vector<std::vector<double>> energy( max_np_types, std::vector<double>(potential_size));

    double potentialStart = config.m_potentialStart;
 

    const std::string destination = config.m_outputDirectory +"/pot-dat/" + npComponents.m_name  ; 
    //const std::string filename = destination + "/" + surfaceData.m_aminoAcid + ".dat";

    boost::filesystem::create_directory("pot-dat") ;
    boost::filesystem::create_directories(destination) ;



    //std::cout << surfaceData.m_aminoAcid << " start: " << potentialStart << "\n";
    
    //if we're pre-summing over all NPs, j needs to run over all the NP beads
    //if not, it needs to run only over types of NP beads
    
    int jLimit = npComponents.m_numBeadTypes;
    if( config.m_sumNPPotentials == true){
    jLimit = npComponents.m_x.size();
    }
    
    
   for(int j = 0; j < jLimit; ++j){
   std::vector<string> aaName;
    aaName.push_back(surfaceData.m_aminoAcid); 
    //std::cout<< "Component: " << j << "\n";
    
    
    int beadFileIndex = j;
      if( config.m_sumNPPotentials == true){
      beadFileIndex = npComponents.m_npBeadType[j] ;
      
      }
    
    
   // std::cout << "Loading Hamaker file: " << npComponents.m_hamakerFile[beadFileIndex] << "\n";
    HamakerConstants  hamakerConstantsComponent(npComponents.m_hamakerFile[beadFileIndex]);

    //The SurfacePMFs object is responsible for reading in a set of PMF files, here set to a single file generated from aaName.
    SurfacePMFs surfaceDataComponentSet(npComponents.m_pmfFile[beadFileIndex], config.m_pmfPrefix, aaName);
    SurfaceData surfacePMFComponent = surfaceDataComponentSet.front();
    SurfacePotential surface_potential(surfacePMFComponent.m_energy, surfacePMFComponent.m_distance);


    double hamaker = hamakerConstantsComponent[surfaceData.m_aminoAcid];
    auto          aminoAcidRadius     = config.AminoAcidRadius(surfaceData.m_aminoAcid);
    double        Z                   = config.AminoAcidCharge(surfaceData.m_aminoAcid);




    double  surface = 0.0, core = 0.0, electrostatic = 0.0;
//    double potentialStart = 0.4; //config.m_potentialStart (when added)
   // const std::string destination = "pot-dat/" + npComponents.m_name  ; 
   // const std::string filename = destination + "/" + surfaceData.m_aminoAcid + ".dat";
  //  std::ofstream handle(filename.c_str());
  double component_centre[3] = {0,0,0};
  int npShape ;
double nanoparticleRadius;

double component_corePrefactor ;
double component_surfacePrefactor ;
double zetaPotential ;
double pmfCutoff ;

int correctionTypeOverride ;
  
  
  
  if( config.m_sumNPPotentials == true){
 // std::cout << "Summing NPs: j = NP index \n";
  //running over all NPs, j = NP index
   component_centre[0] = npComponents.m_x[j];
     component_centre[1] = npComponents.m_y[j];
        component_centre[2] = npComponents.m_z[j];
  int npBeadType = npComponents.m_npBeadType[j];
  
   npShape = npComponents.m_shape[npBeadType];
   nanoparticleRadius = npComponents.m_radius[npBeadType];

 component_corePrefactor = npComponents.m_coreFactor[npBeadType];
 component_surfacePrefactor  = npComponents.m_surfFactor[npBeadType];
 zetaPotential = npComponents.m_zeta[npBeadType];
 pmfCutoff = npComponents.m_pmfCutoff[npBeadType];

  correctionTypeOverride = npComponents.m_correctionType[npBeadType];
   
   
  }
  else{
  //running over NP types, j = NP-type index

 npShape = npComponents.m_shape[j];
  nanoparticleRadius = npComponents.m_radius[j];
  component_corePrefactor = npComponents.m_coreFactor[j];
  component_surfacePrefactor  = npComponents.m_surfFactor[j];
  zetaPotential = npComponents.m_zeta[j];
  pmfCutoff = npComponents.m_pmfCutoff[j];

 correctionTypeOverride = npComponents.m_correctionType[j];
}



double temperatureFactor = 300.0/config.m_temperature;
 //  std::cout<< "NP: " << j << " radius " << nanoparticleRadius << " core-material: " << component_corePrefactor << "*"<< npComponents.m_hamakerFile[j] << " val: " << hamaker*component_corePrefactor   << ", surface-material: " << component_surfacePrefactor << "*"  << npComponents.m_pmfFile[j] << " sample at distance 0.5 " << surface_potential.Value(0.5,10000,1) << "\n";

    for (int i = 0; i < potential_size; ++i) {
        const double r = potentialStart + (i / (potential_size - 1.0)) * (cutoff - potentialStart);
        double U       = 0.0;



//r is the distance from the surface of the NP to the centre of the bead.
//if NP-summation is disabled, then potentials are saved individually for each NP and r is the distance from each component surface to the AA COM
//if NP-summation is enabled, then r is the distance from the nominal NP surface
 
       double rstar  = r;
       if(config.m_sumNPPotentials == true){
       
        //cylindrical distance
       if(npShape == 2 || npShape == 4 || npShape == 5){
rstar = sqrt(   pow(component_centre[1],2) + pow(nanoparticleBoundingRadius + r - component_centre[2],2) ) - nanoparticleRadius;
       }
      else if(npShape == 3){
//for now we approximate the effective distance from a cube by that of the enscribed sphere
 rstar = sqrt( pow(component_centre[0],2) + pow(component_centre[1],2) + pow(nanoparticleBoundingRadius + r - component_centre[2],2) ) - nanoparticleRadius;
     }
else{
//default to spherical
           //the centre of the NP bead is at component_centre[0],cc[1],cc[2]
            //the AA bead is located at npBoundingRadius + r
        rstar = sqrt( pow(component_centre[0],2) + pow(component_centre[1],2) + pow(nanoparticleBoundingRadius + r - component_centre[2],2) ) - nanoparticleRadius; //the actual NP-surface to bead-centre distance for the NP component;
}
       
       
       
       
       }


        bool addShortRangeRepulsion = true;
        if(addShortRangeRepulsion == true){
        if(rstar < 0.1){
        double wcaHS = pow(0.1/rstar, 12.0) - 1.0 ;
        U += wcaHS;
        }
        
        }


        if (config.m_enableSurface) {
            surface = surface_potential.Value(rstar, nanoparticleRadius, pmfCutoff,correctionTypeOverride);
            //std::cout << "bead type "<< j << " distance " << rstar << " surface " << surface << " scale: " << component_surfacePrefactor << "\n";
            U += component_surfacePrefactor*surface;
        }

        if (config.m_enableCore ) {
            if( nanoparticleRadius > pmfCutoff){
            if(npShape == 1){ //sphere
                core = HamakerPotentialV2(hamaker, aminoAcidRadius, nanoparticleRadius, rstar, pmfCutoff);
            }
            else if(npShape == 2 || npShape == 5){ //cylinder, defined in CylinderPotential.h, applicable to solid cylinders and MWCNT
                core = HamakerSphereCylinder(hamaker,aminoAcidRadius,nanoparticleRadius,rstar,pmfCutoff);
             }
            else if(npShape == 3){//cube, defined in CubePotential.h
                core =  HamakerSphereCube(hamaker,aminoAcidRadius,nanoparticleRadius,rstar,pmfCutoff);
            }
            else if(npShape == 4){ //tube, defined in TubePotential.h
                //core =  HamakerSphereTube(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff);
               //approximate the tube as the difference of the outer cylinder and an inner cylinder of radius R - 0.34, where 0.34 is the approx. thickness of one "layer"
              //since r is the NP surface - AA distance we also need to adjust this value
                double wallThickness = 0.34;
                double innerR = nanoparticleRadius - wallThickness;
             
/*
  if(innerR>0){
                core = HamakerSphereCylinder(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff) - HamakerSphereCylinder(hamaker,aminoAcidRadius, innerR,r+wallThickness,pmfCutoff);
}            
else{
  core = HamakerSphereCylinder(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff);
}
*/

core =wallThickness *  HamakerSphereTube(hamaker,aminoAcidRadius,nanoparticleRadius,rstar,pmfCutoff);


}
            else{
                core = HamakerPotentialV2(hamaker, aminoAcidRadius, nanoparticleRadius, rstar, pmfCutoff);
            }

             
 if(core*core > 1e6){
std::cout << "large core at rstar " << rstar << " " << core <<  " from NP component " << j  <<  "\n";
std::cout << "NP surface - AA centre distance: " << rstar << " bounding radius " << nanoparticleBoundingRadius << " NP component radius " << nanoparticleRadius << "AA radius" << aminoAcidRadius << "\n";
}
            }
            else{
                core = SmallHamaker(hamaker, aminoAcidRadius,nanoparticleRadius,rstar,pmfCutoff); // handle this more carefully
            }
            U +=  component_corePrefactor*core;
            
            
        }

        if (config.m_enableElectrostatic) {

        double finalZetaPotential = zetaPotential;
      



            if(npShape == 2 || npShape == 4 || npShape==5){

if( recalculateZetaPotential != 0){
           finalZetaPotential = 2* zetaPotential*bjerumLength*debyeLength * boost::math::cyl_bessel_k(0, nanoparticleRadius/debyeLength)/boost::math::cyl_bessel_k(1,  nanoparticleRadius/debyeLength); //cylinder correction

  // std::cout << " reference zeta potential: " << zetaPotential << " recalculated: " << finalZetaPotential << "\n";
        }


            electrostatic = ElectrostaticCylinderPotential(rstar,finalZetaPotential, Z, nanoparticleRadius, debyeLength);
}
else if(npShape == 3){

if( recalculateZetaPotential  != 0){
           finalZetaPotential = 2*zetaPotential*bjerumLength*debyeLength * 1/( 1 + 0.59692 /(nanoparticleRadius/debyeLength) + 0.00626772/pow(nanoparticleRadius/debyeLength,2)    ); //cube correction
 //  std::cout << " reference zeta potential: " << zetaPotential << " recalculated: " << finalZetaPotential << "\n";
        }


electrostatic = ElectrostaticCubePotential(rstar, finalZetaPotential, Z, nanoparticleRadius, debyeLength);
}
else{

if( recalculateZetaPotential  != 0){
           finalZetaPotential = 2*zetaPotential   *bjerumLength*debyeLength *(nanoparticleRadius/debyeLength)/(1+ nanoparticleRadius/debyeLength); //sphere correction
//   std::cout << " reference zeta potential: " << zetaPotential << " recalculated: " << finalZetaPotential << "\n";
        }

            electrostatic = ElectrostaticPotential(rstar, finalZetaPotential, Z, nanoparticleRadius, debyeLength);
}
            U += electrostatic;
        }

        //If we're explicitly separating NP potentials, store these to individual components. If not save everything to component 0
        if(config.m_sumNPPotentials == false){
        energy[j][i]  = U*temperatureFactor;
         }
         else{
         energy[0][i] += U*temperatureFactor;
         }
         
        //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << r;
       //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << U;
       // handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << surface;
       // handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << core;
        //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << electrostatic;
        //handle << "\n";
    }
//std::cout << " j added, energy at test: " << energy[0] << "\n";

}
    //const std::string destination = "pot-dat";
    const std::string filename = destination + "/" + surfaceData.m_aminoAcid + ".dat";
    //boost::filesystem::create_directory("pot-dat") ;
    //boost::filesystem::create_directory(destination) ;



    std::ofstream handle(filename.c_str());
    for (int i = 0; i < potential_size; ++i) {
        const double r = potentialStart + (i / (potential_size - 1.0)) * (cutoff - potentialStart);

        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << r;
        for(int j = 0; j < npComponents.m_numBeadTypes; ++j){
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << energy[j][i];
        }
        //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << surface;
        //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << core;
        //handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << electrostatic;
        handle << "\n";
}
    handle.close(); 

    return Potential(energy, potentialStart, cutoff);
}
//A = hamaker constant, R1 = AA radius, R2 = NP radius, r = npsurf-aacentre distance, cutoff = r_exclusion
double HamakerPotential(const double A, const double R1, const double R2, const double r, const double cutoff) {

  const double C  = r + R2; // Center to center distance

  if (r > cutoff) { // Classical Hamaker Potential
    return -1.0 * (A / 6.0) * (
      (2.0 * R1 * R2) / (C * C - (R1 + R2) * (R1 + R2)) +
      (2.0 * R1 * R2) / (C * C - (R1 - R2) * (R1 - R2)) +
      std::log((C * C - (R1 + R2) * (R1 + R2)) / (C * C - (R1 - R2) * (R1 - R2)))
    );
  }
  else { // Modified Hamaker Potential (with missing lens piece)
    double top = 2.0 * R1 * (C + R2 - cutoff) * (-1.0 * C * C * C * R1 * R1 + R1 * R1 * (R2 - cutoff) * std::pow(R2 + cutoff, 2.0)
      + C * C * (-1.0 * R1 * R1 * R2 + std::pow(cutoff, 3.0)) + C * (R1 * R1 * R1 * R1 + R2 * std::pow(cutoff, 3.0) + R1 * R1 * R2 * (R2 + cutoff)));
    double bottom = C * (C - R1 + R2) * (C + R1 + R2) * std::pow(R1 - cutoff, 2.0) * std::pow(R1 + cutoff, 2.0);
    double logterm = -1.0 * std::log(C - R1 + R2) + std::log(C + R1 + R2) + std::log(cutoff - R1) - 1.0 * std::log(R1 + cutoff);
    return -1.0 * (A / 6.0) * ((top / bottom) + logterm);
  }
}

double ElectrostaticPotential(const double h, const double zetaS, const double Z, const double R, const double ld) {
    return 38.681727 * (zetaS * Z) * (R/ (R + h)) * exp(-1.0 * h / ld); // The constant is e/kT to handle the unit conversion
}
 

//Defines the approximate debye-huckel potential for an infinitely long cylinder in a solution with a fixed boundary condition at the surface and psi->0 for large distances
//The normalisation convention here is chosen such that at h = 0 (i.e. at the surface of the NP) this should return the same value as the spherical case.
double ElectrostaticCylinderPotential(const double h, const double zetaS, const double Z, const double R, const double ld) {
    double cylinderES = 38.681727 * (zetaS * Z)   * boost::math::cyl_bessel_k(0,  (h+R)/ld ) / boost::math::cyl_bessel_k(0,  R/ld);

    return cylinderES;
}


//calculates the debye-huckel potential for a planar surface in a weakly ionic solution.
//this is a reasonably accurate model for a cubic nanoparticle as long as the side of a cube is a few times larger than the Debye length.
//if it's smaller than this then we use an approximation which is valid for cubes for which R/LD < 10 or so but still requires that the cube is reasonably large as this only returns the potential for the (x=0,y=0,z) line.
double ElectrostaticCubePotential(const double h, const double zetaS, const double Z, const double R, const double ld) {
double cubeES=0;
    if(R/ld > 10){
    cubeES = 38.681727 * (zetaS * Z)   * exp(-h/ld);
    }
    else{
    cubeES = 38.681727 * (zetaS * Z)   * CubeElectrostaticOrder8(R,1/ld, h );
    }
    return cubeES;
}


double SmallHamaker(const double A, const double R1, const double R2, const double r, const double cutoff) {
     const double C  = r + R2; // Center to center distance
double re = cutoff;
double d = C - R1 - R2;
long double res = 0;

//case 1: beads are far apart, return standard hamaker
if(re < 1e-20 || d > re){
res =  -A * 1/6.0 *( 2*R1*R2/(C*C -(R1+R2)*(R1+R2)) + 2*R1*R2/ (C*C -(R1-R2)*(R1-R2)) + std::log(( C*C - (R1+R2)*(R1+R2) )/( C*C - (R1-R2)*(R1-R2) ))   );
}
else if( C   + R1 + R2 < re  ){ // case 2: all of bead 1 and 2 are in the exclusion region, return 0 
res = 0;
    
}
  else{
      //at least some of one sphere is within the exclusion zone
      res = 0;
      
  }
  
  
  
  return res;
}

double HamakerPotentialV2(const double A, const double R1, const double R2, const double r, const double cutoff) {

  double C  = r + R2; // Center to center distance
double re = cutoff;
double d = C - R1 - R2;
long double res = 0;
//If re > min(r1,r2) then the potential should remain finite for all values of C, as there is no divergence. In practice the expression found here diverges for C < R2 even if this condition holds - this may be due to a geometrical issue with how the integral is performed. 
if( C < std::min(R1,R2)){
//std::cout << "warning: overlap is too large! \n";
//a safe C value is then min(R1,R2) + epsilon
//C = r + R2 , so the safe value of r to use going forwards is C- R2

//double hamakerNewCall = HamakerPotentialV2( A, R1, R2, r + R2+0.05, cutoff) ;
//return hamakerNewCall;
C = std::min(R1,R2)+0.01;
d = C - R1-R2;
} 
 
int stepFunc1 = (d > re) ? 1:0;
int stepFunc2 = (d + 2*R1 > re) ? 1:0;
int stepFunc3 = (d + 2*R2 > re) ? 1:0;
int stepFunc4 = (d+2*R1+2*R2 > re) ? 1:0;
//double hamakerA = 3.1415*3.1415; //pi^2 q_1 q_2 lambda

//if d can be negative then we need to be more careful due to the fact it causes some logs to return NaN, and although these are multiplied by 0 from the step functions, 0*NaN is still NaN and this throws errors later on. so we catch these here and manually set them to 0 instead.
//these are superseded when using the if else else else else else else else statement as the cases there ensure the log arguments are non-zeor. 
double logdre = (d > re) ? log(d/re):0.0;
double logdrer1 = (d + 2*R1 > re) ? log((d + 2*R1)/re):0.0;
//double logdrer2 = (d + 2*R2 > re) ? log((d + 2*R2)/re):0.0;
//double logdrer1r2 = (d+2*R1+2*R2 > re) ? log(re/(d + 2*(R1 + R2))):0.0;



int useLarge = 0;
//if re is less than a very small amount or if d > re then we compute the standard Hamaker form as this is quicker
if(re < 1e-20 || d > re){
res =  -A * 1/6.0 *( 2*R1*R2/(C*C -(R1+R2)*(R1+R2)) + 2*R1*R2/ (C*C -(R1-R2)*(R1-R2)) + std::log(( C*C - (R1+R2)*(R1+R2) )/( C*C - (R1-R2)*(R1-R2) ))   );
}
else if(useLarge == 1){
std::cout << "using large NP potential \n";
res =  A * (16*std::pow(R1,3)*(3*d + 3*R1 - 4*re) - (stepFunc1*((d - re)*(3*std::pow(d,4) + std::pow(d,3)*(12*R1 - 13*re) + d*(36*R1 - 25*re)*std::pow(re,2) - 12*R1*std::pow(re,3) + std::pow(d,2)*re*(-36*R1 + 23*re)) + 12*d*std::pow(re,4)*logdre))/d + (stepFunc2*(3*std::pow(d,5) + 2*std::pow(d,4)*(9*R1 - 8*re) + 4*std::pow(d,3)*(6*std::pow(R1,2) - 20*R1*re + 9*std::pow(re,2)) -  48*std::pow(d,2)*(std::pow(R1,3) + 2*std::pow(R1,2)*re - 3*R1*std::pow(re,2) + std::pow(re,3)) + 2*R1*(-48*std::pow(R1,4) + 64*std::pow(R1,3)*re - 48*R1*std::pow(re,3) + 19*std::pow(re,4)) +     d*(-144*std::pow(R1,4) + 64*std::pow(R1,3)*re + 144*std::pow(R1,2)*std::pow(re,2) - 144*R1*std::pow(re,3) + 25*std::pow(re,4)) + 12*(d + 2*R1)*std::pow(re,4)*logdrer1))/(d + 2*R1))/(72.*std::pow(re,4));

}
else{
//if not then we need to use a more complex expression. to reduce the numerical error we can still subdivide this into slightly less awful expressions based on which of the step functions are valid.
/*
res =  -(A*stepFunc1*(3*pow(d,4)*(pow(d,2) + 20*R1*R2 + 5*d*(R1 + R2)) - 20*pow(d,3)*(pow(d,2) + 12*R1*R2 + 4*d*(R1 + R2))*re + 60*pow(d,2)*(pow(d,2) + 6*R1*R2 + 3*d*(R1 + R2))*pow(re,2) - 120*d*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) + 5*(13*pow(d,2) + 12*R1*R2 + 25*d*(R1 + R2))*pow(re,4) +  12*d*pow(re,5) + 60*d*(d + R1 + R2)*pow(re,4)*logdre))/(360.*d*(d + R1 + R2)*pow(re,4)) + (A*stepFunc2*(3*pow(d + 2*R1,4)*((d - 3*R1)*(d + 2*R1) + 5*(d - 2*R1)*R2) - 20*pow(d + 2*R1,3)*(pow(d,2) + 4*d*R2 - 4*R1*(R1 + R2))*re +   60*pow(d + 2*R1,2)*(pow(d,2) - 2*pow(R1,2) + d*(R1 + 3*R2))*pow(re,2) - 120*(d + 2*R1)*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) +   5*(13*pow(d,2) + 2*R1*(R1 + 19*R2) + d*(27*R1 + 25*R2))*pow(re,4) + 12*(d + 2*R1)*pow(re,5) + 60*(d + 2*R1)*(d + R1 + R2)*pow(re,4)*logdrer1))/(360.*(d + 2*R1)*(d + R1 + R2)*pow(re,4)) + (A*stepFunc3* (3*pow(d + 2*R2,4)*(d*(d + 5*R1) - (d + 10*R1)*R2 - 6*pow(R2,2)) - 20*pow(d + 2*R2,3)*(pow(d,2) + 4*d*R1 - 4*R2*(R1 + R2))*re +   60*pow(d + 2*R2,2)*(pow(d,2) - 2*pow(R2,2) + d*(3*R1 + R2))*pow(re,2) - 120*(d + 2*R2)*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) +   5*(13*pow(d,2) + 2*R2*(19*R1 + R2) + d*(25*R1 + 27*R2))*pow(re,4) + 12*(d + 2*R2)*pow(re,5) + 60*(d + R1 + R2)*(d + 2*R2)*pow(re,4)*logdrer2))/ (360.*(d + R1 + R2)*(d + 2*R2)*pow(re,4)) + (A*stepFunc4* (-3*pow(d,6) - 21*pow(d,5)*R1 - 30*pow(d,4)*pow(R1,2) + 120*pow(d,3)*pow(R1,3) + 480*pow(d,2)*pow(R1,4) + 624*d*pow(R1,5) + 288*pow(R1,6) - 21*pow(d,5)*R2 - 120*pow(d,4)*R1*R2 -  120*pow(d,3)*pow(R1,2)*R2 + 480*pow(d,2)*pow(R1,3)*R2 + 1200*d*pow(R1,4)*R2 + 768*pow(R1,5)*R2 - 30*pow(d,4)*pow(R2,2) - 120*pow(d,3)*R1*pow(R2,2) + 480*d*pow(R1,3)*pow(R2,2) +  480*pow(R1,4)*pow(R2,2) + 120*pow(d,3)*pow(R2,3) + 480*pow(d,2)*R1*pow(R2,3) + 480*d*pow(R1,2)*pow(R2,3) + 480*pow(d,2)*pow(R2,4) + 1200*d*R1*pow(R2,4) +   480*pow(R1,2)*pow(R2,4) + 624*d*pow(R2,5) + 768*R1*pow(R2,5) + 288*pow(R2,6) + 20*pow(d + 2*(R1 + R2),3)*(pow(d,2) - 4*(pow(R1,2) - R1*R2 + pow(R2,2)))*re -       60*pow(d + 2*(R1 + R2),2)*(pow(d,2) + d*(R1 + R2) - 2*(pow(R1,2) - R1*R2 + pow(R2,2)))*pow(re,2) + 120*(d + 2*(R1 + R2))*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) -   5*(13*pow(d,2) + 27*d*(R1 + R2) + 2*(pow(R1,2) + 8*R1*R2 + pow(R2,2)))*pow(re,4) - 12*(d + 2*(R1 + R2))*pow(re,5) + 60*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4)*logdrer1r2))/   (360.*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4)) ;

 */

//begin the if else else else else else statement to avoid looking up the syntax for cases
if(stepFunc1&&stepFunc2&&stepFunc3&&stepFunc4){
res = (A*((-4*R1*R2*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2)))/(d*(d + 2*R1)*(d + 2*R2)*(d + 2*(R1 + R2))) + log(((d + 2*R1)*(d + 2*R2))/(d*(d + 2*(R1 + R2))))))/6.; //Case A. This is just the standard Hamaker expression. 
}
else if(stepFunc2&&stepFunc3&&stepFunc4){ //case B
res = (A*(3*pow(d,3)*(pow(d,2) + 20*R1*R2 + 5*d*(R1 + R2)) - 20*pow(d,2)*(pow(d,2) + 12*R1*R2 + 4*d*(R1 + R2))*re + 60*d*(pow(d,2) + 6*R1*R2 + 3*d*(R1 + R2))*pow(re,2) -   120*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) + (5*(d*pow(d + 2*R1,2)*(13*d + 25*R1) + (d + 2*R1)*(77*pow(d,2) + 166*d*R1 + 76*pow(R1,2))*R2 + 8*(d + R1)*(19*d + 32*R1)*pow(R2,2) +  4*(25*d + 38*R1)*pow(R2,3))*pow(re,4))/((d + 2*R1)*(d + 2*R2)*(d + 2*(R1 + R2))) + 12*pow(re,5) + 60*(d + R1 + R2)*pow(re,4)*log(((d + 2*R1)*(d + 2*R2))/((d + 2*(R1 + R2))*re))))/  (360.*(d + R1 + R2)*pow(re,4));
}
else if(stepFunc2&&stepFunc4 &&(!stepFunc3)){ //case C
res = (A*(6*(d + 2*R1)*pow(R2,3)*(5*d*(d + 2*R1) + 10*(d + R1)*R2 + 6*pow(R2,2))*(d + 2*(R1 + R2)) - 80*(d + 2*R1)*pow(R2,3)*(d + R1 + R2)*(d + 2*(R1 + R2))*re +  60*(d + 2*R1)*pow(R2,3)*(d + 2*(R1 + R2))*pow(re,2) + 30*R2*(pow(d,2) + 3*d*R1 + 2*pow(R1,2) + 2*d*R2 + 3*R1*R2)*pow(re,4) +   15*(d + 2*R1)*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4)*log((d + 2*R1)/(d + 2*(R1 + R2)))))/(90.*(d + 2*R1)*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4));
}
else if(stepFunc3&&stepFunc4 &&(!stepFunc2)){ //case D
res = (A*(6*pow(R1,3)*(d + 2*R2)*(d + 2*(R1 + R2))*(5*pow(d,2) + 10*d*(R1 + R2) + 2*R1*(3*R1 + 5*R2)) - 80*pow(R1,3)*(d + R1 + R2)*(d + 2*R2)*(d + 2*(R1 + R2))*re +  60*pow(R1,3)*(d + 2*R2)*(d + 2*(R1 + R2))*pow(re,2) + 30*R1*(d*(d + 2*R1) + 3*(d + R1)*R2 + 2*pow(R2,2))*pow(re,4) +  15*(d + R1 + R2)*(d + 2*R2)*(d + 2*(R1 + R2))*pow(re,4)*log((d + 2*R2)/(d + 2*(R1 + R2)))))/(90.*(d + R1 + R2)*(d + 2*R2)*(d + 2*(R1 + R2))*pow(re,4));
}
else if(stepFunc4){ // case E
res = (A*(-3*pow(d,6) - 21*pow(d,5)*R1 - 30*pow(d,4)*pow(R1,2) + 120*pow(d,3)*pow(R1,3) + 480*pow(d,2)*pow(R1,4) + 624*d*pow(R1,5) + 288*pow(R1,6) - 21*pow(d,5)*R2 -  120*pow(d,4)*R1*R2 - 120*pow(d,3)*pow(R1,2)*R2 + 480*pow(d,2)*pow(R1,3)*R2 + 1200*d*pow(R1,4)*R2 + 768*pow(R1,5)*R2 - 30*pow(d,4)*pow(R2,2) - 120*pow(d,3)*R1*pow(R2,2) +      480*d*pow(R1,3)*pow(R2,2) + 480*pow(R1,4)*pow(R2,2) + 120*pow(d,3)*pow(R2,3) + 480*pow(d,2)*R1*pow(R2,3) + 480*d*pow(R1,2)*pow(R2,3) + 480*pow(d,2)*pow(R2,4) +      1200*d*R1*pow(R2,4) + 480*pow(R1,2)*pow(R2,4) + 624*d*pow(R2,5) + 768*R1*pow(R2,5) + 288*pow(R2,6) + 20*pow(d + 2*(R1 + R2),3)*(pow(d,2) - 4*(pow(R1,2) - R1*R2 + pow(R2,2)))*re -   60*pow(d + 2*(R1 + R2),2)*(pow(d,2) + d*(R1 + R2) - 2*(pow(R1,2) - R1*R2 + pow(R2,2)))*pow(re,2) + 120*(d + 2*(R1 + R2))*(pow(d,2) + 2*R1*R2 + 2*d*(R1 + R2))*pow(re,3) -    5*(13*pow(d,2) + 27*d*(R1 + R2) + 2*(pow(R1,2) + 8*R1*R2 + pow(R2,2)))*pow(re,4) - 12*(d + 2*(R1 + R2))*pow(re,5) + 60*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4)*log(re/(d + 2*(R1 + R2)))))/  (360.*(d + R1 + R2)*(d + 2*(R1 + R2))*pow(re,4));
}
else{
res = 0;
}

}
 
return (double)res;


}

#endif
