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

constexpr int potential_size = 1024;

class Potential {
public:
    double m_energy[potential_size];
    const double m_start;
    const double m_cutoff;
    const double m_dr;

public:
    Potential(const std::vector<double>& energy, const double start, const double cutoff)
        : m_start(start), m_cutoff(cutoff), m_dr((cutoff - start)/ (potential_size - 1.0))
    {
        for (int i = 0; i < (int)energy.size(); ++i) {
            m_energy[i] = energy[i];
        }
    }

    inline double Value(const double r) const noexcept {
        if (r > m_cutoff) {
            return 0.0;
        }
        else if (r < m_start) {
            return m_energy[0];
        }
        const int       index  = static_cast<int>((r - m_start) * (potential_size - 1.0) / (m_cutoff - m_start)); 
        const double    factor = (r - m_start) / m_dr - index;
        return m_energy[index + 1] * factor - m_energy[index] * (factor - 1.0);
    }

    inline double VariableSizeValue(const double r, const int size) const noexcept {
        if (r > m_cutoff) {
            return 0.0;
        }
        const int       index  = static_cast<int>((r - m_start) * (size - 1.0) / (m_cutoff - m_start)); 
        const double    factor = r / m_dr - static_cast<int>((r - m_start)  / m_dr);
        return m_energy[index + 1] * factor - m_energy[index] * (factor - 1.0);
    }

};

Potential GeneratePotential(const SurfaceData&, const HamakerConstants&, const double, const double, const Config&);
double HamakerPotential(const double, const double, const double, const double);
double ElectrostaticPotential(const double, const double, const double, const double, const double);
double HamakerPotentialV2(const double, const double, const double, const double, const double);

class Potentials : public std::vector<Potential> {
public:
    Potentials(const SurfacePMFs& surfacePMFs, const HamakerConstants& hamakerConstants, const double zetaPotential, const double nanoparticleRadius, const Config& config) {
        this->reserve(surfacePMFs.size());
        for (const auto& surfaceData : surfacePMFs) {
            this->emplace_back(GeneratePotential(surfaceData, hamakerConstants, zetaPotential, nanoparticleRadius, config));
        }
    }
};

Potential GeneratePotential(const SurfaceData& surfaceData, const HamakerConstants& hamakerConstants, const double zetaPotential, const double nanoparticleRadius, const Config& config) {

    const double        pmfStart            = surfaceData.m_distance.front();
    const double        pmfCutoff           = surfaceData.m_distance.back();
    const double        cutoff              = config.m_potentialCutoff;
    const double        hamaker             = hamakerConstants[surfaceData.m_aminoAcid];
    //const double        bjerumLength        = config.m_bejerumLength;
    const double        debyeLength         = config.m_debyeLength;
    const auto          aminoAcidRadius     = config.AminoAcidRadius(surfaceData.m_aminoAcid);
    const double        Z                   = config.AminoAcidCharge(surfaceData.m_aminoAcid);
    
    std::vector<double> energy(potential_size);

    SurfacePotential surface_potential(surfaceData.m_energy, surfaceData.m_distance);
    
    double  surface = 0.0, core = 0.0, electrostatic = 0.0;

    const std::string destination = "pot-dat"; 
    const std::string filename = destination + "/" + surfaceData.m_aminoAcid + ".dat";
    std::ofstream handle(filename.c_str());

    for (int i = 0; i < potential_size; ++i) {
        const double r = pmfStart + (i / (potential_size - 1.0)) * (cutoff - pmfStart);
        double U       = 0.0;

        if (config.m_enableSurface) {
            surface = surface_potential.Value(r, nanoparticleRadius, pmfCutoff,config.m_npType);
            U += surface;
        }

        if (config.m_enableCore) {
            if(config.m_npType == 1){ //sphere
                core = HamakerPotentialV2(hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff);
            }
            else if(config.m_npType == 2){ //cylinder, defined in CylinderPotential.h
                core = HamakerSphereCylinder(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff);
            }
            else if(config.m_npType == 3){//cube, defined in CubePotential.h
                core =  HamakerSphereCube(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff);
            }
            else if(config.m_npType == 4){ //tube, defined in TubePotential.h
                core =  HamakerSphereTube(hamaker,aminoAcidRadius,nanoparticleRadius,r,pmfCutoff);
            }
            else{
                core = HamakerPotentialV2(hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff);
            }
            U += core;
        }

        if (config.m_enableElectrostatic) {
            electrostatic = ElectrostaticPotential(r, zetaPotential, Z, nanoparticleRadius, debyeLength);
            U += electrostatic;
        }

        energy[i] = U;
    
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << r;
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << U;
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << surface;
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << core;
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << electrostatic;
        handle << "\n";
    }

    handle.close(); 

    return Potential(energy, pmfStart, cutoff);
}

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

double ElectrostaticPotential(const double h, const double sigma, const double Z, const double R, const double ld) {
    return 38.681727 * (sigma * Z) * (R/ (R + h)) * exp(-1.0 * h / ld); // The constant is e/kT
}

double HamakerPotentialV2(const double A, const double R1, const double R2, const double r, const double cutoff) {

  const double C  = r + R2; // Center to center distance
double re = cutoff;
double d = C - R1 - R2;
long double res = 0;
//If re > min(r1,r2) then the potential should remain finite for all values of C, as there is no divergence. In practice the expression found here diverges for C < R2 even if this condition holds - this may be due to a geometrical issue with how the integral is performed. 
if( C < std::min(R1,R2)){
std::cout << "warning: overlap is too large! \n";
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
