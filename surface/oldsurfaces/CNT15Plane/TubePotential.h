#ifndef TUBEPOTENTIALS__H__
#define TUBEPOTENTIALS__H__

 

#include <cmath>
 

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

 

/*This file defines the functions required to calculate the potential due to a hollow cylinder (a tube) of radius RT and length 2L spanning the interval [-L ,L ] for an atom located at (R1,0,0). 
*/


//the Mathematica version of the elliptic integrals uses a different argument to the boost definition. So here we define wrapper functions to convert from the Mathematica to the Boost argument and also match the naming conventions.
 
 
double EllipticIncompleteE(double phi, double m){
return boost::math::ellint_2(  sqrt(m) , phi);
}
double EllipticIncompleteF(double phi, double m){
return boost::math::ellint_1(  sqrt(m), phi );
}

//here, r1 is the distance from the atom to the centre of the cylinder, RT is the radius of the cylinder and re is the exclusion region.
//the nanotube is infinitely tall because that's a thing that makes sense for a nanoparticle.
 

//this function returns the indefinite integral outside the exclusion region.
// this should be evaluated on the interval [zBound, L] where zBound = sqrt(re^2 - (ra-rc)^2) to get the potential for everything above zbound to the end of the tube, then doubled to take into account the other half of the tube 
double TubeExclusionFree(double RT, double r1,   double zc){
double res = (M_PI*RT*(4*sqrt(pow(r1 - RT,2))*(pow(r1,3) + pow(r1,2)*RT + r1*pow(RT,2) + pow(RT,3))*zc*pow(pow(r1,2) - 2*r1*RT + pow(RT,2) + pow(zc,2),2)*(pow(r1,2) + 2*r1*RT + pow(RT,2) + pow(zc,2)) -    sqrt(pow(pow(r1,2) - pow(RT,2),2))*(zc*(9*pow(r1,8) - 16*pow(r1,7)*RT + 9*pow(RT,8) + 20*pow(RT,6)*pow(zc,2) + 15*pow(RT,4)*pow(zc,4) + 4*pow(RT,2)*pow(zc,6) - 16*r1*pow(RT,3)*pow(pow(RT,2) + pow(zc,2),2) +   4*pow(r1,6)*(3*pow(RT,2) + 5*pow(zc,2)) + 16*pow(r1,5)*(pow(RT,3) - 2*RT*pow(zc,2)) + pow(r1,4)*(-42*pow(RT,4) + 44*pow(RT,2)*pow(zc,2) + 15*pow(zc,4)) + 16*pow(r1,3)*(pow(RT,5) - 4*pow(RT,3)*pow(zc,2) - RT*pow(zc,4)) +     2*pow(r1,2)*(6*pow(RT,6) + 22*pow(RT,4)*pow(zc,2) + 17*pow(RT,2)*pow(zc,4) + 2*pow(zc,6))) +      4*(pow(r1,3) + pow(r1,2)*RT + r1*pow(RT,2) + pow(RT,3))*pow(pow(r1,4) - 2*pow(r1,2)*(pow(RT,2) - pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2),1.5)*   EllipticIncompleteE(asin(zc/sqrt(pow(r1,2) - 2*r1*RT + pow(RT,2) + pow(zc,2))),(4*r1*RT)/pow(r1 + RT,2)) -       pow(r1 - RT,2)*(r1 + RT)*pow(pow(r1,4) - 2*pow(r1,2)*(pow(RT,2) - pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2),1.5)*EllipticIncompleteF(asin(zc/sqrt(pow(r1,2) - 2*r1*RT + pow(RT,2) + pow(zc,2))),(4*r1*RT)/pow(r1 + RT,2)))))/ (4.*pow(r1 - RT,4)*pow(r1 + RT,4)*sqrt(pow(pow(r1,2) - pow(RT,2),2))*pow(pow(r1,4) - 2*pow(r1,2)*(pow(RT,2) - pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2),1.5));
return res;
}

double HamakerAtomCircleUnit(double RT, double r1, double re, double zc){
return (RT*(-2*M_PI*sqrt(pow(r1,4) + 2*pow(r1,2)*(-pow(RT,2) + pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2))*(pow(r1,4) + pow(pow(RT,2) + pow(zc,2),2) + 2*pow(r1,2)*(2*pow(RT,2) + pow(zc,2))) +    ((pow(r1 - RT,2) + pow(zc,2))*(pow(r1 + RT,2) + pow(zc,2))*sqrt(-pow(r1,4) + 2*pow(r1,2)*(pow(re,2) + pow(RT,2) - pow(zc,2)) - pow(-pow(re,2) + pow(RT,2) + pow(zc,2),2))* (pow(r1,4) + (pow(RT,2) + pow(zc,2))*(3*pow(re,2) + pow(RT,2) + pow(zc,2)) + pow(r1,2)*(3*pow(re,2) - 2*pow(RT,2) + 2*pow(zc,2))))/pow(re,4)))/   pow(pow(r1,4) + 2*pow(r1,2)*(-pow(RT,2) + pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2),3) + (4*RT*(pow(r1,4) + pow(pow(RT,2) + pow(zc,2),2) + 2*pow(r1,2)*(2*pow(RT,2) + pow(zc,2)))* atan((2*r1*RT*(pow(r1,2) + 2*r1*RT + pow(RT,2) + pow(zc,2))*sqrt((1 - pow(pow(r1,2) - pow(re,2) + pow(RT,2) + pow(zc,2),2)/(4.*pow(r1,2)*pow(RT,2)))/(pow(r1,4) - 2*pow(r1,2)*(pow(RT,2) - pow(zc,2)) + pow(pow(RT,2) + pow(zc,2),2))))/ (-pow(re,2) + pow(r1 + RT,2) + pow(zc,2))))/pow((pow(r1 - RT,2) + pow(zc,2))*(pow(r1 + RT,2) + pow(zc,2)),2.5) ;
}


//returns the "unit" potential for a cylinder of radius RT interacting with an atom at a distance r+RT away. "unit" here means with the Hamaker constant and density of atoms set to unity.
double HamakerAtomTubeUnit(const double RT, const double r, const double re, int useOpt = 1){
double r1 = r +RT; //this defines the offset between the centre of the cylinder and the centre of the amino acid
 
double rho  = 0;
double deltaRho =  0.0001;
double integratedValue = 0;
double tubeHalfLength  = 1000;
//we first check to see if the point is sufficiently far away that no interactions are excluded. if this is the case we have a straightforward expression in terms of elliptic integrals.
if(r1 - RT > re && useOpt==1 ){
integratedValue = TubeExclusionFree(RT, r1, tubeHalfLength) - TubeExclusionFree(RT,r1,-tubeHalfLength) ;//copy over the appropriate expression
} 
else{ 

//if not then we integrate up from zc = 0 to the point at which the exclusion stops matter, then use the analytical result for the remainder
double zc = 0;
double zbound = sqrt(re*re -(r1 - RT)*(r1-RT) );
double deltaz = 0.05;
//std::cout << RT << " " << r1 << " " << re << " zbound: " << zbound <<  "\n";
while(zc < zbound){
integratedValue += deltaz*HamakerAtomCircleUnit(RT, r1, re, zc);
//std::cout << zc << "  " << integratedValue << "\n";
zc+=deltaz;
}
//std::cout << "done numerical integration " << integratedValue << " , beginning analytical";

double remainingTube = TubeExclusionFree(RT,r1,tubeHalfLength) - TubeExclusionFree(RT,r1,zbound); //add on the section between the end of the tube and zbound
//std::cout << " got " << remainingTube << "\n";
integratedValue = 2*(integratedValue+remainingTube); //add on the other half of the cylinder

}



//std::cout << r1 << " " <<  integratedValue << "\n"; //print out the cylinder-atom potential with q1 lambda = 1 for comparison to mathematica
 if(std::isnan(integratedValue)){
std::cout << " unit: found nan at " << r1 << " " << RT << "\n";
}

 

return integratedValue;
}







//hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff
double HamakerSphereTube(const double A, const double R1, const double RT, const double r, const double re) {
 
double integratedValue = HamakerAtomTubeUnit(RT,r,re);
 
//integratedValue is the potential for an atom with lambda q_1 q_2 = 1. we need to do a final step of rescaling to take into account the Hamaker constant and the fact that the atom is actually an amino acid of radius R1 
//to do this, we multiply by lambda q_1 q_2 * volume of R1, then make use of the fact that A = lambda q_1 _q2 M_PI^2 to absorb some of these factors. (note that implicitly the value of M_PI arising from multiplication by the sphere volume has been cancelled due to dividing by M_PI^2)


if(std::isnan(integratedValue)){
std::cout << "found nan at " << r << " " << RT << "\n";
}
return integratedValue/(M_PI  ) * A * (4.0/3.0) * (R1*R1*R1);
}
 

#endif
