#ifndef CUBEPOTENTIALS__H__
#define CUBEPOTENTIALS__H__

 
#include <cstdio>
#include <random>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iostream>
#include <vector>
#include <cmath>
 
  using namespace std;
 
//this function calculates the potential for a cuboid with part of a sphere subtracted from it. this won't work if the exclusion distance is greater than RT and probably fails in other cases too. 
//the dimensions of the cuboid are [-RT,RT],[-RT,RT],[ZMax,RT], where ZMax is the highest point at which exclusion doesn't matter and is defined by r1 = re + zmax.
double HamakerAtomCuboidExcludedRegion(const double RT, const double r1, const double re){
return -1/(2.*re*pow(RT,2)) + 1/(2.*pow(RT,2)*(r1 + RT)) + 1/(3.*re*(pow(re,2) + pow(RT,2))) + re/(3.*pow(RT,2)*(pow(re,2) + pow(RT,2))) - 1/(3.*(r1 + RT)*(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))) - 
   (r1 + RT)/(3.*pow(RT,2)*(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))) - (5*atan(re/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5)) + 
   (5*atan((r1 + RT)/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5)) + 
   ((-8*M_PI)/pow(re,3) + 2/(re*pow(RT,2)) + (5*sqrt(2)*atan(re/(sqrt(2)*RT)))/pow(RT,3) + 
      (4*(2*pow(r1,4) - 8*pow(r1,3)*(r1 - re) + 2*pow(r1 - re,4) + pow(r1 - re,2)*pow(RT,2) + 2*pow(RT,4) - 2*r1*(r1 - re)*(4*pow(r1 - re,2) + pow(RT,2)) + pow(r1,2)*(12*pow(r1 - re,2) + pow(RT,2)))*
         atan(RT/sqrt(pow(r1,2) - 2*r1*(r1 - re) + pow(r1 - re,2) + pow(RT,2))))/(pow(re,3)*pow(RT,3)*sqrt(pow(r1,2) - 2*r1*(r1 - re) + pow(r1 - re,2) + pow(RT,2))))/12. + 
   atan(RT/sqrt(pow(re,2) + pow(RT,2)))/(2.*re*RT*sqrt(pow(re,2) + pow(RT,2))) - 
   ((2*pow(re,4) + pow(re,2)*pow(RT,2) + 2*pow(RT,4))*atan(RT/sqrt(pow(re,2) + pow(RT,2))))/(6.*pow(re,3)*pow(RT,3)*sqrt(pow(re,2) + pow(RT,2))) - 
   ((3*pow(re,4) + 3*pow(re,2)*pow(RT,2) + 2*pow(RT,4))*atan(RT/sqrt(pow(re,2) + pow(RT,2))))/(6.*pow(re,3)*RT*pow(pow(re,2) + pow(RT,2),1.5)) - 
   ((2*pow(re,4) + 3*pow(re,2)*pow(RT,2) + 3*pow(RT,4))*atan(RT/sqrt(pow(re,2) + pow(RT,2))))/(6.*re*pow(RT,3)*pow(pow(re,2) + pow(RT,2),1.5)) + 
   ((2*M_PI)/pow(r1 - RT,3) + (6*M_PI*(r1 - RT))/pow(re,4) + 4/((r1 - RT)*pow(RT,2)) + 6/(pow(RT,2)*(-r1 + RT)) + (3*sqrt(2)*atan((r1 - RT)/(sqrt(2)*RT)))/pow(RT,3) + 
      (8*sqrt(2)*atan((-r1 + RT)/(sqrt(2)*RT)))/pow(RT,3) + (4*(2*pow(r1,4) - 8*pow(r1,3)*RT + 13*pow(r1,2)*pow(RT,2) - 10*r1*pow(RT,3) + 5*pow(RT,4))*
         atan(RT/sqrt(pow(r1,2) - 2*r1*RT + 2*pow(RT,2))))/(pow(RT,3)*pow(-r1 + RT,3)*sqrt(pow(r1,2) - 2*r1*RT + 2*pow(RT,2))))/12. - 
   atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2)))/(2.*RT*(r1 + RT)*sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))) + 
   ((2*pow(r1,4) + 8*pow(r1,3)*RT + 15*pow(r1,2)*pow(RT,2) + 14*r1*pow(RT,3) + 8*pow(RT,4))*atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))))/
    (6.*pow(RT,3)*(r1 + RT)*pow(pow(r1,2) + 2*r1*RT + 2*pow(RT,2),1.5)) + ((3*pow(r1,4) + 12*pow(r1,3)*RT + 21*pow(r1,2)*pow(RT,2) + 18*r1*pow(RT,3) + 8*pow(RT,4))*
      atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))))/(6.*RT*pow(r1 + RT,3)*pow(pow(r1,2) + 2*r1*RT + 2*pow(RT,2),1.5)) + 
   ((2*pow(r1,4) + pow(r1,2)*pow(RT,2) + 8*r1*pow(RT,3) + 4*pow(RT,4) + 2*r1*RT*(4*pow(r1,2) + pow(RT,2)) + pow(RT,2)*(12*pow(r1,2) + pow(RT,2)))*atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))))/
    (6.*pow(RT,3)*pow(r1 + RT,3)*sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2)));
}


//this function returns the potential integrated over the cube defined by the boundaries ( [-RT,RT], [-RT,RT], [-RT,ZMax] ) for an atom located at (0,0,r1) with no exclusion

double HamakerAtomCuboidNoExclusion(const double RT, const double r1, const double ZMax){
return (RT + ZMax)/(6.*pow(RT,2)*(r1 + RT)*(-r1 + ZMax)) + (5*atan((r1 + RT)/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5)) + 
   ((2*pow(r1,4) + 8*pow(r1,3)*RT + 13*pow(r1,2)*pow(RT,2) + 10*r1*pow(RT,3) + 5*pow(RT,4))*atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))))/
    (3.*pow(RT,3)*pow(r1 + RT,3)*sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))) + ((-2*pow(RT,4) - pow(RT,2)*pow(r1 - ZMax,2) - 2*pow(r1 - ZMax,4))*atan(RT/sqrt(pow(RT,2) + pow(r1 - ZMax,2))))/
    (3.*pow(RT,3)*sqrt(pow(RT,2) + pow(r1 - ZMax,2))*pow(r1 - ZMax,3)) - (5*atan((r1 - ZMax)/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5));
}
 

double HamakerAtomCubeUnit(const double RT, const double r, const double re, int useOpt = 1){
double r1 = r +RT; //this defines the offset between the centre of the cube and the centre of the amino acid

 
double valueFromNonExcluded = 0;

//we first check to see if the point is sufficiently far away that no interactions are excluded. if this is the case we have a straightforward expression:
if(r1 - RT > re && useOpt==1 ){

valueFromNonExcluded =  1/(-3*pow(r1,2)*RT + 3*pow(RT,3)) - (5*atan((r1 - RT)/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5)) + 
   (5*atan((r1 + RT)/(sqrt(2)*sqrt(pow(RT,2)))))/(6.*sqrt(2)*pow(pow(RT,2),1.5)) + 
   ((2*pow(r1,4) - 8*pow(r1,3)*RT + 13*pow(r1,2)*pow(RT,2) - 10*r1*pow(RT,3) + 5*pow(RT,4))*atan(RT/sqrt(pow(r1,2) - 2*r1*RT + 2*pow(RT,2))))/
    (3.*pow(RT,3)*pow(-r1 + RT,3)*sqrt(pow(r1,2) - 2*r1*RT + 2*pow(RT,2))) + ((2*pow(r1,4) + 8*pow(r1,3)*RT + 13*pow(r1,2)*pow(RT,2) + 10*r1*pow(RT,3) + 5*pow(RT,4))*
      atan(RT/sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))))/(3.*pow(RT,3)*pow(r1 + RT,3)*sqrt(pow(r1,2) + 2*r1*RT + 2*pow(RT,2))) ;

} 
else{ 


double ZMax = r1 - re;
if(ZMax > RT){
ZMax = RT;
}
 
valueFromNonExcluded =   HamakerAtomCuboidExcludedRegion(RT,r1,re);

}

//if the particle is nearly on the surface of the cube then the maths breaks down a bit.
//the limit exists and so we take that value.
if((r1-RT)*(r1-RT) < 0.0001 ){
valueFromNonExcluded = (-40*M_PI*pow(RT,3) + pow(re,3)*(5 + 25*sqrt(2)*atan(sqrt(2)) + 19*sqrt(5)*atan(1/sqrt(5))))/(60.*pow(re,3)*pow(RT,3));
}


 if(std::isnan(valueFromNonExcluded)){
std::cout << " unit: found nan at " << r << " " << RT << "\n";
}
return valueFromNonExcluded ;

}







//hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff
double HamakerSphereCube(const double A, const double R1, const double RT, const double r, const double re) {
double piVal = M_PI;
double integratedValue = HamakerAtomCubeUnit(RT,r,re);
 
//integratedValue is the potential for an atom with lambda q_1 q_2 = 1. we need to do a final step of rescaling to take into account the Hamaker constant and the fact that the atom is actually an amino acid of radius R1 
//to do this, we multiply by lambda q_1 q_2 * volume of R1, then make use of the fact that A = lambda q_1 _q2 pi^2 to absorb some of these factors. (note that implicitly the value of pi arising from multiplication by the sphere volume has been cancelled due to dividing by pi^2)

 
if(std::isnan(integratedValue)){
std::cout << "found nan at " << r << " " << RT << "\n";
}
return integratedValue/(piVal ) * A * (4.0/3.0) * (R1*R1*R1);
}
 

#endif
