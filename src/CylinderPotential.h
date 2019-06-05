#ifndef CYLINDERPOTENTIALS__H__
#define CYLINDERPOTENTIALS__H__

 

#include <cmath>
 

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

 

 


//the Mathematica version of the elliptic integrals uses a different argument to the boost definition. So here we define wrapper functions to convert from the Mathematica to the Boost argument and also match the naming conventions.
double EllipticE(double m){
return boost::math::ellint_2(sqrt(m) );
}
double EllipticK(double m){
return boost::math::ellint_1(sqrt(m) );
}
//these are defined elsewhere but are also needed here so we're being lazy and defining another wrapper function to avoid having to set up include protection

double EllipticIncompleteE2(double phi, double m){
return boost::math::ellint_2(  sqrt(m) , phi);
}
double EllipticIncompleteF2(double phi, double m){
return boost::math::ellint_1(  sqrt(m), phi );
}

//here, r1 is the distance from the atom to the centre of the cylinder, RT is the radius of the cylinder and re is the exclusion region.
//the nanotube is infinitely tall because that's a thing that makes sense for a nanoparticle.
double hamakerAtomInfiniteLine(double r1, double rho, double re, double phi, int useOpt=1){
double piVal = M_PI;
double res = 0;
 if(useOpt!=1){
std::cout << "no unoptimised code yet implemented \n";
}
//two cases are possible: the atom is far enough away that the exclusion is meaningless, or it can matter. 
if(r1*r1 + rho*rho - 2 * r1 * rho * cos(phi) > re*re){
res = (-3*piVal*rho*pow(1/(r1*r1 +rho*rho - 2*r1*rho*cos(phi)),2.5))/8.;


if(std::isnan(res)){
std::cout << phi << " " << r1 << " " << rho << " " << re << " from no-exclusion\n";
}

}
else{

if(phi < 0.000001 && rho - (r1-re) < 0.0001){

res = (3*piVal*sqrt(pow(re,-2))*(-r1 + re))/(8.*pow(re,4));
}
else{
res =  2*((-3*piVal*rho*pow(1/(pow(r1,2) + pow(rho,2) - 2*r1*rho*cos(phi)),2.5))/16. + (rho*(3*pow(re,4)*atan(sqrt(-pow(r1,2) + pow(re,2) - pow(rho,2) + 2*r1*rho*cos(phi))/sqrt(pow(r1,2) + pow(rho,2) - 2*r1*rho*cos(phi))) +   (2*pow(r1,2) + 3*pow(re,2) + 2*pow(rho,2) - 4*r1*rho*cos(phi))*sqrt(pow(r1,2) + pow(rho,2) - 2*r1*rho*cos(phi))*sqrt(-pow(r1,2) + pow(re,2) - pow(rho,2) + 2*r1*rho*cos(phi))))/   (8.*pow(re,4)*pow(pow(r1,2) + pow(rho,2) - 2*r1*rho*cos(phi),2.5)));
}


if(std::isnan(res)){
std::cout << phi << " " << r1 << " " << rho << " " << re << " from  exclusion\n";
}

}
return res;
}


//This function calculates the interaction for a tube (hollow cylinder) of radius rho and infinite height.
double hamakerAtomInfiniteTube(double r1, double rho, double re, int useOpt=1){
double piVal = M_PI;

//here we check to see if the distance of closest approach is greater than the cutoff. if so, we have a simple analytical form to use.

if(r1 - rho > re && useOpt==1){ 
double res =  -(piVal*rho*(4*(pow(r1,2) + pow(rho,2))*EllipticE((4*r1*rho)/pow(r1 + rho,2)) - pow(r1 - rho,2)*EllipticK((4*r1*rho)/pow(r1 + rho,2))))/(2.*pow(r1 - rho,4)*pow(r1 + rho,3));
 
return res;

} 




//if that's not true then we numerically integrate around the circle of radius rho. because of the symmetry we only need to integrate around the half-circle and then double the result.
double phi = 0;

double integratedValue = 0;

//The distance between each point considered is approximately given by rho * deltaPhi , so if rho is large then we're stepping over too large a distance at once, which leads to an increased error when the tube is close to the AA.   

int numIt =  std::nearbyint( M_PI * rho / 0.25) ;
if(numIt < 1000){
numIt = 1000;
}

double deltaPhi = piVal/(numIt-1) ;
for(int i = 0; i<numIt; ++i ){
//double contribution =   deltaPhi * hamakerAtomInfiniteLine(r1,rho,re,phi,useOpt);
phi=i*deltaPhi;
integratedValue +=  deltaPhi * hamakerAtomInfiniteLine(r1,rho,re,phi+deltaPhi/2.0,useOpt); //evaluate at phi+deltaPhi/2.0 to take advantage of the midpoint rule.
//std::cout << phi << " : " << integratedValue << "\n";
}

return 2*integratedValue;



}


//returns the "unit" potential for a cylinder of radius RT interacting with an atom at a distance r+RT away. "unit" here means with the Hamaker constant and density of atoms set to unity.

/*
double HamakerAtomCylinderUnit(const double RT, const double r, const double re, int useOpt = 1){
double r1 = r +RT; //this defines the offset between the centre of the cylinder and the centre of the amino acid
double piVal = M_PI;
double rho  = 0;
double deltaRho =  0.0001;
double integratedValue = 0;
//we first check to see if the point is sufficiently far away that no interactions are excluded. if this is the case we have a straightforward expression in terms of elliptic integrals.
if(r1 - RT > re && useOpt==1 ){
integratedValue =  (piVal*((-pow(r1,2) - 7*pow(RT,2))*EllipticE((4*r1*RT)/pow(r1 + RT,2)) + pow(r1 - RT,2)*EllipticK((4*r1*RT)/pow(r1 + RT,2))))/
   (12.*pow(r1 - RT,3)*pow(r1 + RT,2));
} 
else{ //if not then we numerically integrate. fortunately we do not need to numerically integrate over the entire cylinder, just the part that potentially has excluded interactions. to do this we find the radius of the inner cylinder which is entirely interacting and is given by R1 - re in the interval [0,RT] and use the analyical expression for this.

if(r1 - re > 0 && r1 - re < RT && useOpt==1 ){
double RTEff = r1 - re;
integratedValue = (piVal*((-pow(r1,2) - 7*pow(RTEff ,2))*EllipticE((4*r1*RTEff)/pow(r1 + RTEff,2)) + pow(r1 - RTEff,2)*EllipticK((4*r1*RTEff)/pow(r1 + RTEff,2))))/  (12.*pow(r1 - RTEff,3)*pow(r1 + RTEff,2));
rho = RTEff; //set rho to the final value so we can integrate over the rest afterwards.
 
}
//this then integrates over the rest of the cylinder.
while(rho < RT){
integratedValue += deltaRho * hamakerAtomInfiniteTube(r1,rho+deltaRho/2.0,re,useOpt);
rho+=deltaRho;
}
}
//std::cout << r1 << " " <<  integratedValue << "\n"; //print out the cylinder-atom potential with q1 lambda = 1 for comparison to mathematica
 if(std::isnan(integratedValue)){
std::cout << " unit: found nan at " << r << " " << RT << "\n";
}
return integratedValue;
}

*/

//computes the potential for a disk (filled circle) at a height of zc above the origin of radius RC interacting with an atom located a distance d from the surface of the cylinder (i.e. at a distance RC+d from the origin).
//the Mathematica form as used here is technically undefined for zc = 0 but the limit exists so this is used for zc^2 < 0.0001
double HamakerAtomDiskUnit(const double RC, const double d, const double re, double zc){
double res = 0;
//std::cout << RC << " " << d << " " << re << " " << zc << "\n";
if(re*re < d*d + zc*zc){
//this calculates the potential for a full disk, as for this combination of d and zc there's no influence from the exclusion.
res = (M_PI*(pow(d,6) + 6*pow(d,5)*RC + 12*pow(d,4)*pow(RC,2) + 8*pow(d,3)*pow(RC,3) + 3*pow(d,4)*pow(zc,2) + 12*pow(d,3)*RC*pow(zc,2) + 18*pow(d,2)*pow(RC,2)*pow(zc,2) +  12*d*pow(RC,3)*pow(zc,2) + 3*pow(d,2)*pow(zc,4) + 6*d*RC*pow(zc,4) + 2*pow(RC,2)*pow(zc,4) + pow(zc,6) -    pow((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2)),1.5)))/ (4.*pow(zc,4)*pow((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2)),1.5)); 
}
else if(zc > 0.001){
res =-(sqrt(pow(pow(d,2) + 2*d*RC - pow(re,2) + pow(zc,2),2))*sqrt(-((pow(d,2) - pow(re,2) + pow(zc,2))*(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2)))))/  (4.*pow(re,2)*pow(zc,2)*((d + 2*RC - re)*(d + re) + pow(zc,2))*((d - re)*(d + 2*RC + re) + pow(zc,2))) +   (sqrt(-((pow(d,2) - pow(re,2) + pow(zc,2))*(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2))))* (d*(d + 2*RC)*pow(re,2)*(pow(d,2) + 2*d*RC + 4*pow(RC,2) - pow(re,2)) + (4*d*pow(RC,2)*(d + 2*RC) + 2*d*(d + 2*RC)*pow(re,2) - pow(re,4))*pow(zc,2) + (4*pow(RC,2) + pow(re,2))*pow(zc,4)))/  (4.*pow(re,2)*pow(zc,2)*(pow(d,2) + pow(zc,2))*(pow(d + 2*RC,2) + pow(zc,2))*((d + 2*RC - re)*(d + re) + pow(zc,2))*((d - re)*(d + 2*RC + re) + pow(zc,2))) + (M_PI*(pow(d,6) + 6*pow(d,5)*RC + 12*pow(d,4)*pow(RC,2) + 8*pow(d,3)*pow(RC,3) + 3*pow(d,4)*pow(zc,2) + 12*pow(d,3)*RC*pow(zc,2) + 18*pow(d,2)*pow(RC,2)*pow(zc,2) + 12*d*pow(RC,3)*pow(zc,2) +  3*pow(d,2)*pow(zc,4) + 6*d*RC*pow(zc,4) + 2*pow(RC,2)*pow(zc,4) + pow(zc,6) - pow((pow(d,2) + pow(zc,2))*(pow(d + 2*RC,2) + pow(zc,2)),1.5)))/ (4.*pow(zc,4)*pow((pow(d,2) + pow(zc,2))*(pow(d + 2*RC,2) + pow(zc,2)),1.5)) -  ((pow(re,4) - 2*pow(zc,4))*atan(sqrt(-((pow(d,2) - pow(re,2) + pow(zc,2))*(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2))))/(pow(d,2) + 2*d*RC + 2*pow(RC,2) - pow(re,2) + pow(zc,2))))/  (4.*pow(re,4)*pow(zc,4)) + ((pow(d,3)*pow(d + 2*RC,3) + 3*d*(d + 2*RC)*(pow(d,2) + 2*d*RC + 2*pow(RC,2))*pow(zc,2) + (3*pow(d,2) + 6*d*RC + 2*pow(RC,2))*pow(zc,4) + pow(zc,6))*  atan((pow(d,2) - pow(re,2) + pow(zc,2))/sqrt(-(((pow(d,2) + pow(zc,2))*(pow(d,2) - pow(re,2) + pow(zc,2))*(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2)))/(pow(d + 2*RC,2) + pow(zc,2))))))/  (2.*pow(zc,4)*pow((pow(d,2) + pow(zc,2))*(pow(d + 2*RC,2) + pow(zc,2)),1.5)) +  ((pow(re,4) - pow(zc,4))*atan(sqrt(-(((pow(d,2) - pow(re,2) + pow(zc,2))*(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2)))/pow(pow(d,2) + 2*d*RC - pow(re,2) + pow(zc,2),2)))))/(2.*pow(re,4)*pow(zc,4)) -   (sqrt(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2))*((d - re)*(d + 2*RC - re)*(d + re)*(d + 2*RC + re) + 2*(pow(d,2) + 2*d*RC - pow(RC,2) - pow(re,2))*pow(zc,2) + pow(zc,4))*   atan((pow(d,2) + 2*d*RC - pow(re,2) + pow(zc,2))/sqrt(-(((d + 2*RC - re)*(d + re) + pow(zc,2))*((d - re)*(d + 2*RC + re) + pow(zc,2))))))/ (2.*pow(zc,4)*(pow(d,2) - pow(re,2) + pow(zc,2))*pow((((d + 2*RC - re)*(d + re) + pow(zc,2))*((d - re)*(d + 2*RC + re) + pow(zc,2)))/(pow(d,2) - pow(re,2) + pow(zc,2)),1.5)) + (sqrt(pow(d + 2*RC,2) - pow(re,2) + pow(zc,2))*((d - re)*(d + 2*RC - re)*(d + re)*(d + 2*RC + re) + 2*(pow(d,2) + 2*d*RC - pow(RC,2) - pow(re,2))*pow(zc,2) + pow(zc,4))* atan(sqrt(-1 + RC*re*(1/((d + 2*RC - re)*(d + re) + pow(zc,2)) - 1/((d - re)*(d + 2*RC + re) + pow(zc,2))))))/   (2.*pow(zc,4)*(pow(d,2) - pow(re,2) + pow(zc,2))*pow((((d + 2*RC - re)*(d + re) + pow(zc,2))*((d - re)*(d + 2*RC + re) + pow(zc,2)))/(pow(d,2) - pow(re,2) + pow(zc,2)),1.5));
}
else{
res = HamakerAtomDiskUnit(RC, d,  re,  0.01); //disk is essentially at zc = 0 so need the limiting value
}
return res;

}

//returns the indefinite integral for the cylinder potential integrated from a fukk unit disk evaluated at zc. This has no exclusion and so caution needs to be used to make sure you're not evaluating for values of d and zc corresponding to being inside the cylinder
double HamakerAtomCylinderUnitSegmentIntegrandNE(const double RC, const double d, const double zc){
double res = 0;
if(d>0.01){

res =  -(M_PI*(2*pow(d,6) + 16*pow(d,5)*RC + 48*pow(d,4)*pow(RC,2) + 64*pow(d,3)*pow(RC,3) + 32*pow(d,2)*pow(RC,4) + 2*pow(d,4)*pow(zc,2) + 12*pow(d,3)*RC*pow(zc,2) +      28*pow(d,2)*pow(RC,2)*pow(zc,2) + 32*d*pow(RC,3)*pow(zc,2) + 16*pow(RC,4)*pow(zc,2) + pow(d,2)*pow(zc,4) + 4*d*RC*pow(zc,4) + 24*pow(RC,2)*pow(zc,4) +   (40*pow(RC,3)*pow(zc,4))/d + (32*pow(RC,4)*pow(zc,4))/pow(d,2) + pow(zc,6) + (2*RC*pow(zc,6))/d + (8*pow(RC,2)*pow(zc,6))/pow(d,2) -        2*pow(d,4)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))) -  12*pow(d,3)*RC*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))) -   24*pow(d,2)*pow(RC,2)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))) -   16*d*pow(RC,3)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))) -   (pow(d + 2*RC,2)*(pow(d,2) + 2*d*RC + 8*pow(RC,2))*pow(zc,4)*sqrt((pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))/pow(d + 2*RC,2))*    sqrt((pow(d,2)*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2)))/(pow(d + 2*RC,2)*(pow(d,2) + pow(zc,2))))*sqrt(1 + pow(zc,2)/pow(d,2)))/pow(d,2)))/  (24.*d*pow(d + 2*RC,3)*pow(zc,3)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2)))) -   (sqrt(pow(d,-2))*M_PI*(pow(d,2) + 2*d*RC + 8*pow(RC,2))*sqrt((pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))/pow(d + 2*RC,2))*sqrt(1 + pow(zc,2)/pow(d,2))*   EllipticIncompleteE2(atan(sqrt(pow(d,-2))*zc),1 - pow(d,2)/pow(d + 2*RC,2)))/(24.*d*(d + 2*RC)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2)))) +  (sqrt(pow(d,-2))*d*M_PI*sqrt((pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))/pow(d + 2*RC,2))*sqrt(1 + pow(zc,2)/pow(d,2))*  EllipticIncompleteF2(atan(sqrt(pow(d,-2))*zc),1 - pow(d,2)/pow(d + 2*RC,2)))/(24.*(d + 2*RC)*sqrt((pow(d,2) + pow(zc,2))*(pow(d,2) + 4*d*RC + 4*pow(RC,2) + pow(zc,2))));
}
else{
res = HamakerAtomCylinderUnitSegmentIntegrandNE(RC, 0.01, zc); //at this close range we just set it to use a limiting value to avoid the numerical errors caused in the limit d->0.
}
return res;
}

double HamakerAtomCylinderUnit(const double RC, const double d, const double re, int useOpt = 1){
double cylinderHalfLength = 2000;
if(useOpt!=1){
std::cout << "no unoptimised code yet implemented \n";
}
//there are two cases that need to be considered here.
double res = 0;
//case 1: the atom is sufficiently far away that no part of the cylinder is excluded, requiring d > re . In this case we use the indef. integrand of the disk with no exclusions wrt z and exploit symmetry to evaluate this once.
if(d > re){
res = 2*HamakerAtomCylinderUnitSegmentIntegrandNE(RC, d, cylinderHalfLength);
 
}
else{
//case 2 is more complex.  we integrate up from zc = 0 to zc = sqrt(re^2 - d^2) as this is the region in which the exclusion matters. 
double zcMax = sqrt(re*re - d*d);
double zc = 0;
double deltazc = zcMax/100.0;
while(zc < zcMax){
res += HamakerAtomDiskUnit( RC, d, re,   zc+deltazc/2.0 )*deltazc;
 
zc+=deltazc;
}


//next we use the analytical expression for the exclusion-free region spanning zcMax to the end of the cylinder
res += HamakerAtomCylinderUnitSegmentIntegrandNE(RC, d, cylinderHalfLength) -HamakerAtomCylinderUnitSegmentIntegrandNE(RC, d, zcMax);
//finally we use symmetry to evaluate the other half of the cylinder.
res = res*2;
} 
if(std::isnan(res)){
std::cout << "found nan at " << d << " " << RC << "\n";
}
if(std::isinf(res)){
std::cout << "found inf at " << d << " " << RC << "\n";
}
//std::cout << d << " "  << re << " " << res << "\n";
return res;
}


//hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff
double HamakerSphereCylinder(const double A, const double R1, const double RT, const double r, const double re) {
double piVal = M_PI;
double integratedValue = HamakerAtomCylinderUnit(RT,r,re);
// std::cout << r << " "  << re << " " << integratedValue << "\n";
//integratedValue is the potential for an atom with lambda q_1 q_2 = 1. we need to do a final step of rescaling to take into account the Hamaker constant and the fact that the atom is actually an amino acid of radius R1 
//to do this, we multiply by lambda q_1 q_2 * volume of R1, then make use of the fact that A = lambda q_1 _q2 pi^2 to absorb some of these factors. (note that implicitly the value of pi arising from multiplication by the sphere volume has been cancelled due to dividing by pi^2)


if(std::isnan(integratedValue)){
std::cout << "found nan at " << r << " " << RT << "\n";
}
if(std::isinf(integratedValue)){
std::cout << "found inf at " << r << " " << RT << "\n";
}
return integratedValue/(piVal ) * A * (4.0/3.0) * (R1*R1*R1);
}
 

#endif
