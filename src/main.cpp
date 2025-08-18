#define PARALLEL

#include "Config.h"
#include "CommandLineParser.h"
#include "ConfigFileReader.h"
#include "HamakerConstants.h"
#include "TargetList.h"
#include "Surface.h"
#include "PDBFile.h"
#include "LigandFileReader.h"

#include "NPTargetList.h"
#include "NPFile.h"

#include "Potentials.h"

#include "RelaxClasses.h"

#include <omp.h>
#include <cmath>
#include <random>
#include <iomanip>
#include <sys/stat.h>
 
#include <chrono>
#include <ctime>    


//constexpr double        gds             = 0.22;
constexpr double        gds             = 0.0;
constexpr double        delta           = 2.0;

bool isUsingDeltaOverride = false;
//define the angular resolution. These defaults are overrwritten if and only if m_confirmOverrideAngle is set by the appropriate switch AND a new value for angle_delta_deg is supplied in the config file.
double        angle_delta_deg = 5.0;
double        angle_delta     = angle_delta_deg * (M_PI / 180.0);
int           ncols           = 72;
int           nrows           = 36;
int           iterations      = nrows * ncols;
int           samples         = 128; //number of samples per angle bin , overwritten by configuration
int           threadPrintFreq = 1; //print out a symbol every N angles as a crude progress bar

constexpr int           steps           = 512; //this is the resolution along the z-axis used for integration. There is no reason to modify this; the default value of 512 produces a resolution of 0.004 nm and the input potentials are less precise than this (or just thermal fluctuations at this resolution)
constexpr double        dz              = delta / (steps - 1); //for non-uniform NPs this is updated for the integration
//define the Boltzmann and Avogadro constants for energy conversions
constexpr double        kbConst       =  1.380649e-23;
constexpr double        naConst       =  6.02214076e23; 
//m_confirmOverrideAngle


std::random_device randomEngineSeed{};
std::mt19937 randomEngine{ randomEngineSeed()  };
//std::uniform_real_distribution<double> random_angle_offset(0.0, angle_delta); 
std::uniform_real_distribution<double> random_angle_offset(0.0, 1.0);

std::normal_distribution<double> unitNormalDist(0.0, 1.0);

//Define the version number to provide metadata, following the semantic versioning convention https://semver.org 
//In brief: Increment the first number if UA is no longer backwards compatible, i.e. old input files won't work or the output files won't work with pre-existing scripts.
//          Increment the second number if there's new functionality, e.g. new NP shapes, new output, etc.
//          Increment the third number if no functionality has been modified but bug fixes have been applied which may change output.
//When you increment a number, all the following numbers should be reset to zero. E.g. If we're at 1.2.3 and a bug fix is applied, move to 1.2.4 , if we then add new functionality, 1.3.0, then a new version entirely, 2.0.0 

std::string getUAVersion(){
    static std::string uaVersionID("1.5.2"); 
    return uaVersionID;
}




bool CheckForPrecalculated(const double *adsorption_energy, const double *adsorption_error, const double *mfpt_val, const double *minloc_val, const double radius,
        const double zeta, const std::string& name, const std::string& directory, const std::string& npName, int isMFPT=0, int appendAngle = 0, double omegaAngle=0) {
std::string filename;
std::string cylinderFileNameAppend;

if(appendAngle == 0){
cylinderFileNameAppend = "";
} 
else{
cylinderFileNameAppend = "_"  + std::to_string(static_cast<int>(omegaAngle));
}


    if(isMFPT==1){
     filename = directory + "/" + npName+"/"+TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  "_mfpt.uam";
    }
else{
    filename = directory + "/" + npName+"/"+TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  ".uam";
}
    //std::clog << "Info: testing for map at: "<< filename << "\n";
struct stat buffer;
return (stat (filename.c_str(), &buffer) == 0); 

}


// minloc_val  gives the distance between the protein's COM and the surface of the NP, so the com_com min is minlocval + radius
// surf-surf separation is minlocval -  protein_offset[i] 
// and protein surf - np centre is minlocval + radius - protein_offset[i]
void WriteMapFile(const double *adsorption_energy, const double *adsorption_error, const double *mfpt_val, const double *minloc_val, const double *numcontacts_val, const double *protein_offset, const double radius,
        const double zeta, const std::string& name, const std::string& directory,  const std::string& npName ,  double temperature=300.0,  int isMFPT=0, int appendAngle = 0, double omegaAngle=0) {
std::string filename;
std::string cylinderFileNameAppend;


auto outputTimeStamp = std::chrono::system_clock::now();
std::time_t end_time = std::chrono::system_clock::to_time_t(outputTimeStamp);
//If directory doesn't exist, create it.
//boost::filesystem::create_directory(directory);

if(appendAngle == 0){
cylinderFileNameAppend = "";
} 
else{
cylinderFileNameAppend = "_"  + std::to_string(static_cast<int>(omegaAngle));
}


    if(isMFPT==1){
     filename = directory + "/" + npName+"/"+TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  "_mfpt.uam";
    }
else{
    filename = directory + "/" + npName+"/"+ TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  ".uam";
}
    std::clog << "Info: Saving map to: "<< filename << "\n";
    std::ofstream handle(filename.c_str());
    double phi, theta;
    std::string timestamp;
    timestamp = std::ctime(&end_time);
    timestamp.erase( timestamp.end() - 1);
    handle << "#Results generated at: " << timestamp  << " using UA version: " << getUAVersion() << "\n";
    handle << "#" << npName << " - " << name << "\n";
    handle << "#RADIUS: "<<radius<<"\n";





    handle << "#phi-LeftHandEdge theta-LeftHandEdge EAds/kbT=300 SDEV(Eads)/kbT=300 min_surf-surf-dist/nm mfpt*DiffusionCoeff/nm^2 EAds/kJ/mol min_ProtSurf_NPCentre-dist/nm omega NumContacts min_COM_COM-dist/nm\n"; 
    for (int i = 0; i < iterations; ++i) { 
        phi   = (i % ncols) * angle_delta_deg;
        theta = (i / ncols) * angle_delta_deg;
        handle << std::left << std::setw(7) << std::fixed << std::setprecision(1) << phi;
        handle << std::left << std::setw(7) << std::fixed << std::setprecision(1) << theta;

        
        //handle << std::left << std::setw(14) << std::fixed << std::scientific << std::setprecision(5) << adsorption_energy[i];
        
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << adsorption_energy[i];
        
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << adsorption_error[i];

        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << minloc_val[i]  - protein_offset[i];

        handle << std::left << std::setw(14) << std::fixed << std::scientific << std::setprecision(5) << mfpt_val[i];
        
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << adsorption_energy[i] * 300.0 * kbConst * naConst / 1000.0;
        
                handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << minloc_val[i] + radius - protein_offset[i] ; //protein surface is at approx. -protein_offset, add radius to get to NP centre

        handle << std::left << std::setw(7) << std::fixed << std::setprecision(1) << omegaAngle;
        
        
         handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << numcontacts_val[i] ;
         handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << minloc_val[i] + radius ; //minloc is np_surf to protein_com, so adjust for that
        handle << "\n";
    }
    handle.close();
}

void PrintStatistics(const double *adsorption_energy, const double *adsorption_error, const double radius, const std::string& filename) {
    double          sin_theta;
    double          error           = 0.0;
    double          simpleAverage   = 0.0;
    long double     T               = 0.0;
    long double     Z               = 0.0;
    double sinThetaTot = 0.0;
    double minEnergy = 50.0;
    for (int i = 0; i < iterations; ++i) {
        sin_theta       = std::sin((i / ncols) * angle_delta);
        T               += adsorption_energy[i] * sin_theta * std::exp(-1.0 * adsorption_energy[i]);
        Z               += sin_theta * std::exp(-1.0 * adsorption_energy[i]);
        simpleAverage   += sin_theta * adsorption_energy[i];
        error           += sin_theta * adsorption_error[i];
        sinThetaTot     += sin_theta;
        if(adsorption_energy[i] < minEnergy){
        minEnergy = adsorption_energy[i];
        }
    }

    simpleAverage   /= sinThetaTot;
    error           /= iterations;
    double boltzAverage = T/Z;
    if(isnan(boltzAverage)){
    boltzAverage = minEnergy;
    }


    std::cout << std::setw(10) << std::left << TargetList::Filename(filename);
    std::cout << std::setw(10) << std::fixed << std::setprecision(1) << radius;
    std::cout << std::setw(14) << std::fixed << std::setprecision(5) << simpleAverage;
    std::cout << std::setw(14) << std::fixed << std::setprecision(5) << boltzAverage;
    std::cout << std::setw(14) << std::fixed << std::setprecision(5) << error;
    std::cout << std::endl;
}
//Phi,theta here are adjusted from the  phi,theta defining the orientation by phi -> -phi, theta -> pi - theta
//in terms of the angles given here, this corresponds to a rotation around the z axis of phi, followed by a rotation around y of theta.
inline void Rotate (const int size, const double phi, const double theta,
const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz,
double *x, double *y, double *z) {
    double R[3][3];
    
    R[0][0] = std::cos(theta) * std::cos(phi);
    R[0][1] = -1.0 * std::cos(theta) * std::sin(phi);
    R[0][2] = std::sin(theta);
    R[1][0] = std::sin(phi);
    R[1][1] = std::cos(phi);
    R[2][0] = -1.0 * std::sin(theta) * std::cos(phi);
    R[2][1] = std::sin(theta) * std::sin(phi);
    R[2][2] = std::cos(theta);

    for(int i = 0; i < size; ++i) {
        x[i] = cx[i] * R[0][0] + cy[i] * R[0][1] + cz[i] * R[0][2];
        y[i] = cx[i] * R[1][0] + cy[i] * R[1][1];
        z[i] = cx[i] * R[2][0] + cy[i] * R[2][1] + cz[i] * R[2][2];





    }


}
//as the above function, except followed by a rotation  of omega around the z axis
inline void Rotate3 (const int size, const double phi, const double theta, const double omega,
const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz,
double *x, double *y, double *z) {
    double R[3][3];

    R[0][0] = std::cos(theta) * std::cos(phi) * std::cos(omega) - std::sin(omega) * std::sin(phi);
    R[0][1] = -1.0 * std::cos(theta) * std::sin(phi) * std::cos(omega) - std::cos(phi)*std::sin(omega);
    R[0][2] = std::sin(theta)*std::cos(omega);
    R[1][0] = std::sin(phi)*std::cos(omega) + std::cos(phi) * std::cos(theta)*std::sin(omega);
    R[1][1] = std::cos(phi)*std::cos(omega) - std::cos(theta)*std::sin(omega)*std::sin(phi);
    R[1][2] = std::sin(omega) * std::sin(theta);
    R[2][0] = -1.0 * std::sin(theta) * std::cos(phi);
    R[2][1] = std::sin(theta) * std::sin(phi);
    R[2][2] = std::cos(theta);

    for(int i = 0; i < size; ++i) {
        x[i] = cx[i] * R[0][0] + cy[i] * R[0][1] + cz[i] * R[0][2];
        y[i] = cx[i] * R[1][0] + cy[i] * R[1][1] + cz[i] * R[1][2];
        z[i] = cx[i] * R[2][0] + cy[i] * R[2][1] + cz[i] * R[2][2];





    }


}



inline void JitterPDB (const int size, const double jitterMagnitude, double *x, double *y, double *z) {
/*
    R[0][0] = std::cos(theta) * std::cos(phi) * std::cos(omega) - std::sin(omega) * std::sin(phi);
    R[0][1] = -1.0 * std::cos(theta) * std::sin(phi) * std::cos(omega) - std::cos(phi)*std::sin(omega);
    R[0][2] = std::sin(theta)*std::cos(omega);
    R[1][0] = std::sin(phi)*std::cos(omega) + std::cos(phi) * std::cos(theta)*std::sin(omega);
    R[1][1] = std::cos(phi)*std::cos(omega) - std::cos(theta)*std::sin(omega)*std::sin(phi);
    R[1][2] = std::sin(omega) * std::sin(theta);
    R[2][0] = -1.0 * std::sin(theta) * std::cos(phi);
    R[2][1] = std::sin(theta) * std::sin(phi);
    R[2][2] = std::cos(theta);
*/
    //std::normal_distribution<double> jitterDist(0.0, jitterMagnitude);


    for(int i = 0; i < size; ++i) {
        x[i] = x[i]  + unitNormalDist(randomEngine)*jitterMagnitude;
        y[i] = y[i] + unitNormalDist(randomEngine)*jitterMagnitude;
        z[i] = z[i] + unitNormalDist(randomEngine)*jitterMagnitude;

    }


}





//rotate the PDB anticlockwise by omega around the z-axis
inline void RotateZ (const int size, const double omega,
const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz,
double *x, double *y, double *z) {
    double R[3][3];

    R[0][0] =  std::cos(omega);
    R[0][1] = - std::sin(omega);
    R[0][2] = 0;
    R[1][0] = std::sin(omega);
    R[1][1] = std::cos(omega);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;

    for(int i = 0; i < size; ++i) {
        x[i] = cx[i] * R[0][0] + cy[i] * R[0][1] + cz[i] * R[0][2];
        y[i] = cx[i] * R[1][0] + cy[i] * R[1][1];
        z[i] = cx[i] * R[2][0] + cy[i] * R[2][1] + cz[i] * R[2][2];





    }


}


void SquareXY (const int size, double *x, double *y) {
    for (int i = 0; i < size; ++i) {
        x[i] *= x[i];
        y[i] *= y[i];
    }
}

void ShiftZ (const int size, double *z) {
    const double z_shift = *std::min_element(z, z + size);
    for (int i = 0; i < size; ++i) {
         // cout << z[i] << "->" <<  z[i] - z_shift << "\n";
        z[i] -= z_shift;
    }
}

void Sum (const int size, double* result, double *arr) {
    *result = 0.0;
    for (int i = 0; i < size; ++i) {
        *result += arr[i];
    }
}
//note that for this integration routine and all others, "ssd" was the distance between the centre of the NP and the closest point of the protein.
//as of Feb 2024, it now refers to the COM-COM distance for consistency with small NPs/large biomolecules
void Integrate(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption,double temperature) {
    long double area = 0.0;
    double minEnergy = 500.0;
    for (int i = 0; i < size; ++i) {
        double energyDiff = energy[i] - init_energy;
        if(energyDiff < minEnergy){
        minEnergy = energyDiff;
        }
        area += static_cast<long double>(ssd[i] * ssd[i] * dz) * std::exp(static_cast<long double>(-1.0 * (energyDiff))); 
    }
    const double factor = 3.0 / std::fabs(std::pow(ssd[0], 3.0) - std::pow(ssd[size - 1], 3.0));
   

    double adsorptionEnergy = -1.0 * (temperature/300.0) * std::log(factor * area);
    if(isnan(adsorptionEnergy)){
    adsorptionEnergy = minEnergy;
    }
    *adsorption = adsorptionEnergy;


}


/*


*/

void IntegrateMFPT(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *mfpt, int npType, int calculateMFPT,double temperature) {
    long double res = 0.0;
    long double innerRes = 0.0;
    long double outerRes = 0.0;

if(calculateMFPT!=0){

//the ssd list contains the positions, note that the ssd[0] corresponds to the furthest distance and ssd[size] is the closest. 

int minimumRange = 0; //number of steps to go away from the minimum when averaging over the initial position
int stepSize = 5;
int correctionFactorPower = 0;
if(npType == 1){
correctionFactorPower = 2;
}
else if(npType == 2 || npType == 4 || npType==5){
correctionFactorPower = 1;
}

//find the location of the deepest minima. the integral defining the MFPT is most strongly peaked around here.
int minI = 0;
double minEnergy = energy[0];
for(int i =0;i<size;++i){
if(energy[i] < minEnergy){
minEnergy  = energy[i];
minI = i;
}
}



for(int k = std::max(minI-minimumRange*stepSize,0); k < std::min(size,minI+minimumRange*stepSize+1); k+=stepSize){
//integrate from the end point at ssd[0] to r0 at ssd[k]
    res = 0.0;
    for (int i = 0; i < k; ++i) { //integrate over "y", runs from r0 to the escape point
        //here ssd[i] = y and ssd[j] = rprime
        innerRes = 0;
        for(int j =size; j > i; --j){ //integrate over "x": from touching the NP to the  value of y
            innerRes +=  dz*std::exp(static_cast<long double>(energy[i]  -1.0*(  energy[j] )))  * pow(ssd[j]/ ssd[i], correctionFactorPower ); 
        }
        res+= dz    * innerRes;
    } //at this point res contains the MFPT for a given initial position
    outerRes += dz*res* pow(ssd[k], correctionFactorPower )*std::exp(static_cast<long double>(  -1.0*(  energy[k] ))); //weight the MFPT by the probability for that initial position
}
double z = 0;
for(int i=std::max(minI-minimumRange*stepSize,0);i<std::min(size,minI+minimumRange*stepSize+1);i+=stepSize){
 z+= dz* pow(ssd[i], correctionFactorPower )*std::exp(static_cast<long double>(  -1.0*(  energy[i] )));
 //calculate the partition function, used to normalise the average MFPT
}
     *mfpt = outerRes/z;

}
else{
*mfpt = -1;
}





}



void IntegrateCylinder(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption,double temperature) {
    long double area = 0.0;
    for (int i = 0; i < size; ++i) {

if(std::isnan(energy[i])){
std::cout << "Cylinder integration: found nan at " << i << " " <<   "\n";
}
if(std::isinf(energy[i])){
std::cout << "Cylinder integration: found inf at " << i << " " <<   "\n";
}

        area += static_cast<long double>(ssd[i]  * dz) * std::exp(static_cast<long double>(-1.0 * (energy[i] - init_energy))); 
    }
    const double factor = 2.0 / std::fabs(std::pow(ssd[0], 2.0) - std::pow(ssd[size - 1], 2.0));

    if(factor*area < 0){ 
 std::cout << "Cylinder integration:  warning: factor*area < 0, unphysical result " << factor << " " << area << "\n";    
}



    *adsorption = -1.0 * (temperature/300.0)  * std::log(factor * area);
}


void IntegrateCube(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption,double temperature) {
    long double area = 0.0;
    for (int i = 0; i < size; ++i) {
        area += static_cast<long double>(dz  ) * std::exp(static_cast<long double>(-1.0 * (energy[i] - init_energy))); 
    }
    const double factor = 1.0 / std::fabs(std::pow(ssd[0], 1.0) - std::pow(ssd[size - 1], 1.0));
    *adsorption = -1.0  * (temperature/300.0) * std::log(factor * area);
}



void MeanAndSD(const int size, double *mean, double *sd, double *arr) {
    double lmean = 0;
    double lsd   = 0;
    for (int i = 0; i < size; ++i) {
        lmean += arr[i];
    }
    lmean /= size;
    for (int i = 0; i < size; ++i) {
        lsd += std::pow(arr[i] - lmean, 2.0);
    }
    lsd /= size;
    
    *mean   = lmean;
    *sd     = std::sqrt(lsd);
 
}



void BoltzMeanAndSD(const int size, double *mean, double *sd, double *arr) {
    double lmeanN = 0;
    double lmeanDN = 0.000001;
    double lmean = 0;
    double lsd   = 0;
    double referenceEnergy =  (*std::min_element(arr,arr + size));
    for (int i = 0; i < size; ++i) {
        double BoltzFactor = std::exp(static_cast<long double>(-1.0 * (arr[i] - referenceEnergy))) ; //use local weighting for the factor such that the exp argument remains finite.
        lmeanN += arr[i] * BoltzFactor  ;
        lmeanDN += BoltzFactor;
    }
    lmean =  lmeanN/lmeanDN;
    for (int i = 0; i < size; ++i) {
        lsd += std::pow(arr[i] - lmean, 2.0);
    }
    lsd /= size;
    
    *mean   = lmean;
    *sd     = std::sqrt(lsd);
 
}



//geometry -specific radius for a bead at (sqrt(x),sqrt(y),z) at an NP-protein distance start.
double getDistance(double x, double y, double z, double start, double radius,int npType){
double distance;


                if(npType == 2 || npType == 4 || npType==5){
                //for cylinder NPs we only take into consideration the radial distance with the cylinder aligned along the x-axis.
                distance    = std::sqrt(y +  (z + start) * (z + start)) - radius; // AA Center To NP Surface Distance (originally stop, changed to start)
                }
                else if(npType == 3){
                //for cubes we consider only the "vertical" distance, i.e. we assume that each bead is approximately at the centre of the cube
                distance = std::sqrt(  (z +start) * (z + start)) - radius;
                }
                else{
                distance    = std::sqrt(x + y + (z + start) * (z + start)) - radius; // Center To Surface Distance
//               std::cout << i << " "  << phi << " "  << theta << " " <<  x[i] << " " << y[i] << " " << z[i] << " " <<distance << "\n";

                }

return distance;
}

//get the distance between the centre of an AA bead and the surface of an NP
double getDistance3D(double aaX, double aaY, double aaZ, double npX, double npY, double npZ, double npRadius, int npType){
double distance = 0;

                if(npType == 2 || npType == 4 || npType==5){
                //for cylinder NPs we only take into consideration the radial distance with the cylinder aligned along the x-axis.
                distance    = std::sqrt((aaY - npY)*(aaY-npY) +     (aaZ - npZ)*(aaZ-npZ) ) - npRadius; // AA Center To NP Surface Distance (originally stop, changed to start)
                }
                else if(npType == 3){
                //previously: for cubes we consider only the "vertical" distance, i.e. we assume that each bead is approximately at the centre of the cube
                //distance = std::sqrt( (aaZ - npZ)*(aaZ-npZ)   ) - npRadius;
                //double cubeDelta = std::min(npR, 1.0);
                //updated algorithm: x,y components only count if they are further than edge - 1nm 
                //double ydist =     std::max( 0, aaY - (npY + npRadius - cubeDelta ) ) ;
                //double ydist1 = std::max(0,  aaY - (npY + npRadius - cubeDelta) );
                //double ydist2 = std::min(0,  aaY - (npY - npRadius + cubeDelta) );
                //double ycomponentsq = std::max(  ydist1*ydist1, ydist2*ydist2);

                //double xdist1 = std::max(0,  aaX - (npX + npRadius - cubeDelta) );
                //double xdist2 = std::min(0,  aaX - (npX - npRadius + cubeDelta) );
                //double xcomponentsq = std::max(  xdist1*xdist1, xdist2*xdist2);

                distance  = std::sqrt(  std::pow(max(0.0, std::abs(aaX-npX ) - npRadius )  ,2) +  std::pow(max(0.0, std::abs(aaY-npY ) - npRadius )  ,2) +   std::pow(max(0.0, std::abs(aaZ-npZ ) - npRadius )  ,2)   ) ;

                }
                else{
                distance    = std::sqrt((aaX - npX)*(aaX-npX) + (aaY - npY)*(aaY-npY) +     (aaZ - npZ)*(aaZ-npZ)   ) - npRadius; // sphere surface to x,y,z distance
 

                }

return distance;

}


double GetEnergy(double *x, double *y, double *z, int j, double dr, bool sumNPs, int numNPs, int size, const PDB& pdb,const NP& np,const Potentials& potentials, double radius, int npType){
       double totalEnergy = 0;
       for(int i = 0; i< size; ++i){
             for(int k=0; k<numNPs; ++k){
              if(sumNPs == false){
                  double dist = getDistance3D(x[i],y[i],z[i] ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
                  totalEnergy += ( static_cast<double>(potentials[ pdb.m_id[i]].Value(dist, np.m_npBeadType[k] ))) ;
              }
              else{
                  double dist = getDistance3D(x[i],y[i],z[i] ,   0 , 0 , 0 , radius,npType);
                  totalEnergy += ( static_cast<double>(potentials[ pdb.m_id[i]].Value(dist, np.m_npBeadType[k] ))) ;
              }
             }
       }
    return totalEnergy;
}


double GetPotGrad(double *x, double *y, double *z, int j, double dr, bool sumNPs, int numNPs, int dir, const PDB& pdb,const NP& np,const Potentials& potentials, double radius, int npType)
{
           double dudx = 0;
           double dx = (dir == 0) ? dr/2:0.0;
           double dy = (dir == 1) ? dr/2:0.0;
           double dz = (dir == 2) ? dr/2:0.0;
           double dx1 = 0;
           double dx0 = 0;

           for(int k=0; k<numNPs; ++k){

              // ForceX  = -d U( h(x,y,z) )/dx
              // ForceX  = -  dU/dh * dh/dx
              if(sumNPs == false){
              dx1 = getDistance3D(x[j]+dx,y[j]+dy,z[j]+dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              dx0 = getDistance3D(x[j]-dx,y[j]-dy,z[j]-dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);

              dudx += ( static_cast<double>(potentials[ pdb.m_id[j]].Value(dx1, np.m_npBeadType[k] )) - static_cast<double>(potentials[ pdb.m_id[j]].Value(dx0, np.m_npBeadType[k] )) )/dr ;
              }
              else{

              dx1 = getDistance3D(x[j]+dx,y[j]+dy,z[j]+dz ,   0 , 0 , 0 , radius,npType);
              dx0 = getDistance3D(x[j]-dx,y[j]-dy,z[j]-dz ,   0 , 0 , 0 , radius,npType);
              dudx += ( static_cast<double>(potentials[ pdb.m_id[j]].Value(dx1, np.m_npBeadType[k] )) - static_cast<double>(potentials[ pdb.m_id[j]].Value(dx0, np.m_npBeadType[k] )) )/dr ;
               }
            }
            return dudx;
}

double GetPotGrad4(double *x, double *y, double *z, int j, double dr, bool sumNPs, int numNPs, int dir, const PDB& pdb,const NP& np,const Potentials& potentials, double radius, int npType)
{
           double dudx = 0;
           double dx = (dir == 0) ? dr/2:0.0;
           double dy = (dir == 1) ? dr/2:0.0;
           double dz = (dir == 2) ? dr/2:0.0;
           double dxp1 = 0;
           double dxp2 = 0;
           double dxm1 = 0;
           double dxm2 = 0;
           double ep2 = 0;
           double ep1 = 0;
           double em1 = 0;
           double em2 = 0;

           for(int k=0; k<numNPs; ++k){

              // ForceX  = -d U( h(x,y,z) )/dx
              // ForceX  = -  dU/dh * dh/dx
              if(sumNPs == false){
              dxp2 = getDistance3D(x[j]+2*dx,y[j]+2*dy,z[j]+2*dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              dxp1 = getDistance3D(x[j]+dx,y[j]+dy,z[j]+dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              dxm1 = getDistance3D(x[j]-dx,y[j]-dy,z[j]-dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              dxm2 = getDistance3D(x[j]-2*dx,y[j]-2*dy,z[j]-2*dz ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              ep2 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxp2, np.m_npBeadType[k] ));
              ep1 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxp1, np.m_npBeadType[k] ));
              em1 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxm1, np.m_npBeadType[k] ));
              em2 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxm2, np.m_npBeadType[k] ));


              dudx += ( -1 *ep2 + 8*ep1 - 8*em1 + em2  )/(12.0 * dr) ;
              }
              else{
              dxp2 = getDistance3D(x[j]+2*dx,y[j]+2*dy,z[j]+2*dz ,   0 , 0 , 0 , radius,npType);
              dxp1 = getDistance3D(x[j]+dx,y[j]+dy,z[j]+dz ,   0 , 0 , 0 , radius,npType);
              dxm1 = getDistance3D(x[j]-dx,y[j]-dy,z[j]-dz ,   0 , 0 , 0 , radius,npType);
              dxm1 = getDistance3D(x[j]-2*dx,y[j]-2*dy,z[j]-2*dz ,   0 , 0 , 0 , radius,npType);

              ep2 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxp2, np.m_npBeadType[k] ));
              ep1 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxp1, np.m_npBeadType[k] ));
              em1 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxm1, np.m_npBeadType[k] ));
              em2 = static_cast<double>(potentials[ pdb.m_id[j]].Value(dxm1, np.m_npBeadType[k] ));


              dudx += ( -1 *ep2 + 8*ep1 - 8*em1 + em2  )/(12.0 * dr) ;



              //dudx += ( static_cast<double>(potentials[ pdb.m_id[j]].Value(dx1, np.m_npBeadType[k] )) - static_cast<double>(potentials[ pdb.m_id[j]].Value(dx0, np.m_npBeadType[k] )) )/dr ;
               }
            }
            return dudx;
}



//This function is used to relax the initial structure onto the surface of the NP. 
//It's essentially a very approximate molecular dynamics engine in a Brownian dynamics regime, so we don't even keep track of particle velocities

//unit convention: all distances are in NM for consistency, energies are KBT as they must interface directly with NP potentials
inline void RelaxPDB (const int size, double *x, double *y, double *z, const Config& config, const Potentials& potentials ,const PDB& pdb,const NP& np, double radius, int npType,  bool applyFinalShift) {
   //int numSteps = 10000;
   bool bSavePDB=false; //this is separate to the global one and used just for testing relaxation code


   int numFinalAverageSteps = 200; 
    double xFin[size];
    double yFin[size];
    double zFin[size]; 

    //store the force associated with timestep N  - only needed for the Langevin dynamics, which is not ready yet
    /*
    double xForce[size];
    double yForce[size];
    double zForce[size];
    double vx[size];
    double vy[size];
    double vz[size];
    */

    bool doCentralPull = false; //if false, pull is in -z only

    //partition the total number of steps: the initial pulling force is used for the first half and then equilibrium positions are sampled for the final 10% of steps 
    int numSteps = std::max(10, config.m_relaxSteps);
    int numPullSteps =  int(numSteps/2);
    numFinalAverageSteps = std::max(5, int(numSteps/10) ); 
   
    //fewer than 5 relaxation steps is interpreted as meaning "actually just do rigid dynamics"
    if(config.m_relaxSteps <= 5){
    numSteps = 0;
    numPullSteps = 0;
    numFinalAverageSteps = 0;
    }
    //int numPullSteps = std::min(10,  int( ceil(numSteps/10)) );
    double dr = 0.001; //step size for numerical gradients - we use central differencing so it ends up being half of this in either direction
    
    double initialZGradient = config.m_relaxZPull; //force in kbT per nm to apply  initially



    double kbTVal = 1.0; 
    double dtVal = 1e-11;

    //GJF alpha is equal to gamma*mass
    

    double nominalDiffusionCoeff = 1e5; //kbT/(gamma * mass) 
    double updateScaleSize = dtVal*nominalDiffusionCoeff/kbTVal; 
    bool addNoise  =   true;

    double noiseScale = 1.00 ;
    double noiseMagnitude = sqrt( 2 * noiseScale* nominalDiffusionCoeff * dtVal);

    //calculate parameters for GJF method
    /*
    double gjfMass = 10e-3;

    double gjfDampingRate = 1.0/(200.0*dtVal) ; 
    double gjfAlpha =  gjfDampingRate*gjfMass ; //kbTVal/nominalDiffusionCoeff ; //friction coefficient - this is equal to kbT/diffusion  


    noiseMagnitude = sqrt(2 * kbTVal * gjfAlpha *dtVal); //sdev of noise
    double gjfB = 1.0/(1.0 + gjfAlpha*dtVal / (2*gjfMass) );
    double gjfA = (1.0 - gjfAlpha*dtVal/(2*gjfMass))/  (1.0 + gjfAlpha*dtVal/(2*gjfMass)) ;
    //std::cout << "noise magnitude " << noiseMagnitude << " A:" << gjfA << " B: " << gjfB << "\n";
    //exit(1);

    */
    double maxDisplacement = 0.2; //set the maximum a bead can displace in a single step - this is mostly to stop explosions 
    //note that this implies updateScaleSize*force should ideally be less than 


    //for stability in the absence of noise, dt*nominalDiffusion =1e-5 is stable

    double proteinMaxDisplace = 1.0;
    double protMinZ = 0;
    for(int j=0; j<size; ++j){
        xFin[j] = x[j];
        yFin[j] = y[j];
        zFin[j] = z[j];

        if(z[j] < protMinZ){
        protMinZ = z[j];
        }
    }


    //writePDB("test_initial", "pdbouttest", size,  pdb, x, y, z);


    proteinMaxDisplace = 0.1 - protMinZ;

    

    //step 1: place the initial offset so that the biomolecule is just in contact with the NP
    




    std::vector< std::vector<int>  > nbNeighbours;
    int indexAtom = 0;
    double indexDist = 500;
       for(int j = 0; j < size; ++j){
          std::vector<int> iExclusionsNB;
           nbNeighbours.emplace_back(iExclusionsNB);
          double newDist = getDistance3D(x[j],y[j],z[j]+radius ,   0 , 0 , 0 , radius,npType);
          if(newDist < indexDist){
          indexAtom = j;
            indexDist = newDist;
          }
       }


    double x0init = x[indexAtom];
    double y0init = y[indexAtom];
    double z0init = z[indexAtom];
 int numNPs = np.m_npBeadType.size() ;
     if (config.m_sumNPPotentials == true){
      numNPs = 1;
        }


    double calcOffset = 0; 

            double currentOffset  = 0;  // -1.0 * (  *std::max_element(z, z + size) );
            double appliedOffset = 0.0; //store the offset applied in the z-direction relative to the COM
            bool shiftToSeparation = true; //this should only be set to false if you want to test for compatability with old-style UA
            bool scanFullProtein = config.m_enableFullScan; //if true (this will be settable via config file) then the start/stop parameters are set to ensure the full protein gets scanned. make sure there is a repulsive potential.
            bool foundShift = false;
            double beadDelta = 0.01;


            bool useCubeDist = false;
            int xfactor = 1; //used to define if x/y contribute to the distance
            int yfactor = 1;
            if(npType == 2 || npType == 4 || npType==5){
            xfactor = 0;
            }
            if(npType == 3){
            xfactor = 0;
            yfactor = 0;
            useCubeDist = true;
            }

            double rClosest = radius;
            if(config.m_zshiftToPlane==true){
            scanFullProtein = false;
            shiftToSeparation =  false;

            }



            //algorithm for finding the offset at closest approach
            
            double currentRCOM0 = 0;
            //double beadDelta = 0.4;
            //double rInner = radius + appliedOffset; is the closest approach for the RCC value . This means that setting appliedOffset = 0 corresponds to the COM of the protein lies on the NP radius
            //Next apply a further shift so that the protein can get closer for small NPs, long proteins or concave proteins for which the zmin plane results in a large distance from the NP
            //the algorithm: for each bead, compute the vertical distance necessary to bring that bead into contact with the NP (distSq < 0 implies no shift will achieve this)
            currentOffset = -radius; //set a default value such that RCC can take a value of zero
            if(shiftToSeparation == true){
            for(int i = 0; i < size; ++i){
             double distSq =  (radius + beadDelta)*(radius + beadDelta)   - (xfactor* x[i]*x[i] + yfactor*y[i]*y[i] );
             
             if(useCubeDist == true){
             //cubes are a special case 
             if(  std::abs( x[i]) > radius || std::abs( y[i] ) > radius){
             
             distSq = -5 ;
             }
             }
             
            double dzNeededI = 0;
             
             
            if(distSq >0){
            
            
             if(scanFullProtein == true){
              double dz1 = - std::sqrt(distSq) - z[i] ; 
              double dz2 =   std::sqrt(distSq) - z[i] ; 
              double newRCOM0 = 0;
              if(dz1 > 0){
               newRCOM0 = dz1;
              }
              else{
               if(dz2 > 0){
                newRCOM0 = dz2;
               }
            
              }
            
              //check to see if the bead is outside the NP radius anyway in which case we just skip it - for cubes using getDistance, noting that this returns 0 if the bead is inside the cube and allowing a small epsilon because floats
              if(   (xfactor*x[i]*x[i] + yfactor*y[i]*y[i] + z[i]*z[i] > radius*radius && useCubeDist == false) || ( useCubeDist == true && getDistance3D(x[i], y[i], z[i], 0.0,0.0,0.0, radius, 3)  > 0.01)    ){
               newRCOM0 = 0;
              }
            
              currentRCOM0 = std::max( currentRCOM0, newRCOM0);
              currentOffset = currentRCOM0 - radius;
            
            
            
            }
            
            
            else{
             dzNeededI = std::sqrt( distSq) -( z[i] + radius ) ; 
             currentOffset = std::max(currentOffset, dzNeededI);
            }
            
            foundShift = true;

            }
            
            
            }
            
            if(foundShift == false){
            //std::cout << "No distance could be found to bring the protein in contact \n";
            //std::cout << phi << " " << theta << "\n";
            currentOffset = -radius; //let the COMs overlap to handle ring-like proteins and small NPs
            }
            
            
            appliedOffset = currentOffset;
            }
            else{
            appliedOffset = -1.0 * (*std::min_element(z, z + size));   //classic mode: RCC at closest approach will put the lowest point of the protein on the NP radius-plane
            
            //ShiftZ(size, z);
            }

    appliedOffset += 0.0; //this puts the closest residue close to a minima
    //std::cout << " computed offset: " << appliedOffset << "\n";
    //std::cout <<" shifting all beads to +z by " << appliedOffset + radius << "\n"; 


    //std::cout << "index bead before at: " << x[indexAtom]<<","<<y[indexAtom]<<","<<z[indexAtom]<<": "<<  getDistance3D(x[indexAtom],y[indexAtom],z[indexAtom] ,   0 , 0 , 0 , radius,npType) << "\n";

    //std::cout << "before\n";
    for(int i = 0; i < size; ++i) {

        z[i] = z[i] + appliedOffset + radius;
         //std::cout << x[i] << "," << y[i] <<  "," << z[i] << "\n";

    }
       indexDist = 500;
       for(int j = 0; j < size; ++j){
          double newDist = getDistance3D(x[j],y[j],z[j] ,   0 , 0 , 0 , radius,npType);
           //std::cout << j << ":"<<z[j]<<" " << newDist << " ";
          if(newDist < indexDist){
          indexAtom = j;
            indexDist = newDist;
          }
       }
    x0init = x[indexAtom];
    y0init = y[indexAtom];
    z0init = z[indexAtom];

   double minEnergy = 200;
   double minEnergyDZI = 0;
   double initScanDZ = 0.005;
   double distFromInit = 1.0;
   int numDZIPoints = int( distFromInit/initScanDZ );
 
   for(int dzi =0; dzi < numDZIPoints; ++dzi){
       double totalEnergy = 0;
       for(int i = 0; i< size; ++i){
             
             for(int k=0; k<numNPs; ++k){
              if(config.m_sumNPPotentials == false){
                  double dist = getDistance3D(x[i],y[i],z[i]+dzi*initScanDZ ,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
                  totalEnergy += ( static_cast<double>(potentials[ pdb.m_id[i]].Value(dist, np.m_npBeadType[k] ))) ;
              }
              else{
                  double dist = getDistance3D(x[i],y[i],z[i]+dzi*initScanDZ ,   0 , 0 , 0 , radius,npType);
                  totalEnergy += ( static_cast<double>(potentials[ pdb.m_id[i]].Value(dist, np.m_npBeadType[k] ))) ;
              }
             }
       }
       //std::cout << dzi << " " << dzi*initScanDZ << " " << totalEnergy << "\n";
       if(totalEnergy < minEnergy){
           minEnergy = totalEnergy;
           minEnergyDZI = dzi;
       }
   }
   //exit(1);
   //std::cout <<"minima located at " << minEnergyDZI << "setting as initial offset \n";

   //we set it very slightly to the right of the discovered minima - that way we ensure the particle is coming in from infinity
   for(int i = 0; i< size; ++i){
     z[i] += (initScanDZ) * (1+minEnergyDZI) ; 
    }


   //step 1b -  optimise the initial  separation

    //writePDB("test_init", "pdbouttest", size,  pdb, x, y, z);

  
   //step 1c - rigid body optimisation
    //std::cout << "starting rigid optimisation \n";
     double initRigidEnergy =  GetEnergy(x, y,z, 0,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType);
     double lastInitRigid = initRigidEnergy + 1.0;
     double rigidEnergyEps = 1e-6;
     bool fixRigidCOMAxis = false;
    int numRigidSteps =   config.m_rigidSteps;
    double rigidDT = 5e-4; 

     bool rigidDone = false;

     //if all the particles lie on a plane or a line then the inertia matrix is non-invertible and we can't get accelerations from torque
     //so we add four virtual particles in a tetrahedron around the COM to prevent this. These particles are etremely lightweight  so they only add trivial damping
     double fakeTetraSize = 0.05; 
     double tetraPointsX[4] = { fakeTetraSize*sqrt(8.0/9.0), -1*fakeTetraSize* sqrt(2.0/9.0), -1*fakeTetraSize*sqrt(2.0/9.0), 0.0 };
     double tetraPointsY[4] = { fakeTetraSize*0, sqrt(2.0/3.0)*fakeTetraSize, -1*fakeTetraSize*sqrt(2.0/3.0), 0.0 };
     double tetraPointsZ[4] = { -1.0/3.0*fakeTetraSize, -1.0/3.0*fakeTetraSize, -1.0/3.0 * fakeTetraSize, fakeTetraSize*1 };

     double ixxf = 0;
     double iyyf = 0;
     double izzf = 0;
     double ixyf = 0;
     double ixzf = 0;
     double iyzf = 0;

      double fakemass = 0.01;
      for( int j=0; j< 4; ++j){
          double fakex = tetraPointsX[j];
          double fakey = tetraPointsY[j];
          double fakez = tetraPointsZ[j];
          ixxf += fakemass*(fakey*fakey + fakez*fakez);
          iyyf += fakemass*(fakex*fakex + fakez*fakez);
          izzf += fakemass*(fakey*fakey + fakex*fakex);
          ixyf -=  fakemass*(fakex*fakey);
          ixzf -=  fakemass*(fakex*fakez);
          iyzf -=  fakemass*(fakey*fakez);
     }



    for( int s=0; s< numRigidSteps; ++s){
      if(rigidDone == true){
      //std::cout << " rigid done \n ";
      break;
      }

      double taux = 0;
      double tauy = 0;
      double tauz = 0;

      double fcomx = 0;
      double fcomy = 0;
      double fcomz = 0;
      double comx = 0;
      double comy = 0;


      double comz = 0;
      double ixx = ixxf;
      double iyy = iyyf;
      double izz = izzf;
      double ixy = ixyf;
      double ixz = ixzf;
      double iyz = iyzf;


      for (int j=0; j< size; ++j){
      comx += x[j]/size;
      comy += y[j]/size;
      comz += z[j]/size;
  
      }






      for (int j=0; j< size; ++j){
      //comz += z[j]/size;
      ixx += (y[j]-comy)*(y[j]-comy) + (z[j]-comz)*(z[j]-comz);
      iyy += (x[j]-comx)*(x[j]-comx) + (z[j]-comz)*(z[j]-comz);
      izz += (y[j]-comy)*(y[j]-comy) + (x[j]-comx)*(x[j]-comx);
      ixy -=  (x[j]-comx)*(y[j]-comy);
      ixz -=  (x[j]-comx)*(z[j]-comz);
      iyz -=  (y[j]-comy)*(z[j]-comz);
      }


      for (int j=0; j< size; ++j){
        double fxj = -1.0*GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 0,  pdb,np, potentials, radius,npType);
        double fyj = -1.0*GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 1,  pdb,np, potentials, radius,npType);
        double fzj = -1.0*GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 2,  pdb,np, potentials, radius,npType);
        /*
        bool addRigidNoise = true;
        double rigidNoiseCoeff = sqrt(2*1e5* rigidDT);
        if(addRigidNoise == true){
            fxj += rigidNoiseCoeff*unitNormalDist(randomEngine) ;
            fyj += rigidNoiseCoeff*unitNormalDist(randomEngine) ;
            fzj += rigidNoiseCoeff*unitNormalDist(randomEngine) ;

        }   
        */
      

        fcomx += fxj;
        fcomy += fyj;
        fcomz += fzj;
        //tau is torque = r x F
        taux += (y[j]-comy) * fzj - (z[j]-comz)*fyj;
        tauy += (z[j]-comz) * fxj - (x[j]-comx)*fzj;
        tauz += (x[j]-comx) * fyj - (y[j]-comy)*fxj;
      }
     //std::cout << "Tau: " << taux << "," << tauy << "," << tauz << "\n";
     //std::cout << "z: " << comz <<   "FZCom: " << fcomz << "\n"; 
            //double inertiaEpsilon  = 1e-5;
     double inertiaDenom = (std::pow(ixz,2)*iyy - 2*ixy*ixz*iyz + std::pow(ixy,2)*izz + ixx*(std::pow(iyz,2) - iyy*izz));
     for (int j =0; j<size; ++j){
      //double axjold =  -tauz * (y[j]-comy)/izz + tauy * (z[j] - comz)/iyy;

      
      //double axjold = -tauz * (y[j]-comy)/izz + tauy * (z[j] - comz)/iyy;
      //double ayjold = tauz * (x[j]-comx)/izz - taux*(z[j]-comz)/ixx;
      //double azjold = -tauy * (x[j]-comx)/iyy + taux * (y[j]-comy)/ixx;
      
     double xrel = x[j] - comx;
    double yrel =y[j]-comy;
    double zrel = z[j]-comz;
       //double inertiaEpsilon  = 1e-5;
      //conversion from the resultant torque to per-particle components. 
      double axj =  ((-(ixz*iyy*taux) + ixy*iyz*taux + ixy*ixz*tauy - ixx*iyz*tauy - std::pow(ixy,2)*tauz + ixx*iyy*tauz)*yrel + 
     (ixy*izz*taux + std::pow(ixz,2)*tauy - ixx*izz*tauy + ixx*iyz*tauz - ixz*(iyz*taux + ixy*tauz))*zrel)/inertiaDenom;

     double ayj = ((ixz*iyy*taux - ixy*iyz*taux - ixy*ixz*tauy + ixx*iyz*tauy + std::pow(ixy,2)*tauz - ixx*iyy*tauz)*xrel + 
     (-(std::pow(iyz,2)*taux) + iyy*izz*taux + ixz*iyz*tauy - ixy*izz*tauy - ixz*iyy*tauz + ixy*iyz*tauz)*zrel)/inertiaDenom;


     double azj = ((ixz*iyz*taux - ixy*izz*taux - std::pow(ixz,2)*tauy + ixx*izz*tauy + ixy*ixz*tauz - ixx*iyz*tauz)*xrel + 
     (std::pow(iyz,2)*taux - iyy*izz*taux + ixy*izz*tauy + ixz*iyy*tauz - iyz*(ixz*tauy + ixy*tauz))*yrel)/inertiaDenom;

      azj += fcomz/size; //

      if(fixRigidCOMAxis==false){
      axj += fcomx/size;
      ayj += fcomy/size;
      }
      /*
      if(j==indexAtom){
          std::cout <<"acceleration: "<< axj << "  " << ayj << " " << azj << "\n";
          std::cout << "old acceleration  " << axjold  << " " << ayjold << " " << azjold<< "\n";

       std::cout << x[j] << "," << y[j] << "," << z[j] << " : " << inertiaDenom << "\n";
       }
       */

      
      //if( x[j]*x[j] > 1000 || y[j]*y[j] > 1000 || z[j]*z[j] > 1000 || axj*axj > 1000 || ayj*ayj>1000 || azj*azj>1000){
      /*
      if(j==indexAtom){
          //std::cout << "broke on atom " << j << "\n";
          std::cout << "before: " << x[j] << " " << y[j] << " " << z[j] << "\n";
          //std::cout << "old acceleration  " << axjold  << " " << ayjold << " " << azjold<< "\n";

           //for(int k =0; k<size; ++k){
           //std::cout << k<< ":" << x[k] << " " << y[k] << " " << z[k] << "\n";
           //}
          //std::cout << "x denom: " << (std::pow(ixz,2)*iyy - 2*ixy*ixz*iyz + std::pow(ixy,2)*izz + ixx*(std::pow(iyz,2) - iyy*izz)) <<"\n"; 
          //std::cout << "y denom: " << (std::pow(ixz,2)*iyy - 2*ixy*ixz*iyz + std::pow(ixy,2)*izz + ixx*(std::pow(iyz,2) - iyy*izz)) <<"\n"; 
          //std::cout << "z denom: " << (std::pow(ixz,2)*iyy - 2*ixy*ixz*iyz + std::pow(ixy,2)*izz + ixx*(std::pow(iyz,2) - iyy*izz)) <<"\n";

          //std::cout << ixx << " " << iyy << " " << izz << "\n";
          //std::cout << ixy << " " << ixz << " " << iyz << "\n";

          std::cout <<"acceleration: "<< axj << "  " << ayj << " " << azj << "\n";
          //std::cout << "old acceleration  " << axjold  << " " << ayjold << " " << azjold<< "\n";


          //exit(1);
      }
       
       */
      //for brownian motion: x' = x + dt*diffusion/kbT 
      

      x[j] += axj * rigidDT;
      y[j] += ayj * rigidDT;
      z[j] += azj * rigidDT;
      if(j==indexAtom){
          double currentTotalEnergy = GetEnergy(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType);
           double eDelta = lastInitRigid - currentTotalEnergy; //positive if energy is decreasing
          //std::cout << "updated to: " << x[j] << " " << y[j] << " " << z[j] << "current energy: " << currentTotalEnergy << " current delta: " <<  eDelta << "\n"; 
          //std::cout << "COMS: " << comx << " " << comy << " " << comz << "\n";
          if( eDelta < rigidEnergyEps  ){
          //std::cout << " probably converged, saving \n";
          rigidDone = true;
          }
          lastInitRigid  = std::min(lastInitRigid, currentTotalEnergy);
      }
     }


    }
    double postRigidEnergy = GetEnergy(x, y,z, 0,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType) ; 
    //std::cout << "Rigid mechanics: Initial energy: " << initRigidEnergy << " updated to " << GetEnergy(x, y,z, 0,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType) << "\n"; 
    //std::cout << "index bead now at: " << x[indexAtom]<<","<<y[indexAtom]<<","<<z[indexAtom]<<": "<<  getDistance3D(x[indexAtom],y[indexAtom],z[indexAtom] ,   0 , 0 , 0 , radius,npType) << "\n";
    //writePDB("test_post_rigid", "pdbouttest", size,  pdb, x, y, z);
    //exit(1);
   
    //step 2: identify "bonded pairs" of residues  - these are ones we treat as if they are harmonically bonded
    double bondCutoff = config.m_bondCutoffNM; 
    //double bondCutoffSq = bondCutoff*bondCutoff;
    int numBonds = static_cast<int>(pdb.m_bondSet.size());

   //step 3: overdamped relaxation - no velocity just position updates towards the local minima
   //the basic equation of interest here is:
   //xprime = x + deltat * ( sqrt(2 * kbT/gamma*mass) - 1/(gamma*mass) dU/dx )
   // units:
   // nanometers = nanometers + seconds*(  sqrt( 2 * kbT * /(kg*seconds)  )- 1/(seconds*kg)  * kbT * energy[kbT]/ nanometers  )


   //int numSteps = 1000;
   //double dt = 1e-12;
   //double frictionCoeff = 1e-12;
   //double updateScaleSize = 1e-5 ;  //  dt/frictionCoeff; 

   std::vector<double> xUpdates(size);
   std::vector<double> yUpdates(size);
   std::vector<double> zUpdates(size);
   double wcaEpsilon = 0.25; //kBT
   double wcaSigma = 0.4; //nm
       //writePDB("test_pre", "pdbouttest", size,  pdb, x, y, z);

   int numAverageStepsDone = 0;




    //set initial forces to zero and velocities to thermal noise

    /*for(int i=0; i<size; i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    xForce[i] = 0.0;
    yForce[i] = 0.0;
    zForce[i] = 0.0;
    }
   */
   for(int s = 0; s < numSteps ; ++s){
       //reset all updates for initialisation
       //xUpdates[:] = 0;
       //yUpdates[:] = 0;
       //zUpdates[:] = 0;
      //std::cout << "cleaning updates \n"; 
                  std::fill( xUpdates.begin(), xUpdates.end(),0.0 );
                  std::fill( yUpdates.begin(), yUpdates.end(),0.0 );
                  std::fill( zUpdates.begin(), zUpdates.end(),0.0 );

            //std::cout << "NP forces \n";


       //phase 1: get updates due to the NP-biomolecule potential(s)
       for(int j = 0; j < size; ++j){
          
              double dudx = GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 0,  pdb,np, potentials, radius,npType);
              double dudy = GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 1,  pdb,np, potentials, radius,npType);
              double dudz = GetPotGrad4(x, y,z, j,  dr, config.m_sumNPPotentials,  numNPs, 2,  pdb,np, potentials, radius,npType);
               xUpdates[j] -=  dudx; //potGradX(i,k); 
               yUpdates[j] -=  dudy; //potGradY(i,k);
               zUpdates[j] -=  dudz; //potGradZ(i,k);
           }
       
        //std::cout << "z updates for index atom after NP" << zUpdates[indexAtom] << "\n";


       //phase 2: get updates due to bonded beads
       for( int k=0; k<numBonds; ++k){
          
           int i = pdb.m_bondSet[k].i;
           int j = pdb.m_bondSet[k].j;
           //std::cout << "updating bond " << i << " to " << j << "\n";
           double ijDist =  sqrt( pow(x[i]-x[j],2 )+pow(y[i]-y[j],2 )+pow(z[i]-z[j],2 ) ); 
           double bondLength = pdb.m_bondSet[k].length; 
           double bondMag = pdb.m_bondSet[k].bondk;
           ijDist = std::max( 0.001, ijDist) ; 
           double xForceB = -bondMag*( x[i] - x[j] )*( ijDist - bondLength)/ijDist;
           double yForceB = -bondMag*( y[i] - y[j] )*( ijDist - bondLength)/ijDist;
           double zForceB = -bondMag*( z[i] - z[j] )*( ijDist - bondLength)/ijDist;
           //std::cout << " bond " << k << "eq. length: " << bondLength << " current length  " << ijDist << "bond constant " << bondMag << "\n";
           //std::cout << x[i] << " to " << x[j] <<     "," <<     y[i] << " to " << y[j] << "," <<z[i] << " to " << z[j] <<           "\n";
           //std::cout << "forces on atom "<< i <<  xForce << "," << yForce << "," << zForce << "\n";

          
           xUpdates[i] += xForceB;
           yUpdates[i] += yForceB;
           zUpdates[i] += zForceB;

           xUpdates[j] -= xForceB;
           yUpdates[j] -= yForceB;
           zUpdates[j] -= zForceB;
           

       }

        //std::cout << "z updates for index atom after bond "  << zUpdates[indexAtom] << "\n";


      //phase 3: nonbonded terms - this is just to prevent overlap

      //first update the neighbour list - once a bead comes within 2nm of another bead, we permanently switch on the NB term
      int neighbourUpdateFreq = 10;
       if( s % neighbourUpdateFreq == 0){
           for(int i = 0; i< size; ++i){
               for(int j = i+1; j < size; ++j){
                   if( std::find( pdb.m_nbExclusions[i].begin(),pdb.m_nbExclusions[i].end(), j) == pdb.m_nbExclusions[i].end() ){
                   //particle is not in the permanent exclusions
                        if( std::find( nbNeighbours[i].begin(),nbNeighbours[i].end(), j) == nbNeighbours[i].end() ){
                        //particle is not already a neighbour
                            double resDist = sqrt( pow(x[i]-x[j],2 )+pow(y[i]-y[j],2 )+pow(z[i]-z[j],2 ) );
                            if(resDist < 1.4){
                                nbNeighbours[i].emplace_back(j);
                                //std::cout << "step: " << s << ":" << i << "," << j << " have become neighbours \n";
                            }
                        }
                   }
               }
           }
       }

       for(int i = 0; i< size; ++i){
        for(auto & j: nbNeighbours[i] ){ //only loop over neighbours
        //std::cout << "computing NB for " << i << "," << j << "\n";
        //for(int j = i+1; j < size; ++j){
                         //std::cout << s <<  " computing nonbond " << i << " to " << j << "\n";

            double resDist = sqrt( pow(x[i]-x[j],2 )+pow(y[i]-y[j],2 )+pow(z[i]-z[j],2 ) );
             //std::cout << "distance is: " << resDist << "\n";
            //double wcaCutoff = 
           //nbExclusions[i]
              //if find returns the end then it means j is not found in the exclusion list for atom i
            //if( std::find( pdb.m_nbExclusions[i].begin(),pdb.m_nbExclusions[i].end(), j) == pdb.m_nbExclusions[i].end() &&     resDist < 1.2246 * wcaSigma){
       
             //std::cout << s << " applying nonbond " << i << " to " << j << "\n";
              //std::cout << " getting radii for i,j " << pdb.m_resTag[i] <<  ":" << config.AminoAcidRadius(pdb.m_resTag[i]) << " " <<pdb.m_resTag[j] << ":" << config.AminoAcidRadius(pdb.m_resTag[j]) <<    "\n";
               wcaSigma =  config.AminoAcidRadius(pdb.m_resTag[i]) +  config.AminoAcidRadius(pdb.m_resTag[j]) ;
                if(resDist < 1.4){

                //double wcaMag = 4 * wcaEpsilon * pow(wcaSigma,12);
                 //double wcaForceMag  = wcaMag * 12.0/pow(resDist*resDist,7);
                 double wcaMag = 4 * wcaEpsilon * pow(wcaSigma,6);
                 double wcaForceMag = wcaMag * 6.0 / pow(resDist*resDist,4);
                 double xForceLJ = (x[i] - x[j]) *wcaForceMag;
                 double yForceLJ = (y[i] - y[j]) *wcaForceMag;
                 double zForceLJ = (z[i] - z[j]) *wcaForceMag;

 
                 xUpdates[i] += xForceLJ;
                 yUpdates[i] += yForceLJ;
                 zUpdates[i] += zForceLJ;
                 xUpdates[j] -= xForceLJ;
                 yUpdates[j] -= yForceLJ;
                 zUpdates[j] -= zForceLJ;

            }
         }
       }


    







   //phase 4: apply an extra pull for the initial few steps (now merged into stage 5) 
          //std::cout << "z updates for index atom after NP+bond+nonbond " << zUpdates[indexAtom] << "\n";

   //phase 5: apply updates
   //bool addNoise = true;
   //if(addNoise == true and s > numSteps){
   // addNoise = false; //disable the noise once we reach the final stage
   //}

   for(int i = 0; i< size; ++i){
        //std::cout << s<< " " << i << "applying updates " << xUpdates[i]<<","<<yUpdates[i]<<","<<zUpdates[i]<< "\n";

            if( s < numPullSteps){
             if(doCentralPull == true){
              double centralDist = sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i] ) ;
              if(centralDist > 1e-5){ 
               xUpdates[i] -= initialZGradient*x[i]/centralDist;
               yUpdates[i] -= initialZGradient*y[i]/centralDist;
               zUpdates[i] -= initialZGradient*z[i]/centralDist;
              }


              }
              else{
             zUpdates[i] -= initialZGradient; //apply a nominal force of 1 kbT/nm in the negative z direction
             }

            }


        
        //bool addNoise = true;
         //add thermal noise
        double xNoiseNP1=0;
        double yNoiseNP1=0;
        double zNoiseNP1=0;

        if(addNoise==true && s < numPullSteps ){
        xNoiseNP1 = noiseMagnitude*unitNormalDist(randomEngine);
        yNoiseNP1 = noiseMagnitude*unitNormalDist(randomEngine);
        zNoiseNP1 = noiseMagnitude*unitNormalDist(randomEngine);

        
        }
         /*
         if(i==indexAtom){
        std::cout << "Moving bead " << i <<  "("<<x[i]<<","<<y[i]<<","<<z[i]<<") : " << xUpdates[i] <<"," << yUpdates[i] << "," << zUpdates[i] ;
           if( s < numPullSteps){
             std::cout << " (pull on) \n";
           }
           else{
             std::cout  << "\n"; 
            }
        } */

       //GJF algorithm for  velocity verlet



        //if(i == indexAtom){
        //std::cout << s << "before x: " << x[i] << " " << vx[i] <<" " <<  xForce[i] << "\n";
        //std::cout << s << "y: " << y[i] << " " << vy[i] << yForce[i] << "\n";
        //std::cout << s << "z: " << z[i] << " " << vz[i] << " " <<  zForce[i] << "\n";
        //std::cout << gjfA << "\n";
        //std::cout <<"vel size: " <<  gjfB << " " << gjfB * dtVal * vx[i] << "\n";
        //std::cout << xNoiseNP1 << "\n";
        //std::cout << "noise contribution to position: " << gjfB * dtVal/(2*gjfMass) *xNoiseNP1 << "\n"; 
        //}

       /*if(s>0){
        x[i] = x[i] + gjfB*dtVal * vx[i] + (gjfB*dtVal*dtVal/(2*gjfMass))*xForce[i] + (gjfB * dtVal/(2*gjfMass) )*xNoiseNP1;
        vx[i] = gjfA * vx[i] + dtVal/(2*gjfMass)*(gjfA * xForce[i] + xUpdates[i]) + gjfB/gjfMass * xNoiseNP1 ; 
        //xForce[i] = xUpdates[i];   


        y[i] = y[i] + gjfB*dtVal * vy[i] + (gjfB*dtVal*dtVal/(2*gjfMass))*yForce[i] + (gjfB * dtVal/(2*gjfMass)) *yNoiseNP1;
        vy[i] = gjfA * vy[i] + dtVal/(2*gjfMass)*(gjfA * yForce[i] + yUpdates[i]) + gjfB/gjfMass * yNoiseNP1 ;
        //yForce[i] = yUpdates[i];


        z[i] = z[i] + gjfB*dtVal * vz[i] + (gjfB*dtVal*dtVal/(2*gjfMass))*zForce[i] + (gjfB * dtVal/(2*gjfMass)) *zNoiseNP1;
        vz[i] = gjfA * vz[i] + dtVal/(2*gjfMass)*(gjfA * zForce[i] + zUpdates[i]) + gjfB/gjfMass * zNoiseNP1 ;
        }
        xForce[i] = xUpdates[i];
        yForce[i] = yUpdates[i];
        zForce[i] = zUpdates[i];


        if(i == indexAtom){
        //std::cout << s << "x: " << x[i] << " " << vx[i] <<" " <<  xForce[i] << "\n";
        //std::cout << s << "y: " << y[i] << " " << vy[i] << yForce[i] << "\n";
        //std::cout << s << "z: " << z[i] << " " << vz[i] << " " <<  zForce[i] << "\n";


        }

        */


        //basic old method - note that the scaling is no longer correct
        x[i] += std::min( std::max(-maxDisplacement, updateScaleSize*xUpdates[i] + xNoiseNP1), maxDisplacement)  ;
        y[i] += std::min( std::max(-maxDisplacement, updateScaleSize*yUpdates[i] + yNoiseNP1), maxDisplacement);
        z[i] += std::min( std::max(-maxDisplacement, updateScaleSize*zUpdates[i] + zNoiseNP1), maxDisplacement); 

        if(i == indexAtom){
        //std::cout << s << "x: " << x[i] <<  "\n";
        //std::cout << s << "y: " << y[i] << "\n";
        //std::cout << s << "z: " << z[i] << " E: " <<  GetEnergy(x, y,z, 0,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType) <<"\n" ;


        }



         /*
        if(i==indexAtom){
        std::cout << std::min( std::max(-maxDisplacement, zUpdates[i]), maxDisplacement) << "\n";
        }
         */

        //std::cout << x[i] << "," << y[i] << "," << z[i] << "\n"; 
     } 
    

   //if we're in the final stage, compute average positions
     if(s > numSteps -numFinalAverageSteps){
      for(int i =0; i<size; ++i){
        xFin[i] += x[i];
        yFin[i] += y[i];
        zFin[i] += z[i];
      } 
        numAverageStepsDone += 1;
        //std::cout << "current: " << x[indexAtom] << "," << y[indexAtom] << "," << z[indexAtom] << "\n";
        //std::cout << xFin[indexAtom]/numAverageStepsDone << "," <<  yFin[indexAtom]/numAverageStepsDone  << "," << zFin[indexAtom]/numAverageStepsDone  << "\n";

 
     }


   }
  // std::cout << "completed relaxation \n";
    //std::cout << "writing pdb \n";
    //save the protein out in the final position

    /*
    if(bSavePDB == true){
    writePDB("test_post", "pdbouttest", size,  pdb, x, y, z);
    }
    */
    //std::cout <<"pdb done, exiting \n";
    //std::cout <<  x0init<<","<<y0init<<","<<z0init<<"\n";
    //exit(1);
    //writePDB("test_post_md_preaverage", "pdbouttest", size,  pdb, x, y, z);

    //replace coordinates by the average over the last set of steps to reduce thermal fluctuations
    
    bool replaceAverage = false;
    if(numAverageStepsDone > 0 && replaceAverage==true){
           for(int i =0; i<size; ++i){
               x[i] = xFin[i]/numAverageStepsDone;
               y[i] = yFin[i]/numAverageStepsDone;
               z[i] = zFin[i]/numAverageStepsDone;

              //yFin[i] += y[i];
              //zFin[i] += z[i];
      }

    }
   //end the stepping - final post-processing
   //finally shift the protein back to its original z location

    //std::cout << "Post relaxation: Initial energy: " << initRigidEnergy <<  "  Rigid optimisation: " << postRigidEnergy << " MD updated to " << GetEnergy(x, y,z, 0,  dr, config.m_sumNPPotentials,  numNPs, size,  pdb,np, potentials, radius,npType) << "\n";


    //writePDB("test_post_md", "pdbouttest", size,  pdb, x, y, z);
    //exit(1);

    double newZCom = 0;
    //std::cout << "after\n";
    for(int i = 0; i < size; ++i) {
      //   std::cout << x[i] << "," << y[i] <<  "," << z[i] << "\n";


        newZCom += z[i]/size;
    }
    //exit(0);
    if(applyFinalShift == true){
    for(int i = 0; i < size; ++i) {
        z[i] = z[i] - newZCom;

        //std::cout << x[i] << "," << y[i] <<  "," << z[i] << "\n";

    }
    }

    //writePDB("test_post_backshift", "pdbouttest", size,  pdb, x, y, z);
    //std::cout <<"pdb done, exiting \n";

     //exit(1);
    /*
    std::cout << "shifted: " << x0init << " to " << x[indexAtom] << "\n";
    std::cout << "shifted: " << y0init << " to " << y[indexAtom] << "\n";
    std::cout << "shifted: " << z0init << " to " << z[indexAtom] << "\n";
    */

}

void SaveSpecificOrientation(double phi, double theta, const int size,  const Config& config, const Potentials& potentials ,const PDB& pdb,
const NP& np, double radius, int npType, double cylinderAngle, const std::string& pdbname, const std::string& outputdirectory, const std::string& npName, double ccdAtMin) {

    double x[size];
    double y[size];
    double z[size];
    double phi_adjusted   = -1.0 * (phi  + angle_delta/2);
    double theta_adjusted = M_PI - (theta + angle_delta/2);

     Rotate3(size, phi_adjusted, theta_adjusted,cylinderAngle, pdb.m_x, pdb.m_y, pdb.m_z, x, y, z);

     if(config.m_relaxPDB == true){
     writePDB(pdbname+"-prelax", outputdirectory+"/"+npName+"/opt_pdbs",  size,  pdb, x, y, z);

     RelaxPDB (size, x, y, z, config, potentials ,pdb, np, radius, npType,  false) ; 
     }
     else{
       for(int i=0; i<size; ++i){ 
          z[i] += ccdAtMin;  
     }  
     } 
     writePDB(pdbname+"-relax", outputdirectory+"/"+npName+"/opt_pdbs",  size,  pdb, x, y, z);

}


void AdsorptionEnergies(const PDB& pdb,const NP& np, const Config& config, const Potentials& potentials, const double radius, const double outerRadius,const int angle_offset, const int n_angles, double *adsorption_energy, double *adsorption_error, double *mfpt_val, double *mfpt_err, double *minloc_val, double *minloc_err,  double *numcontacts_val,   double *numcontacts_err,  double *protein_offset_val, double *protein_offset_err,  const std::string& pdbname, const std::string& outputdirectory, const std::string& npName,      int npType = 1, double imaginary_radius = -1, int calculateMFPT = 0, double cylinderAngleDeg=0,double zeta =0,int savePotentials = 0,double temperature=300.0, double overlapPenalty=0.0) { 

    // Decleare all variables at the begining . Integration runs from the start at the largest value of r and proceeds inwards to the stop value
    const int               size            = pdb.m_id.size();
    //const double            stop            = imaginary_radius < 0 ? radius + gds : imaginary_radius + gds; //gds = 0 and imaginary radius is unused so this should always just return radius (inner bounding radius)
    

    int                     i;
    int                     j;
    int			k;
    double                  energy;
    double                  theta;
    double                  phi;
    double                  theta_adjusted;
    double                  phi_adjusted;
    double                  init_energy;
    double                  distance;
    double                  ssd;
    double                  rcc; //R com-com distance
     
    double stop = radius;
    double x[size];
    double y[size];
    double z[size];
    
    double total_energy[steps];
    double SSD[steps];
    double sample_energy[samples];
    double sample_mfpt[samples];
    double sample_minloc[samples];
    double sample_proteinoffset[samples];
    double sample_numcontacts[samples];
    double closestSCD[size]; //closest allowed NP_nominal_surface to aa_centre distance allowed for each AA
   //NP component properties are looked up via e.g.: np.m_radius[m_npBeadType[j]]
   //where j is the index of that bead
   //get the number of NP components
   int numNPBeads = np.m_npBeadType.size() ;
   //decide how many NP beads to sum over. By default, only the isotropic-averaged NP is used unless the beads are explicitly treated separately, in which case we sum over all beads
               int numNPBeadSummation = 1;
            if(config.m_sumNPPotentials == false){
            numNPBeadSummation = numNPBeads;
            }
     // if(angle_offset == 0){
   //std::cout << "thread 0 starting "  << angle_offset << "\n" << std::flush;
    //}
             int maxFlexPoints = ceil( 1.0/config.m_flexResolution);
            std::vector<double> flexEnergyTerms(maxFlexPoints*2 + 1 );
            std::vector<double> flexEnergyDenomTerms(maxFlexPoints*2 + 1 );




    // Iterate over angles
    for (int angle = 0; angle < n_angles; ++angle) {
        
        phi   = ((angle_offset + angle) % ncols) * angle_delta;
        theta = ((angle_offset + angle) / ncols) * angle_delta;
           double cylinderAngle = cylinderAngleDeg * M_PI/180.0;
           //std::cout << "starting: " << phi << " " << theta << "\n";

	    // Sample a angle multiple times.The sample = samples run is not included in averaging but instead used to generate output potentials  
        for (int sample = 0; sample < samples+1 ; ++sample) {
  //         double cylinderAngle = 90 * M_PI/180.0;
             //std::cout << sample << "/" << samples << "\n";
    	    // Rotation step - random sampling is enabled if the number of samples to generate is > 1. If not, it just calculates it at the bin center

            if(sample <  samples && samples > 1){
    	    phi_adjusted   = -1.0 * (phi +  angle_delta*random_angle_offset(randomEngine));
    	    theta_adjusted = M_PI - (theta +  angle_delta*random_angle_offset(randomEngine));
             }
             else{
            phi_adjusted   = -1.0 * (phi  + angle_delta/2);
            theta_adjusted = M_PI - (theta + angle_delta/2);
             //std::cout << "final sample, fixing angle to midpoint" << sample << " " <<  samples << "\n";
              }

            Rotate3(size, phi_adjusted, theta_adjusted,cylinderAngle, pdb.m_x, pdb.m_y, pdb.m_z, x, y, z);
            
            if(config.m_pdbJitterMag > 0.001){
            JitterPDB(size, config.m_pdbJitterMag, x, y, z);
             }
            
            //Prepare the initial offset of the protein.
            //There are three modes for finding the minimum value of R
            //1) Classic mode ( config.m_zshiftToPlane == true) - the closest approach is defined by the plane of the lowest z-coordinate. this gives odd results for elongated proteins.
            //2) Regular mode ( config.m_enableFullScan = false, m_zshiftToPlane == false) - put the protein at infinity, move it inwards until first contact (bead centre to NP radius), the R value at which this occurs is RMin = stop.
            //3) Full-scan mode (config.m_enableFullScan = true, m_zshiftToPlane == false) - put the protein's COM on the surface of the NP, move outwards until no overlap. 
            
            //For each of these, we must find a suitable closest approach RStop = NPRadius + appliedOffset to shift the protein COM to the correct position for the above case
            //classic mode: RStop = NPRadius - min(z) -> appliedOffset = -min(z)
            //regular mode: 
            //full-scan mode: 
            
            double finalRDelta = delta;
            double currentOffset  = 0;  // -1.0 * (  *std::max_element(z, z + size) );
            double appliedOffset = 0.0; //store the offset applied in the z-direction relative to the COM
            bool shiftToSeparation = true; //this should only be set to false if you want to test for compatability with old-style UA
            bool scanFullProtein = config.m_enableFullScan; //if true (this will be settable via config file) then the start/stop parameters are set to ensure the full protein gets scanned. make sure there is a repulsive potential.
            bool foundShift = false;
            double beadDelta = 0.01;
            
       
            bool useCubeDist = false;
            int xfactor = 1; //used to define if x/y contribute to the distance
            int yfactor = 1;
            if(npType == 2 || npType == 4 || npType==5){
            xfactor = 0;
            }
            if(npType == 3){
            xfactor = 0;
            yfactor = 0;
            useCubeDist = true;
            }

            double rClosest = radius;
            if(config.m_zshiftToPlane==true){
            scanFullProtein = false;
            shiftToSeparation =  false;

            }
      
      



            bool bRelaxProtein = config.m_relaxPDB;
            if(bRelaxProtein == true){
          //RelaxPDB (const int size, double *x, double *y, double *z, const Config& config, const Potentials& potentials ,const PDB& pdb,const NP& np)
                RelaxPDB( size, x, y, z, config, potentials, pdb, np, radius,npType,true); //relaxes the protein onto the surface of the NP
                //std::cout << "relaxation done, returning to main UA \n"; 
            }
            


            
            double currentRCOM0 = 0;
 
            //double rInner = radius + appliedOffset; is the closest approach for the RCC value . This means that setting appliedOffset = 0 corresponds to the COM of the protein lies on the NP radius
            //Next apply a further shift so that the protein can get closer for small NPs, long proteins or concave proteins for which the zmin plane results in a large distance from the NP
            //the algorithm: for each bead, compute the vertical distance necessary to bring that bead into contact with the NP (distSq < 0 implies no shift will achieve this)
            currentOffset = -radius; //set a default value such that RCC can take a value of zero

            if(shiftToSeparation == true){
            for( i = 0; i < size; ++i){
             double distSq =  (radius + beadDelta)*(radius + beadDelta)   - (xfactor* x[i]*x[i] + yfactor*y[i]*y[i] );
             
             if(useCubeDist == true){
             //cubes are a special case 
             if(  std::abs( x[i]) > radius || std::abs( y[i] ) > radius){
             
             distSq = -5 ;
             }
             }
             
            double dzNeededI = 0;
             
             
            if(distSq >0){
            
            
             if(scanFullProtein == true){
              double dz1 = - std::sqrt(distSq) - z[i] ; 
              double dz2 =   std::sqrt(distSq) - z[i] ; 
              double newRCOM0 = 0;
              if(dz1 > 0){
               newRCOM0 = dz1;
              }
              else{
               if(dz2 > 0){
                newRCOM0 = dz2;
               }
            
              }
            
              //check to see if the bead is outside the NP radius anyway in which case we just skip it - for cubes using getDistance, noting that this returns 0 if the bead is inside the cube and allowing a small epsilon because floats
              if(   (xfactor*x[i]*x[i] + yfactor*y[i]*y[i] + z[i]*z[i] > radius*radius && useCubeDist == false) || ( useCubeDist == true && getDistance3D(x[i], y[i], z[i], 0.0,0.0,0.0, radius, 3)  > 0.01)    ){
               newRCOM0 = 0;
              }
            
              currentRCOM0 = std::max( currentRCOM0, newRCOM0);
              currentOffset = currentRCOM0 - radius;
            
            
            
            }
            
            
            else{
             dzNeededI = std::sqrt( distSq) -( z[i] + radius ) ; 
             currentOffset = std::max(currentOffset, dzNeededI);
            }
            
            foundShift = true;

            }
            
            
            }
            
            if(foundShift == false){
            //std::cout << "No distance could be found to bring the protein in contact \n";
            //std::cout << phi << " " << theta << "\n";
            currentOffset = -radius; //let the COMs overlap to handle ring-like proteins and small NPs
            }
            
            
            appliedOffset = currentOffset;
            }
            else{
            appliedOffset = -1.0 * (*std::min_element(z, z + size));   //classic mode: RCC at closest approach will put the lowest point of the protein on the NP radius-plane
            
            //ShiftZ(size, z);
            }

            //std::cout << "Sample " << sample << " of " << samples  << " set offset to: " <<appliedOffset << "\n";
            //finalRDelta = delta - (*std::min_element(z,z+size)) ; 
            //Find the closest approach
            
            //the outermost point is previously that which the minimum point of the protein is at a distance delta=2 nm from the outer radius-plane of the NP 
            //we keep this for consistency: the RCC value is equal to the outerRadius of the NP + delta  + the protein's surface to COM width
            finalRDelta = delta  -1.0 * (*std::min_element(z,z+size));
            
            
          



           //appliedOffset = currentRCOM0 - radius;
           double rInner = radius + appliedOffset;
           double rOuter = outerRadius + finalRDelta;
        //  std::cout << rInner << " " << rOuter << "\n";
            //this overlap code is augmented by the WCA potential and is needed only if the NP potential is summed over initially, with the goal of identifying RCCs which cause overlap.
            
            for( i = 0; i< size; ++i){ //loop over AA beads
             double closestAllowedSSD = 0;
            
             closestSCD[i] = 0;
             
                  if(config.m_sumNPPotentials == true){
             
             if(  overlapPenalty > 0){
              double aaBeadRadius = config.m_aminoAcidRadii[ pdb.m_id[i] ] * config.m_overlapRadiusFactor ;
              //std::cout << pdb.m_id[i] <<  " " << config.m_aminoAcidRadii[ pdb.m_id[i] ]   <<"\n";
                   for( j = 0; j < numNPBeads  ; ++j){
                       //std::cout << i << ":" << j << "\n";
                       double npBeadRadius = np.m_radius[np.m_npBeadType[j]] * config.m_overlapRadiusFactor;

                       double distTerm1 = pow( aaBeadRadius + npBeadRadius, 2) - pow( x[i] - np.m_x[j],2 ) - pow( y[i] - np.m_y[j],2 );
                      //std::cout << distTerm1 << "\n";
                       if(distTerm1 > 0){
                                  
                            double ssdAtTouching1 = np.m_z[j] -aaBeadRadius  - z[i] - radius +  sqrt( distTerm1)  ;
                            //double ssdAtTouching2 = np.m_z[j] -aaBeadRadius  - z[i] - radius - sqrt( distTerm1 )  ;
                             
                            if( ssdAtTouching1 > closestAllowedSSD){ 
                            closestAllowedSSD = ssdAtTouching1;

                            }

                       }

                    }


               if(  closestAllowedSSD + aaBeadRadius > closestSCD[i]  ){
                     closestSCD[i] = closestAllowedSSD + aaBeadRadius  ; //correct for the fact that here SSD is the surface-surface distance while everywhere else ssd means surface-centre distance because David 
               }

            }
           }
           
           }
           
           
           
           
           //once we've established the protein offset, define the "start" location, i.e. the outermost r value
           double            start           =  rOuter;  //stop + delta;
           double            stop            = rInner;
           double actualDZ = (start - stop)/(steps-1);
           
              //x,y no longer squared
            
            init_energy = 0.0;



            for (i = 0; i < size; ++i) {
            

            
            
            for(k = 0; k<numNPBeadSummation; ++k) {
            
              if(config.m_sumNPPotentials == false){
              distance = getDistance3D(x[i],y[i],z[i] + start,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              }
              else{
              distance = getDistance3D(x[i],y[i],z[i] + start,  0 , 0 , 0 , radius,npType);
              }
                //std::cout << i << " " << k << " " << distance << "\n"; 
                double energyAtDist = pdb.m_occupancy[i] *  static_cast<double>(potentials[pdb.m_id[i]].Value(distance, np.m_npBeadType[k] ))  ;
                if(energyAtDist > 1){
                //std::cout << "large shift (" << energyAtDist << ") for distance " << distance << "\n";
                }
                //std::cout << pdb.m_occupancy[i] << "\n";
                init_energy +=  pdb.m_occupancy[i] *  static_cast<double>(potentials[pdb.m_id[i]].Value(distance, np.m_npBeadType[k] ));



              //find the distance of closest approach



         }             
            }
 



//std::cout << "total shift/num. amino acid: " <<  init_energy/size << "\n";
            // build the distance-potential array
              // 
            double minLoc = start;
            double minEnergy = 0;
            int numContactsAtMin = 0;
            int numContactsAtStep = 0;
            int resHasContacted = 0;
            double flexEnergyDenom = 0;
            double flexEnergyTot = 0;
               
                            int flexMethod = config.m_flexMethod;
            // int maxFlexPoints = 2;
            //if( config.m_flexMethod==4){
            //std::vector<double> flexEnergyTerms(maxFlexPoints*2 + 1 ); 
            //std::vector<double> flexEnergyDenomTerms(maxFlexPoints*2 + 1 );

            //}
            for (i = 0; i < steps; ++i) {
                rcc     = start - i * actualDZ; //was SSD - renamed because it didn't actually fit that description and this was making code maintenance difficult.
               // std::cout << rcc << "\n";
                
                energy  = 0;
                flexEnergyTot = 0;
                flexEnergyDenom = 0;
                 numContactsAtStep = 0;
                //loop over each residue. loop variables: i = R index, j = residue index. The ssd value itself is given by the radius of the NP plus an offset distance of i*dz


                for (j = 0; j < size; ++j) {
               flexMethod = config.m_flexMethod; //reset the flex method in case we previously dealt with an extremely well-defined residue
               resHasContacted = 0;
 
               double appliedOverlapPenalty = 0.0;


                //int flexMethod = config.m_flexMethod;
                 //flex methods: these are ways to allow a small amount of residue flexibility
                 //0: classic UA, all beads are fixed at their given coordinates
                 //1: Gaussian smoothing 
                 //2: minimum in interval with the penalty due to displacement
                 //3: local smoothing of potentials by free-energy
                 //4: free-energy averaging over noise co-ordinate

 
                 double rmsd = pdb.m_rmsd[j] ;
                 double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                 int numSDs = config.m_flexNumSDev;
                 double flexScanRes = config.m_flexResolution; //step size for flexibility
                 int numPoints =  ceil( numSDs*rmsd/flexScanRes) ;
                 //double totalFlexRange = 2*numPoints*flexRes;
                 numPoints = std::max( 1, numPoints) ;
                 numPoints = std::min(numPoints, maxFlexPoints); //cap how far the flexibility can go 
                 //set up how many points to loop over for this particular central location
                 if(rmsd < config.m_flexResolution){
                 flexMethod = 0; //for very well-defined residues, just return the central value
                 }
                 if(flexMethod == 0){
                 numPoints = 0;
                 }
                 if(flexMethod == 4){
                 //energyTerms = vector(float, 2*numPoints+1) ;
                  std::fill( flexEnergyTerms.begin(), flexEnergyTerms.end(),0.0 );
                  std::fill( flexEnergyDenomTerms.begin(), flexEnergyDenomTerms.end(),0.0 );

                 }
                 double totalFlexRange = 2*numPoints*flexScanRes;



                   //loop over each NP component present
                    for(k = 0; k<numNPBeadSummation; ++k) {
              //distance is AA centre to NP (Nominal) surface
                            if(config.m_sumNPPotentials == false){
              distance = getDistance3D(x[j],y[j],z[j] + rcc,  np.m_x[k] , np.m_y[k] , np.m_z[k] , np.m_radius[np.m_npBeadType[k]]    ,np.m_shape[np.m_npBeadType[k]]);
              }
              else{
              distance = getDistance3D(x[j],y[j],z[j] + rcc,  0 , 0 , 0 , radius,npType);
              
              
                             if(distance < closestSCD[j]){
                //std::cout <<"applying overlap penalty "<<overlapPenalty << " at scd " << distance << " due to restriciton " << closestSCD[j] << "\n";
                appliedOverlapPenalty = overlapPenalty;
                }
              
              
              }

               if(distance < 0.5){ //count close-range interactions for later summary statistics
               resHasContacted = 1;
               }

                double noiseEnergy =   pdb.m_occupancy[j]*appliedOverlapPenalty ;   
                 
                
                 //std::cout << " sampling " << numPoints << "for residue with RMSD " << rmsd << "\n" ;

                if( flexMethod == 0){
                 //std::cout << "using basic fixed model \n";
                  noiseEnergy +=   pdb.m_occupancy[j] * static_cast<double>(potentials[ pdb.m_id[j]].Value(distance, np.m_npBeadType[k] ));
   
                }
                else if( flexMethod == 1){

                //double rmsd = 0.05; //in nanometers, typical bfactor of 25 = 0.5 A = 0.05 nm , rmsd^2 = bfactor/8 pi^2
                //double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                double flexEnergy = 0;
                //double flexEnergyDenom = 0.00000001;
                for( int di = -numPoints; di < numPoints+1; di++){
                double dOff = di*flexScanRes; 
                 double eWeight = exp(-dOff*dOff/(sigmaSq*2) )/ sqrt( 2 * M_PI * sigmaSq) ;
                 flexEnergy +=  eWeight * static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] ));
                 flexEnergyDenom += eWeight; 
                }
                noiseEnergy += pdb.m_occupancy[j] * flexEnergy/flexEnergyDenom;

                }
                else if( flexMethod == 2){

                //double rmsd = 0.05; //in nanometers, typical bfactor of 25 = 0.5 A = 0.05 nm , rmsd^2 = bfactor/8 pi^2
                //double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                double flexEnergy = 500;
                //double flexEnergyDenom = 0.00000001;
                for( int di = -numPoints; di < numPoints+1; di++){
                double dOff = di*flexScanRes;
                double flexPenalty = dOff*dOff/( 2 * sigmaSq);
                 //double eWeight = exp(-dOff*dOff/(sigmaSq*2) )/ sqrt( 2 * M_PI * sigmaSq) ;
                 double trialEnergy = static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) + flexPenalty;
                 flexEnergy = min(flexEnergy, trialEnergy) ;
                }
                noiseEnergy += pdb.m_occupancy[j] * flexEnergy;

                }
                /*
                else if( flexMethod == 3){
                double midpoint = 0;
                //double rmsd = 0.05; //in nanometers, typical bfactor of 25 = 0.5 A = 0.05 nm , rmsd^2 = bfactor/8 pi^2
                double dr = flexScanRes;
                //double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                double flexEnergy = 0;
                //flexEnergyDenom = 0
                //double flexEnergyDenom = 0.00000001;
                for( int di = -numPoints; di < numPoints+1; di++){
                double dOff = di*dr;
                //if(di == 0){
                //   midpoint = static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] ));
                //}
                double flexPenalty = dOff*dOff/( 2 * sigmaSq);
                 //double eWeight = exp(-dOff*dOff/(sigmaSq*2) )/ sqrt( 2 * M_PI * sigmaSq) ;

                //std::cout << "distance" << distance << "offset: " << dOff << " component: " << static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) << "penalty: " << flexPenalty << "\n" ; 
                 double trialEnergy = exp(-1.0 * pdb.m_occupancy[j] * ( static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) + flexPenalty+appliedOverlapPenalty));
                 //flexEnergy = min(flexEnergy, trialEnergy) ;
                   flexEnergy += trialEnergy*dr;
                   flexEnergyDenom += exp(-1.0 * pdb.m_occupancy[j] *flexPenalty)* dr;  
                   //std::cout << flexEnergyDenom << "\n"; 
                }
                //std::cout << flexEnergy << "/" << flexEnergyDenom << "\n";
                flexEnergyTot += flexEnergy; 
                flexEnergy = -log( flexEnergy/flexEnergyDenom) ;
                //std::cout <<"final: " << distance << ":" <<    flexEnergy << "\n";
                noiseEnergy +=  flexEnergy;

                }
                */

                else if( flexMethod == 3){
                double midpoint = 0;
                //double rmsd = 0.05; //in nanometers, typical bfactor of 25 = 0.5 A = 0.05 nm , rmsd^2 = bfactor/8 pi^2
                double dr = rmsd;
                //double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                double flexEnergy = 0;
                double flexEnergyDenom = 0.00000001;
                for( int di = -numPoints; di < numPoints+1; di++){
                double dOff = di*dr;
                //if(di == 0){
                //   midpoint = static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] ));
                //}
                double flexPenalty = dOff*dOff/( 2 * sigmaSq);
                 //double eWeight = exp(-dOff*dOff/(sigmaSq*2) )/ sqrt( 2 * M_PI * sigmaSq) ;

                //std::cout << "distance" << distance << "offset: " << dOff << " component: " << static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) << "penalty: " << flexPenalty << "\n" ;
                 double trialEnergy = exp(-1.0 *( static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) + flexPenalty));
                 //flexEnergy = min(flexEnergy, trialEnergy) ;
                   flexEnergy += trialEnergy*dr;
                   flexEnergyDenom += exp(-1.0 * flexPenalty)* dr;
                   //std::cout << flexEnergyDenom << "\n";
                }
                //std::cout << flexEnergy << "/" << flexEnergyDenom << "\n";
                flexEnergy = -log( flexEnergy/flexEnergyDenom) ;
                //std::cout <<"final: " << distance << ":" <<    flexEnergy << "\n";
                noiseEnergy += pdb.m_occupancy[j] * flexEnergy;

                }



                else if( flexMethod == 4){
                double midpoint = 0;
                //double rmsd = 0.05; //in nanometers, typical bfactor of 25 = 0.5 A = 0.05 nm , rmsd^2 = bfactor/8 pi^2
                double dr = flexScanRes;
                //double sigmaSq = rmsd*rmsd;   //approximate a Gaussian distribution
                double flexEnergy = 0;
                //flexEnergyDenom = 0
                //double flexEnergyDenom = 0.00000001;
                for( int di = -numPoints; di < numPoints+1; di++){
                double dOff = di*dr;
                double flexPenalty = dOff*dOff/( 2 * sigmaSq);
                 double trialEnergy =  pdb.m_occupancy[j] * ( static_cast<double>(potentials[ pdb.m_id[j]].Value(distance+dOff, np.m_npBeadType[k] )) + flexPenalty+appliedOverlapPenalty);
                 //std::cout << " central " << distance << " offset: " << dOff << " index: " << di+numPoints << " contribtion " << trialEnergy << "\n";
                 flexEnergyTerms[di+numPoints] += trialEnergy;  
                 flexEnergyDenomTerms[di+numPoints] += pdb.m_occupancy[j] * ( flexPenalty );
                }
                noiseEnergy += 0; // pdb.m_occupancy[j] * flexEnergy;

                }



                else{
                   //fall back to default method
                  noiseEnergy +=   pdb.m_occupancy[j] * static_cast<double>(potentials[ pdb.m_id[j]].Value(distance, np.m_npBeadType[k] ));

                }


   /*
                if(config.m_potNoiseMag > 100){
                 //noiseEnergy +=     unitNormalDist(randomEngine)*  (  std::min( std::fabs(noiseEnergy*0.1) , config.m_potNoiseMag)); //gaussian noise with zero mean and magnitude of 10% of the energy or the value specified in the config file, whichever is smaller

                 noiseEnergy +=  pdb.m_occupancy[j] *unitNormalDist(randomEngine) * config.m_potNoiseMag;
                }
                */
                
                energy += noiseEnergy;
                
                 if(std::isnan(pdb.m_occupancy[j] * static_cast<double>(potentials[pdb.m_id[j]].Value(distance, np.m_npBeadType[k] )))){
                  std::cout << "NaN generated by res. " << j << " type: "  << pdb.m_id[j]  <<   " at " << distance << "with NP bead " << k << " of type " << np.m_npBeadType[k] << "\n";
                    }

                }
                 

               //loop over k is ended, still have access to the individual residue here
               if(flexMethod==4){
                 double flexEnergyNumerator = 0;
                 flexEnergyDenom = 1e-10;
                  for( int di = -numPoints; di < numPoints+1; di++){
                     flexEnergyNumerator += exp( -flexEnergyTerms[di+numPoints] );
                     flexEnergyDenom += exp( -flexEnergyDenomTerms[di+numPoints] ) ;
                  }

               //flexEnergyDenom = totalFlexRange;
                 //if(distance<2){
                //std::cout << " residue " << j << " at " << distance <<" terms: " << flexEnergyNumerator << "/" << flexEnergyDenom << " final: " << -log( flexEnergyNumerator/flexEnergyDenom) << "\n";
                //}
                energy += -log( flexEnergyNumerator/flexEnergyDenom) ;
               }

                
                numContactsAtStep += resHasContacted;
                

                }//loop over the residue index j ends here
                


                 if(config.m_potNoiseMag > 0.01){
                 //std::cout << "adding potential noise \n";
                 energy = energy + unitNormalDist(randomEngine) *config.m_potNoiseMag;
                 //std::cout << "added potential noise \n";
                 }
                



                SSD[i]          = rcc;
                total_energy[i] = energy;
    
               //std::cout << theta*180/M_PI<< " " << phi*180/M_PI << " " << SSD[i] << " " << total_energy[i] << "\n";

               if(energy < minEnergy){
               minEnergy = energy;
               minLoc = rcc;
               numContactsAtMin = numContactsAtStep;
               }

            }

           if(sample == samples){
          //for the special run at the bin-center, write out the potential and exit the loop to avoid overloading the results arrays
           //std::cout << "preparing potential output \n";
if(savePotentials==1){
std::string filename;
std::string cylinderFileNameAppend;
if(npType == 2 || npType == 4 || npType==5){
cylinderFileNameAppend = "_"  + std::to_string(static_cast<int>(cylinderAngleDeg));
}
else{
cylinderFileNameAppend = "";
}

  //std::cout<<phi << " " << theta << "\n";
          filename = outputdirectory + "/" + npName + "/"+ pdbname + "_"+std::to_string(static_cast<int>(radius)) + "_" + std::to_string(static_cast<int>(1000 * zeta))+ "_" + std::to_string(static_cast<int>(phi*180/M_PI))+ "_" + std::to_string(static_cast<int>(theta*180/M_PI))+ cylinderFileNameAppend  +".uap";
    std::ofstream handle(filename.c_str());
             handle<<"#rcc,E(kbT)\n";
            for (i = 0; i < steps; ++i) {
                handle << (SSD[i] ) << ", " << total_energy[i] << "\n";
            }
    handle.close();
}
          continue;
          }



           sample_minloc[sample] = minLoc - radius; //this gives the distance between the protein's COM and the surface of the NP
           sample_numcontacts[sample] = numContactsAtMin;
           sample_proteinoffset[sample] =  appliedOffset;

         /*   for (i = 0; i < steps; ++i) {
                std::cout << (SSD[i] - radius) << " " << total_energy[i] << "\n"; 
            }*/
           // exit(0);

            // Integrate the results
                if(npType == 2 || npType == 4 || npType==5){
            IntegrateCylinder(steps, actualDZ, init_energy, total_energy, SSD, &(sample_energy[sample]),temperature);
            IntegrateMFPT(steps, actualDZ, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType,calculateMFPT,temperature);
            }
             else if(npType == 3){
             IntegrateCube(steps, actualDZ, init_energy, total_energy, SSD, &(sample_energy[sample]),temperature);
            IntegrateMFPT(steps, actualDZ, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType,calculateMFPT,temperature);
          }
           else{


            Integrate(steps, actualDZ, init_energy, total_energy, SSD, &(sample_energy[sample]), temperature);
           //std::cout << "integration done" << sample <<   "\n" ;
         //and also integrate to get the product of the MFPT and the diffusion constant, MFPT*D. 


            IntegrateMFPT(steps, actualDZ, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType, calculateMFPT,temperature);

 


        }



        }
  
        // Mean of all the samples
        
 
        if(config.m_enableLocalBoltz == true){
        BoltzMeanAndSD(samples, &(adsorption_energy[angle_offset + angle]), &(adsorption_error[angle_offset + angle]), sample_energy); 
        }
        else{
        MeanAndSD(samples, &(adsorption_energy[angle_offset + angle]), &(adsorption_error[angle_offset + angle]), sample_energy); 
        }
        
        if(angle_offset == 0 && (angle + 1) % threadPrintFreq == 0){
        std::cout << "*"  << std::flush ;
        }
        MeanAndSD(samples, &(mfpt_val[angle_offset + angle]), &(mfpt_err[angle_offset + angle]), sample_mfpt); 
         MeanAndSD(samples, &(minloc_val[angle_offset + angle]), &(minloc_err[angle_offset + angle]), sample_minloc);
        MeanAndSD(samples, &(numcontacts_val[angle_offset + angle]), &(numcontacts_err[angle_offset + angle]), sample_numcontacts);
         MeanAndSD(samples, &(protein_offset_val[angle_offset + angle]), &(protein_offset_err[angle_offset + angle]), sample_proteinoffset);
 
    }
}

void SurfaceScan(const PDB& pdb, const Potentials& potentials, const double zeta, const double radius, const double outerRadius,const Config& config, double omegaAngle, const NP& np) {
    std::clog << "\nInfo: Processing '" << pdb.m_name << "' (R = " << radius << ") NBeads: " << pdb.m_id.size() << "\n";

    const double imaginary_radius = config.m_imaginary_radius;

    double adsorption_energy[iterations] = {};
    double mfpt_val[iterations] = {};
    double adsorption_error[iterations]  = {};
        double mfpt_err[iterations] = {};
   double minloc_val[iterations] = {};
   double minloc_err[iterations] = {};
      double protein_offset_val[iterations] = {};
   double protein_offset_err[iterations] = {};
    double numcontacts_val[iterations] = {};
   double numcontacts_err[iterations] = {};  
   
int appendAngle = 0;   
 if(config.m_npType == 2 || config.m_npType == 4 || config.m_npType == 5){
 appendAngle = 1;
   }  
   if(config.m_npType == 1 && omegaAngle > 0.1){
   appendAngle = 1;
   }
   
    bool checkPrecalc = CheckForPrecalculated(adsorption_energy, adsorption_error, mfpt_val,minloc_val,radius, zeta, pdb.m_name, config.m_outputDirectory, np.m_name, 0,appendAngle,omegaAngle);
     //std::clog << checkPrecalc << "\n";
    if(checkPrecalc == true){
        std::clog << "Info: Target '" <<np.m_name << ":" << pdb.m_name << "' (R = " << radius << ") already calculated, skipping. \n";
    }
    else{
    #ifdef PARALLEL  
    const int n_threads         =    omp_get_max_threads();
    #else
    const int n_threads         = 1;
    #endif
  
    const int n_per_thread      = iterations / n_threads;
    const int n_remaining       = iterations % n_threads; 

    int nMaxSymbols = n_per_thread;
    //int nMaxSymbolsPerThread = ceil( nMaxSymbols/n_threads) ;

    //int iterationsPerSymbol = ceil( n_threads * n_per_thread/nMaxSymbols);
    threadPrintFreq =   1; // iterationsPerSymbol ;


    for(int p =0; p< nMaxSymbols; ++p){
    std::cout << "-" ;
    }
    std::cout << "\n" << std::flush;
    //std::cout << " total iterations: " << iterations << " threads selected: " << n_threads << "\n";
    //std::cout << "each thread handles " << n_per_thread << " such that one symbol should be printed every " << iterationsPerSymbol << "\n";
    int thread;
    #ifdef PARALLEL  
    #pragma omp parallel shared(pdb, potentials, adsorption_energy, adsorption_error), private(thread)
    #endif
    {
        #ifdef PARALLEL  
        #pragma omp for schedule(dynamic, 1)
        #endif
        for (thread = 0; thread < n_threads; ++thread) {
            AdsorptionEnergies(
                    pdb, 
                    np,
                    config,
                    potentials, 
                    radius, 
                    outerRadius,
                    thread * n_per_thread  + (thread < n_remaining ? thread : n_remaining), 
                    n_per_thread + (thread < n_remaining), 
                    adsorption_energy, 
                    adsorption_error,
                    mfpt_val,
                    mfpt_err,
                    minloc_val,
                    minloc_err,
                    numcontacts_val,
                    numcontacts_err,
                    protein_offset_val,
                    protein_offset_err,
                    pdb.m_name,
                    config.m_outputDirectory,
                    np.m_name,
                    config.m_npType,
                    imaginary_radius,
                    config.m_calculateMFPT,
                    omegaAngle,
                    zeta,
                    config.m_savePotentials,
                    config.m_temperature,
                    config.m_overlapPenalty

            ); 
        }
    }
    


    //find minimum

    double thetaMin = 0;
    double phiMin = 0;
    double minEnergy = 50;
    double ccdAtMin = 0;
    for (int i = 0; i < iterations; ++i) {
        if(adsorption_energy[i] < minEnergy){
        minEnergy = adsorption_energy[i];
       ccdAtMin = minloc_val[i] + radius;
        thetaMin = (i / ncols) * angle_delta;
        phiMin = (i % ncols) * angle_delta;
        }
    }
    std::cout << "\nMinimum located at phi=" << phiMin *180.0/M_PI << ", theta="<< thetaMin*180.0/M_PI << " E="<<minEnergy <<"\n";
    bool bSaveMinimum = config.m_saveOptimumPDB;
    if(bSaveMinimum==true){
    SaveSpecificOrientation(phiMin, thetaMin, pdb.m_id.size() ,  config, potentials , pdb,  np, radius, config.m_npType, omegaAngle* M_PI/180.0, pdb.m_name, config.m_outputDirectory, np.m_name, ccdAtMin) ;
    }
    WriteMapFile(adsorption_energy, adsorption_error, mfpt_val,minloc_val, numcontacts_val, protein_offset_val, radius, zeta, pdb.m_name, config.m_outputDirectory, np.m_name, config.m_temperature,     0,appendAngle,omegaAngle); 
    //WriteMapFile(mfpt_val, mfpt_err, radius, zeta, pdb.m_name, config.m_outputDirectory,1,isCylinder,cylinderAngle); 
    PrintStatistics(adsorption_energy, adsorption_error, radius, pdb.m_name);
    
    
    }
}

int main(const int argc, const char* argv[]) {
    CommandLineParser commandLine(argc, argv);
    ConfigFileReader  configFile(commandLine.ConfigFileName());

    Config config;
    configFile.UpdateConfig(config);
    commandLine.UpdateConfig(config);

    HamakerConstants  hamakerConstants(config.m_hamakerFile);
    TargetList        targetList(config.m_pdbTargets);
    SurfacePMFs       surfaces(config.m_pmfDirectory, config.m_pmfPrefix, config.m_aminoAcids);


    //if(config.m_readLigands == true){
    LigandMap         ligandMap(config.m_ligandFile);
    //}


    PDBs              pdbs(targetList.m_paths, config.AminoAcidIdMap()  ,   ligandMap.m_ligandAALookup, config.m_backboneReplaceSet, config.m_disorderStrat , config.m_disorderMinBound, config.m_disorderMaxBound, config.m_readLigands, config.m_bondCutoffNM, config.m_enableBackbone, config.m_backboneTag, config.m_backboneScale);



       samples = config.m_numSamples;



        //angle_delta_deg = 2.0;
        //angle_delta     = angle_delta_deg * (M_PI / 180.0);
       if(config.m_confirmOverrideAngle == true){
        //attempt override
         int intAngleDeg = int(  config.m_angleDelta ) ;
         if( 360 % intAngleDeg == 0 && 180 % intAngleDeg == 0 ){
        
        angle_delta_deg = config.m_angleDelta; 
        angle_delta     = angle_delta_deg * (M_PI / 180.0);
        ncols           = 360/intAngleDeg;
        nrows           = 180/intAngleDeg;



       
       std::cout << "Overriding angular resolution. New ncols: " << ncols << " new nrows: "<< nrows << "\n";
        iterations      = nrows * ncols;
        }
        else{
       std::cout << "Warning: you have chosen an override value of angle-delta with non-zero modulus w.r.t theta or phi \n";
       std::cout << "Suggested value is 5. \n";
        }
        }

    std::cout << "Total number of calculations per PDB per NP: " << iterations*samples << "\n"; 

    std::vector<double>          omegaAngleSet;
    // config.m_npTargets is a list of strings - targets to look at for .np files 
    
    
    int omegaDelta = 180;
    if(config.m_npType == 2 || config.m_npType == 4 || config.m_npType == 5){
    omegaDelta = 45;
     }
     //if final rotation angles omega have been manually set in the config file, use these. Else generate them based on the NP type.
     
    if( config.m_omegaAngles.size() == 0){
    for(int i = 0; i < 180; i = i + omegaDelta){
    omegaAngleSet.push_back(i);
    std::cout << "Adding automatic rotation angle: " << i << "\n";
     }
    }
    else{
    for( auto & omegaAngle : config.m_omegaAngles){
    omegaAngleSet.push_back(  omegaAngle ) ;
        std::cout << "Adding rotation angle: " << omegaAngle << "\n";
    }
    }
     
    boost::filesystem::create_directory(config.m_outputDirectory);
    boost::filesystem::create_directory(config.m_outputDirectory+"/nps");
    boost::filesystem::create_directory(config.m_outputDirectory+"/opt_pdbs");

    
    std::string configFileIn = commandLine.ConfigFileName();
    //Uncomment this line to re-enable config file saving
    //boost::filesystem::copy_file(configFileIn, config.m_outputDirectory+"/"+configFileIn, boost::filesystem::copy_options::overwrite_existing);
    
    //Or uncomment this line if you're on Ubuntu 18.04 and Boost is locked to a version with the slightly different naming convention that they then changed
    // boost::filesystem::copy_file(configFileIn, config.m_outputDirectory+"/"+configFileIn, boost::filesystem::copy_option::overwrite_if_exists);


     //generate some NPs 
        std::string npFileDir = config.m_outputDirectory+"/nps/"+config.m_configFile +  "_NPs";
              boost::filesystem::create_directory(npFileDir);
     //manually generate NPs , save them to a folder for storage, read this folder in as input to make sure generated NPs are handled identically to read-in NPs and so that their structure is kept
    if(   config.m_multiNP == 0){   
         int npID = 1;
    for (const double nanoparticleRadius : config.m_nanoparticleRadii) {
        for (const double zetaPotential : config.m_zetaPotential) {
            
             //std::to_string(static_cast<int>(radius))    + "_" + std::to_string(static_cast<int>(1000 * zeta))
            //x,y,z,radius,zeta,coreFactor=1,surfFactor=1,shape, hamakerFile,pmfFile,pmfCutoff,correctionType
            std::string npTypeOutString = "TYPE," + to_string(nanoparticleRadius)+"," + to_string(zetaPotential) + "," + "1,1," + to_string(config.m_npType) +","+ config.m_hamakerFile +","+ config.m_pmfDirectory+"," + to_string(config.m_PMFCutoff) + ","+to_string( config.m_npType )+"\n";
            std::cout << npTypeOutString;
            
           std::string npBeadOutString = "BEAD,0,0,0,0\n";
            std::cout << npBeadOutString;
            
            std::string npIDString = "np"+to_string(npID)+"R_"+to_string(static_cast<int>(nanoparticleRadius))+"_ZP_"+to_string(  static_cast<int>( zetaPotential*1000));
            std::string npOutLoc = npFileDir+"/"+npIDString+".np";
            npID++;
            
             std::ofstream handle(npOutLoc.c_str());
             handle<<"#NP file generated by UnitedAtom\n";
             handle << npTypeOutString;
             handle << npBeadOutString;
             handle.close();
            
        }
    }
    config.m_npTargets.emplace_back(npFileDir); 
    }
    // 
        NPTargetList npTargetList(config.m_npTargets);
        NPs nps(npTargetList.m_paths, config.AminoAcidIdMap());
    
   //scan over multiple NP files 

    for( const auto& np: nps){
    std::string npEnergyDir = config.m_outputDirectory+"/"+np.m_name;
    boost::filesystem::create_directory(npEnergyDir);
    double zetaPotential = np.m_zetaName;
    double nanoparticleBoundingRadius;
    double nanoparticleOuterBoundingRadius;
    if(config.m_boundingRadius < 0){
    nanoparticleBoundingRadius = np.m_boundRadius;
    nanoparticleOuterBoundingRadius = np.m_outerBoundRadius;
    }
    else{
    nanoparticleBoundingRadius = config.m_boundingRadius;
    nanoparticleOuterBoundingRadius = config.m_boundingRadius+0.01;
    }
    std::cout << "Scanning from (approx):" << nanoparticleBoundingRadius << " to " << nanoparticleOuterBoundingRadius + 2.0 << "\n";
    Potentials potentials(surfaces, hamakerConstants, zetaPotential, nanoparticleBoundingRadius, config,np);
    for (const auto& pdb : pdbs){
    

    //SpecificPDBPotentials 
 

    for( auto & omegaAngle : omegaAngleSet){
    SurfaceScan(pdb,potentials,zetaPotential,nanoparticleBoundingRadius,nanoparticleOuterBoundingRadius,config,omegaAngle,np);
    }
    
    
    }
    }
    

    return 0;
}

