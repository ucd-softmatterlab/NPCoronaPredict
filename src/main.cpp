#define PARALLEL

#include "Config.h"
#include "CommandLineParser.h"
#include "ConfigFileReader.h"
#include "HamakerConstants.h"
#include "TargetList.h"
#include "Surface.h"
#include "PDBFile.h"
#include "Potentials.h"

#include <omp.h>
#include <cmath>
#include <random>
#include <iomanip>

//constexpr double        gds             = 0.22;
constexpr double        gds             = 0.0;
constexpr double        delta           = 2.0;
constexpr double        angle_delta     = 5.0 * (M_PI / 180.0);
constexpr int           ncols           = 72;
constexpr int           nrows           = 36;
constexpr int           iterations      = nrows * ncols;
constexpr int           samples         = 128;
constexpr int           steps           = 512;
constexpr double        dz              = delta / (steps - 1);

std::random_device randomEngine;
std::uniform_real_distribution<double> random_angle_offset(0.0, angle_delta);

void WriteMapFile(const double *adsorption_energy, const double *adsorption_error, const double *mfpt_val, const double *minloc_val, const double radius,
        const double zeta, const std::string& name, const std::string& directory, int isMFPT=0, int isCylinder = 0, double cylinderAngle=0) {
std::string filename;
std::string cylinderFileNameAppend;
if(isCylinder == 0){
cylinderFileNameAppend = "";
} 
else{
cylinderFileNameAppend = "_"  + std::to_string(static_cast<int>(cylinderAngle));
}

    if(isMFPT==1){
     filename = directory + "/" + TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  "_mfpt.uam";
    }
else{
    filename = directory + "/" + TargetList::Filename(name) + "_" + std::to_string(static_cast<int>(radius))
        + "_" + std::to_string(static_cast<int>(1000 * zeta)) + cylinderFileNameAppend  +  ".uam";
}
    std::clog << "Info: Saving map to: "<< filename << "\n";
    std::ofstream handle(filename.c_str());
    double phi, theta;
    handle << "#phi theta EAds/kbT Error(Eads)/kbT min_surf-surf-dist/nm mfpt*DiffusionCoeff/nm^2 \n"; 
    for (int i = 0; i < iterations; ++i) { 
        phi   = (i % ncols) * 5.0;
        theta = (i / ncols) * 5.0;
        handle << std::left << std::setw(7) << std::fixed << std::setprecision(1) << phi;
        handle << std::left << std::setw(7) << std::fixed << std::setprecision(1) << theta;

        
        //handle << std::left << std::setw(14) << std::fixed << std::scientific << std::setprecision(5) << adsorption_energy[i];
        
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << adsorption_energy[i];
        
        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << adsorption_error[i];

        handle << std::left << std::setw(14) << std::fixed << std::setprecision(5) << minloc_val[i];

        handle << std::left << std::setw(14) << std::fixed << std::scientific << std::setprecision(5) << mfpt_val[i];

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
    for (int i = 0; i < iterations; ++i) {
        sin_theta       = std::sin((i / ncols) * angle_delta);
        T               += adsorption_energy[i] * sin_theta * std::exp(-1.0 * adsorption_energy[i]);
        Z               += sin_theta * std::exp(-1.0 * adsorption_energy[i]);
        simpleAverage   += sin_theta * adsorption_energy[i];
        error           += sin_theta * adsorption_error[i];
        sinThetaTot     += sin_theta;
    }

    simpleAverage   /= sinThetaTot;
    error           /= iterations;

    std::cout << std::setw(10) << std::left << TargetList::Filename(filename);
    std::cout << std::setw(10) << std::fixed << std::setprecision(1) << radius;
    std::cout << std::setw(14) << std::fixed << std::setprecision(5) << simpleAverage;
    std::cout << std::setw(14) << std::fixed << std::setprecision(5) << (T/Z);
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
//note that for this integration routine and all others, "ssd" is actually the distance between the centre of the NP and the closest point of the protein.
void Integrate(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption) {
    long double area = 0.0;
    for (int i = 0; i < size; ++i) {
        area += static_cast<long double>(ssd[i] * ssd[i] * dz) * std::exp(static_cast<long double>(-1.0 * (energy[i] - init_energy))); 
    }
    const double factor = 3.0 / std::fabs(std::pow(ssd[0], 3.0) - std::pow(ssd[size - 1], 3.0));
    *adsorption = -1.0 * std::log(factor * area);
}


/*


*/

void IntegrateMFPT(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *mfpt, int npType, int calculateMFPT) {
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



void IntegrateCylinder(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption) {
    long double area = 0.0;
    for (int i = 0; i < size; ++i) {

if(std::isnan(energy[i])){
std::cout << "found nan at " << i << " " <<   "\n";
}
if(std::isinf(energy[i])){
std::cout << "found inf at " << i << " " <<   "\n";
}

        area += static_cast<long double>(ssd[i]  * dz) * std::exp(static_cast<long double>(-1.0 * (energy[i] - init_energy))); 
    }
    const double factor = 2.0 / std::fabs(std::pow(ssd[0], 2.0) - std::pow(ssd[size - 1], 2.0));

    if(factor*area < 0){ 
 std::cout << "warning: factor*area < 0, unphysical result " << factor << " " << area << "\n";    
}



    *adsorption = -1.0 * std::log(factor * area);
}


void IntegrateCube(const int size, const double dz, const double init_energy, const double *energy, const double *ssd, double *adsorption) {
    long double area = 0.0;
    for (int i = 0; i < size; ++i) {
        area += static_cast<long double>(dz  ) * std::exp(static_cast<long double>(-1.0 * (energy[i] - init_energy))); 
    }
    const double factor = 1.0 / std::fabs(std::pow(ssd[0], 1.0) - std::pow(ssd[size - 1], 1.0));
    *adsorption = -1.0 * std::log(factor * area);
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

void AdsorptionEnergies(const PDB& pdb, const Potentials& potentials, const double radius, const int angle_offset, const int n_angles, double *adsorption_energy, double *adsorption_error, double *mfpt_val, double *mfpt_err, double *minloc_val, double *minloc_err, int npType = 1, double imaginary_radius = -1, int calculateMFPT = 0, double cylinderAngleDeg=0) { 

    // Decleare all variables at the begining 
    const int               size            = pdb.m_id.size();
    const double            stop            = imaginary_radius < 0 ? radius + gds : imaginary_radius + gds;
    const double            start           = stop + delta;

    int                     i;
    int                     j;
    double                  energy;
    double                  theta;
    double                  phi;
    double                  theta_adjusted;
    double                  phi_adjusted;
    double                  init_energy;
    double                  distance;
    double                  ssd;

    double x[size];
    double y[size];
    double z[size];
    
    double total_energy[steps];
    double SSD[steps];
    double sample_energy[samples];
    double sample_mfpt[samples];
    double sample_minloc[samples];
    // Iterate over angles
    for (int angle = 0; angle < n_angles; ++angle) {
        
        phi   = ((angle_offset + angle) % ncols) * angle_delta;
        theta = ((angle_offset + angle) / ncols) * angle_delta;
           double cylinderAngle = cylinderAngleDeg * M_PI/180.0;

	    // Sample a angle multiple times 
        for (int sample = 0; sample < samples; ++sample) {
  //         double cylinderAngle = 90 * M_PI/180.0;

    	    // Rotation step
    	    phi_adjusted   = -1.0 * (phi + random_angle_offset(randomEngine));
    	    theta_adjusted = M_PI - (theta + random_angle_offset(randomEngine));
            Rotate3(size, phi_adjusted, theta_adjusted,cylinderAngle, pdb.m_x, pdb.m_y, pdb.m_z, x, y, z);
               
//           double cylinderAngle = 90 * M_PI/180.0;
          // if(npType == 2 || npType == 4 || npType==5){
             //
         //  Rotate(size, cylinderAngle,0, xtemp, ytemp, ztemp, x, y, z);

           // }

            // Convert to SSD
            ShiftZ(size, z);
   //         std::cout << z[0] << " " << phi << " " << theta << "\n";               
            // Pre-square x and y
            SquareXY (size, x, y);
            
            // Get bulk energy at the STARTING point
            init_energy = 0.0;

            for (i = 0; i < size; ++i) {
                if(npType == 2 || npType == 4 || npType==5){
                //for cylinder NPs we only take into consideration the radial distance with the cylinder aligned along the x-axis.
                distance    = std::sqrt(y[i] +  (z[i] + start) * (z[i] + start)) - radius; // AA Center To NP Surface Distance (originally stop, changed to start)
                }
                else if(npType == 3){
                //for cubes we consider only the "vertical" distance, i.e. we assume that each bead is approximately at the centre of the cube
                distance = std::sqrt(  (z[i] +start) * (z[i] + start)) - radius;
                }



                else{
                distance    = std::sqrt(x[i] + y[i] + (z[i] + start) * (z[i] + start)) - radius; // Center To Surface Distance
//               std::cout << i << " "  << phi << " "  << theta << " " <<  x[i] << " " << y[i] << " " << z[i] << " " <<distance << "\n";

                }
                double energyAtDist = pdb.m_occupancy[i] *  static_cast<double>(potentials[pdb.m_id[i]].Value(distance))  ;
                if(energyAtDist > 1){
                //std::cout << "large shift (" << energyAtDist << ") for distance " << distance << "\n";
                }
                //std::cout << pdb.m_occupancy[i] << "\n";
                init_energy +=  pdb.m_occupancy[i] *  static_cast<double>(potentials[pdb.m_id[i]].Value(distance));
            }
 
//std::cout << "total shift/num. amino acid: " <<  init_energy/size << "\n";
            // build the distance-potential array
            double minLoc = start;
            double minEnergy = 0;
            for (i = 0; i < steps; ++i) {
                ssd     = start - i * dz;
                energy  = 0;
                //loop over each residue
                for (j = 0; j < size; ++j) {
                if(npType == 2 || npType == 4 || npType==5){
                //for cylinder NPs we only take into consideration the radial distance z as the NP is assumed to be sufficiently long that the edge effects can be neglected.
                distance    = std::sqrt(y[j] +  (z[j] + ssd) * (z[j] + ssd)) - radius; // Center To Surface Distance
                }
                else if(npType == 3){
                //for cubes we consider only the "vertical" distance, i.e. we assume that each bead is approximately at the centre of the cube
                distance = std::sqrt(  (z[j] + ssd) * (z[j] + ssd)) - radius;
                }
                else{
                    distance  = std::sqrt(x[j] + y[j] + (z[j] + ssd) * (z[j] + ssd)) - radius; // Center To Surface Distance
                }
                energy +=  pdb.m_occupancy[j] * static_cast<double>(potentials[pdb.m_id[j]].Value(distance));
                }
                
                SSD[i]          = ssd;
                total_energy[i] = energy;
               if(energy < minEnergy){
               minEnergy = energy;
               minLoc = ssd;
               }

            }
           sample_minloc[sample] = minLoc - radius;
         /*   for (i = 0; i < steps; ++i) {
                std::cout << (SSD[i] - radius) << " " << total_energy[i] << "\n"; 
            }*/
           // exit(0);

            // Integrate the results
                if(npType == 2 || npType == 4 || npType==5){
            IntegrateCylinder(steps, dz, init_energy, total_energy, SSD, &(sample_energy[sample]));
            IntegrateMFPT(steps, dz, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType,calculateMFPT);
}
else if(npType == 3){
     IntegrateCube(steps, dz, init_energy, total_energy, SSD, &(sample_energy[sample]));
            IntegrateMFPT(steps, dz, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType,calculateMFPT);
}
else{


            Integrate(steps, dz, init_energy, total_energy, SSD, &(sample_energy[sample]));

         //and also integrate to get the product of the MFPT and the diffusion constant, MFPT*D. 


            IntegrateMFPT(steps, dz, init_energy, total_energy, SSD, &(sample_mfpt[sample]), npType, calculateMFPT);

 


}
        }
  
        // Mean of all the samples
        MeanAndSD(samples, &(adsorption_energy[angle_offset + angle]), &(adsorption_error[angle_offset + angle]), sample_energy); 
        MeanAndSD(samples, &(mfpt_val[angle_offset + angle]), &(mfpt_err[angle_offset + angle]), sample_mfpt); 
         MeanAndSD(samples, &(minloc_val[angle_offset + angle]), &(minloc_err[angle_offset + angle]), sample_minloc);



 
    }
}

void SurfaceScan(const PDB& pdb, const Potentials& potentials, const double zeta, const double radius, const Config& config, double cylinderAngle) {
    std::clog << "Info: Processing '" << pdb.m_name << "' (R = " << radius << ")\n";

    const double imaginary_radius = config.m_imaginary_radius;

    double adsorption_energy[iterations] = {};
    double mfpt_val[iterations] = {};
    double adsorption_error[iterations]  = {};
        double mfpt_err[iterations] = {};
   double minloc_val[iterations] = {};
   double minloc_err[iterations] = {};


    #ifdef PARALLEL  
    const int n_threads         =   omp_get_max_threads();
    #else
    const int n_threads         = 1;
    #endif
  
    const int n_per_thread      = iterations / n_threads;
    const int n_remaining       = iterations % n_threads; 

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
                    potentials, 
                    radius, 
                    thread * n_per_thread  + (thread < n_remaining ? thread : n_remaining), 
                    n_per_thread + (thread < n_remaining), 
                    adsorption_energy, 
                    adsorption_error,
                    mfpt_val,
                    mfpt_err,
                    minloc_val,
                    minloc_err,
                    config.m_npType,
                    imaginary_radius,
                    config.m_calculateMFPT,
                    cylinderAngle
            ); 
        }
    }
    
int isCylinder = 0;   
 if(config.m_npType == 2 || config.m_npType == 4 || config.m_npType == 5){
 isCylinder = 1;
   }  
  
    WriteMapFile(adsorption_energy, adsorption_error, mfpt_val,minloc_val,radius, zeta, pdb.m_name, config.m_outputDirectory,0,isCylinder,cylinderAngle); 
    //WriteMapFile(mfpt_val, mfpt_err, radius, zeta, pdb.m_name, config.m_outputDirectory,1,isCylinder,cylinderAngle); 
    PrintStatistics(adsorption_energy, adsorption_error, radius, pdb.m_name);
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
    PDBs              pdbs(targetList.m_paths, config.AminoAcidIdMap());
    int omegaDelta = 180;
    if(config.m_npType == 2 || config.m_npType == 4 || config.m_npType == 5){
    omegaDelta = 45;
 }

    for (const double nanoparticleRadius : config.m_nanoparticleRadii) {
        for (const double zetaPotential : config.m_zetaPotential) {
            Potentials potentials(surfaces, hamakerConstants, zetaPotential, nanoparticleRadius, config);
            for (const auto& pdb : pdbs) {
                for(int i = 0; i < 180; i = i + omegaDelta){
                SurfaceScan(pdb, potentials, zetaPotential, nanoparticleRadius, config,i);
                }
            }
        }
    }

    return 0;
}

