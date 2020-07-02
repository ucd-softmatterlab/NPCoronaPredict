#include <cmath>
#include <iostream>

 double HamakerPotential(const double A, const double R1, const double R2, const double c) {
     return -1.0 * (A / 6.0) * (
       (2.0 * R1 * R2) / (c * c - (R1 + R2) * (R1 + R2)) +
       (2.0 * R1 * R2) / (c * c - (R1 - R2) * (R1 - R2)) +
       std::log((c * c - (R1 + R2) * (R1 + R2)) / (c * c - (R1 - R2) * (R1 - R2)))
     );
 }

 double HamakerPotential(const double A, const double R1, const double R2, const double c, const double cutoff) {
 
   if ((c - R2) > cutoff) { // Classical Hamaker Potential
     return -1.0 * (A / 6.0) * (
       (2.0 * R1 * R2) / (c * c - (R1 + R2) * (R1 + R2)) +
       (2.0 * R1 * R2) / (c * c - (R1 - R2) * (R1 - R2)) +
       std::log((c * c - (R1 + R2) * (R1 + R2)) / (c * c - (R1 - R2) * (R1 - R2)))
     );
   }
   else { // Modified Hamaker Potential (with missing lens piece)
     double top = 2.0 * R1 * (c + R2 - cutoff) * (-1.0 * c * c * c * R1 * R1 + R1 * R1 * (R2 - cutoff) * std::pow(R2 + cutoff, 2.0)
       + c * c * (-1.0 * R1 * R1 * R2 + std::pow(cutoff, 3.0)) + c * (R1 * R1 * R1 * R1 + R2 * std::pow(cutoff, 3.0) + R1 * R1 * R2 * (R2 + cutoff)));
     double bottom = c * (c - R1 + R2) * (c + R1 + R2) * std::pow(R1 - cutoff, 2.0) * std::pow(R1 + cutoff, 2.0);
     double logterm = -1.0 * std::log(c - R1 + R2) + std::log(c + R1 + R2) + std::log(cutoff - R1) - 1.0 * std::log(R1 + cutoff);
     return -1.0 * (A / 6.0) * ((top / bottom) + logterm);
   }
 }


int main() {

    int     size    = 500;
    double  start   = 0.001;
    double  stop    = 1.5;
    double  A       = 1;
    double  R1      = 0.5;
    double  R2      = 20.0;
    double  cutoff  = 1.0;


    for (int i = 0; i < size; ++i) {

        double d = start + (stop - start) * (i / (size - 1.0));
        double c = d + R1 + R2; 
        
        double U0 = HamakerPotential(A, R1, R2, c);
        double U1 = HamakerPotential(A, R1, R2, c, cutoff);
            
        std::cout << d << " " << U0 << " " << U1 << "\n";
    }

    return 0;
}
