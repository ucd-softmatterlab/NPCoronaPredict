#ifndef POTENTIALS__H__
#define POTENTIALS__H__

#include "Surface.h"
#include "HamakerConstants.h"
#include "Config.h"

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
double HamakerPotential(const double, const double, const double, const double, const double);
double ElectrostaticPotential(const double, const double, const double, const double, const double, const double);

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
    const double        bjerumLength        = config.m_bejerumLength;
    const double        debyeLength         = config.m_debyeLength;
    const auto          aminoAcidRadius     = config.AminoAcidRadius(surfaceData.m_aminoAcid);
    const double        Z1                  = config.AminoAcidCharge(surfaceData.m_aminoAcid);
    const double        Z2                  = 38.681 * zetaPotential * nanoparticleRadius * (1.0 + 4.0 * M_PI * bjerumLength) / bjerumLength;
    
    std::vector<double> energy(potential_size);

    SurfacePotential surface_potential(surfaceData.m_energy, surfaceData.m_distance);
    
    double  surface = 0.0, core = 0.0, electrostatic = 0.0;

    const std::string destination = "/home/dpower/Git/unitedatom/pot-dat";
    const std::string filename = destination + "/" + surfaceData.m_aminoAcid + ".dat";
    std::ofstream handle(filename.c_str());

    handle << "# Distance    Total         Surface       Core          Electrostatic\n";

    for (int i = 0; i < potential_size; ++i) {
        const double r = pmfStart + (i / (potential_size - 1.0)) * (cutoff - pmfStart);
        double U       = 0.0;

        if (config.m_enableSurface) {
            surface = surface_potential.Value(r, nanoparticleRadius, pmfCutoff);
            U += surface;
        }

        if (config.m_enableCore) {
            core = HamakerPotential(hamaker, aminoAcidRadius, nanoparticleRadius, r, pmfCutoff);
            U += core;
        }

        if (config.m_enableElectrostatic) {
            electrostatic = ElectrostaticPotential(r, Z1, Z2, nanoparticleRadius, bjerumLength, debyeLength);
            U += electrostatic;
        }

        energy[i] = (double) U;
    
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

double ElectrostaticPotential(const double r, const double Z1, const double Z2, const double nanoparticleRadius, const double bjerumLength, const double debyeLength) {
    return (Z1 * Z2 * bjerumLength * std::exp((-1.0 * r) / debyeLength)) / (r + nanoparticleRadius);
}

#endif
