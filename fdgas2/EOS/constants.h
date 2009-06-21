/***************************************************************************
                          constants.h  -  description
                             -------------------
    begin                : Mon Apr 16 2001
    copyright            : (C) 2001 by Yasser Rathore
    email                : yasser@caltech.edu
 ***************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace CGS_Constants {   // All constants in SI units
  static const double G = 6.6726e-8;   // Gravitational constant
  static const double c = 2.99792458e10;   // Speed of light in vacuum
  static const double h = 6.6260755e-27;   // Planck constant
  static const double hbar = 1.05457267e-27;   // Reduced Planck constant
  static const double e = 4.8032e-10;   // Electron charge magnitude
  static const double me = 9.1093898e-28;   // Electron mass
  static const double mp = 1.6726231e-24;   // Proton mass
  static const double mu = 1.6605402e-24;   // Unified atomic mass unit
  static const double k = 1.380658e-16;   // Boltzmann constant
}

namespace MKS_Constants {   // All constants in SI units
  static const double G = 6.6726e-11;   // Gravitational constant
  static const double c = 2.99792458e8;   // Speed of light in vacuum
  static const double h = 6.6260755e-34;   // Planck constant
  static const double hbar = 1.05457267e-34;   // Reduced Planck constant
  static const double e = 1.60217733e-19;   // Electron charge magnitude
  static const double me = 9.1093898e-31;   // Electron mass
  static const double mp = 1.6726231e-27;   // Proton mass
  static const double mu = 1.6605402e-27;   // Unified atomic mass unit
  static const double epsilon0 = 8.854187817e-12;   // Permittivity of free space
  static const double mu0 = 1.2566370614e-6;   // Permeability of free space
  static const double k = 1.380658e-23;   // Boltzmann constant
  static const double a = 7.56591e-16;   // Stefan-Boltzmann radiation constant
  static const double SolarMass = 1.98892e30;   // kg
  static const double SolarRadius = 6.96e8;   // m
  static const double SolarDensity = (3.0*SolarMass/(4.0*M_PI*std::pow(SolarRadius,3)));
}

#endif
