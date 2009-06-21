/*************************************************************/
/*** RETURNS THE DENSITY AND PRESSURE AND ENERGY IN CGS    ***/
/*                                                           */
/* Takes the chemical potential and temperature in units     */
/* of the rest mass.                                         */
/*                                                           */
/*************************************************************/


#ifndef FERMIDIRAC_H
#define FERMIDIRAC_H

#include <cmath>
#include "constants.h"
#include "fermidiracintegrator.h"

class FermiDirac {
 public:
  FermiDirac(double m, double s, double eps=1.0e-8, double infty=1e4);

  double number_density(double mu, double T); // number density
  double energy_density(double mu, double T); // energy per unit volume
  double pressure(double mu, double T);       // pressure

  inline double mass() const { return _m; }
  inline double spin() const { return (_g-1.0)/2.0; }

 private:
  double _m, _g;        // mass and degeneracy
  double _eps, _infty;  // tunable parameters for integrations

  FermiDiracIntegrator _fdi;

  class number_density_integrand : public Func
  {
    virtual double operator()(double E) {
      return 1.0;
    }
  };
  class energy_density_integrand : public Func
  {
    virtual double operator()(double E) {
      return E;
    }
  };
  class pressure_integrand : public Func
  {
    virtual double operator()(double E) {
      return (E*E-1.0)/E;
    }
  };
};

FermiDirac::FermiDirac(double m, double s, double eps, double infty)
  : _m(m), _g(2.0*s+1.0), _eps(eps), _infty(infty)
{
}

double FermiDirac::number_density(double mu, double T)
{
  number_density_integrand integrand;
  double val = _fdi.integrate(integrand,mu,T,_eps,_infty);
  val *= _g * 4.0*M_PI * std::pow(_m * CGS_Constants::c/CGS_Constants::h, 3.0); // Normalization to real units
  return val;
}

double FermiDirac::energy_density(double mu, double T)
{
  energy_density_integrand integrand;
  double val = _fdi.integrate(integrand,mu,T,_eps,_infty);
  val *= _g * 4.0*M_PI * std::pow(_m * CGS_Constants::c/CGS_Constants::h, 3.0) * _m * CGS_Constants::c*CGS_Constants::c;
  return val;
}

double FermiDirac::pressure(double mu, double T)
{
  pressure_integrand integrand;
  double val = _fdi.integrate(integrand,mu,T,_eps,_infty);
  val *= _g * (4.0*M_PI/3.0) * std::pow(_m * CGS_Constants::c/CGS_Constants::h, 3.0) * _m * CGS_Constants::c*CGS_Constants::c;
  return val;
}


#endif
