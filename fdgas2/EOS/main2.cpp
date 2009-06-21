#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "constants.h"
#include "numericalrecipes.h"
#include "fermidirac.h"


class ChargeEquilFunc
{
 public:
  ChargeEquilFunc(double m, double s, double eps=1.0e-8, double infty=1.0e4) : _FD(m,s,eps,infty) {};

  void set_T(double T) { _T=T; };
  void set_n(double n) { _n=n; };

  double operator()(double mum1) {
    double mu = mum1*_T + 1;  // Using mu-1 in units of temperature

    double nm = _FD.number_density(mu,_T);
    double na = _FD.number_density(-mu,_T);

    /*
    std::cout << std::setw(15) << mum1
	      << std::setw(15) << (_n - nm + na)/(_n+nm-na)
	      << '\n';
    */

    return (_n - nm + na)/(_n+nm-na);
  }

 private:
  FermiDirac _FD;
  double _T; // Temperature in units of rest mass
  double _n; // Number density in cgs units
};


int main(int argc, char *argv[])
{
  FermiDirac fd(CGS_Constants::me,0.5,1.0e-9,1.0e4);
  /*
  double T = 100 * (CGS_Constants::k/(CGS_Constants::me*CGS_Constants::c*CGS_Constants::c));
  std::cerr << 2.0 * CGS_Constants::mu * fd.number_density(1.02086,T) << '\n';
  */

#if 1
  ChargeEquilFunc cef(CGS_Constants::me,0.5,1.0e-9,1.0e4);

  size_t Nrhostps = 100;
  double lrho_min = 2.0; // log(rho) in cgs
  double lrho_max = 11.0;
  double lrho_stp = (lrho_max-lrho_min)/(double(Nrhostps)-1);

  size_t NTstps = 8;
  double lT_min = 2.0;
  double lT_max = 9.0;
  double lT_stp = (lT_max - lT_min)/(double(NTstps)-1);
  
  // Temperature in K
  //double Te = 1.0e7;
  
  // Density in cgs
  //double rho = 1122;
  //double rho = 1e5;

  for (size_t j=0; j<NTstps; ++j) {
    for (size_t i=0; i<Nrhostps; ++i) {
      double Te = std::pow(10.0, lT_min + j*lT_stp);
      double rho = std::pow(10.0, lrho_min + i*lrho_stp);
      
      
      double T = Te * (CGS_Constants::k/(CGS_Constants::me*CGS_Constants::c*CGS_Constants::c)); // Temperature in rest mass    
      double n = rho / (2.0*CGS_Constants::mu); // proton number density
      
      cef.set_T(T); // Set temperature
      cef.set_n(n); // Set proton number density
      
#if 0
      for (size_t j=0; j<Nstps; ++j) {
	double mum1 = -1 + 2*j/(double(Nstps)-1);
	std::cout << std::setw(15) << Te
		  << std::setw(15) << mum1
		  << std::setw(15) << cef(mum1)
		  << "\n";
      }
      std::cout << '\n' << std::endl;
#endif
      
#if 1
      // Define root finder
      NumericalRecipes::RootFinder<ChargeEquilFunc> rf(cef);
      
      // Get chemical potential (in units of rest mass)
      double muemin = std::max(-3.0,-1.0/T);
      double muemax = std::min(1.0,2.0/T);

      // Check to see if the root is bracketed
      double fa = cef(muemin);
      double fb = cef(muemax);
      for (size_t nbrac=0; nbrac<1000 && fa*fb>0.0; ++nbrac) {
	muemin = std::max(muemin-1.0,-1.0/T);
	muemax = muemax*2.0;
	fa = cef(muemin);
	fb = cef(muemax);
	/*
	  std::cerr << std::setw(15) << muemin
	  << std::setw(15) << muemax
	  << '\n';
	*/
      }
      
      std::cerr << "Done with bracketing! "
		<< std::setw(15) << fa << " @" << std::setw(15) << muemin
		<< std::setw(15) << fb << " @" << std::setw(15) << muemax
		<< "\n";
      
      double mue = rf.zbrent(muemin,muemax,1.0e-7);
      
      mue = 1 + mue*T;
      
      // Use this to find pressure
      double ne = fd.number_density(mue,T) - fd.number_density(-mue,T);
      double P = fd.pressure(mue,T) + fd.pressure(-mue,T);
      
      // Density from E_f = mue
      double pf = CGS_Constants::me * CGS_Constants::c * std::sqrt( mue*mue - 1.0 );
      if (isnan(pf))
	pf = 0.0;
      double rhof = (8.0*M_PI/3.0) * std::pow( pf/CGS_Constants::h , 3.0 ) * 2.0*CGS_Constants::mu;
      
      // Degeneracy pressure estimates from density
      double Pnr = (8.0*M_PI*CGS_Constants::h*CGS_Constants::h/(15.0*CGS_Constants::me)) * std::pow( 3.0*n/(8.0*M_PI), 5.0/3.0 );
      double Pur = (2.0*M_PI*CGS_Constants::h*CGS_Constants::c/3.0) * std::pow( 3.0*n/(8.0*M_PI), 4.0/3.0 );
      double kf = std::pow(3.0*n/(8.0*M_PI),1.0/3.0) * CGS_Constants::h/(CGS_Constants::me*CGS_Constants::c);
      double ef = std::sqrt(kf*kf+1.0);
      
      double Pdeg = (8*M_PI/3.0) * std::pow(CGS_Constants::me*CGS_Constants::c/CGS_Constants::h,3) * CGS_Constants::me*CGS_Constants::c*CGS_Constants::c
	* (  (1.0/4.0)*kf*ef*ef*ef   -   (5.0/8.0)*kf*ef   +   (3.0/8.0)*std::log(kf+ef)  );
      
      
      // Output
      std::cout << std::setw(15) << Te
		<< std::setw(15) << rho
		<< std::setw(15) << rhof
		<< std::setw(15) << ne * 2.0 * CGS_Constants::mu
		<< std::setw(15) << P
		<< std::setw(15) << Pnr
		<< std::setw(15) << Pur
		<< std::setw(15) << Pdeg
		<< std::setw(15) << mue
		<< std::endl;
#endif
    }
    std::cout << '\n' <<  std::endl;
  }
#endif

  return 0;
}
