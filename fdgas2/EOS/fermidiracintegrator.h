/**************************************************************/
/*** TAKES CARE OF INTEGRATING THE FERMI-DIRAC DISTRIBUTION ***/
/*                                                            */
/* This takes care of making sure that the integrals are      */
/* properly done near the transition in the Fermi-Dirac       */
/* distribution function.                                     */
/*                                                            */
/* Expects a function f(E) to be integrated (note that this   */
/* is a function of energy!, not x)                           */
/*                                                            */
/**************************************************************/

#ifndef FERMI_DIRAC_INTEGRATOR_H
#define FERMI_DIRAC_INTEGRATOR_H

#include <algorithm>
#include <valarray>
#include <cmath>

// Function Class Definition
class Func {
 public:
  virtual double operator()(double) = 0;
};

#define NDIM_INTS 1 // Number of integrals

class FermiDiracIntegrator
{
 public:
  double integrate(Func& f, double mu, double T, double eps=1.0e-10, double infty=1.0e2);

 private:
  Func* _f; // Function of energy in UNITS OF THE REST MASS!!!!
  double _T, _mu; // Temperature and chemical potential in units of the rest mass!

  // The Fermi-Dirac distribution function in terms of x=(E-mu)/T
  double fermi_dirac_dist(double x);

  // The integrand to be integrated in terms of x=(E-mu)/T (the integral is y[0])
  inline void derivs(double x, double y[], double dydx[]);

  // Error scales
  inline void get_yscal(double,double,double[],double[],double[]);

  // Numerics
  int rkqs(double[],double[],int,double&,double,double,double[],double&,double&);
  void rkck(double [],double [],int,double,double,double [],double []);
};


double FermiDiracIntegrator::fermi_dirac_dist(double x)
{
  return 1.0/(1.0 + std::exp(x));
}

#define NSTPS_MAX 100000
double FermiDiracIntegrator::integrate(Func& f, double mu, double T, double eps, double infty)
{
  _f = &f;
  _mu = mu;
  _T = T;

  double x = (1.0-_mu)/(_T);   // Initial (E-mu)/T
  double dydx[1], y[] = {0.0}; // Initial integral value
  double yscal[1];             // Error scale
  double hnext, hdid, h;       // Some stepsizes

  double xstepping = 10.0;     // Region to carefully step through
  double hstepping = 1.0;      // Step size in carefully-stepping region

  h = std::max(std::fabs(x)/xstepping,hstepping);

  for (size_t nstps=0; x<infty && nstps<NSTPS_MAX; ++nstps) {

    // Get derivatives for first rkqs
    derivs(x,y,dydx);

    // Get yscal for rkqs
    get_yscal(h,x,y,dydx,yscal);

    // Check if we are entering carefully-stepping region
    if ( std::min(std::fabs(x),std::fabs(x+h)) < xstepping ) {
      if ( x < -xstepping ) // If first entering
	h = -xstepping - x + hstepping; // let step into gingerly
      else
	h = std::min(h,hstepping);  // reduce to proper stepsize if necessary
    }

    /*
    // DEBUGGING
    std::cout << std::setw(15) << x
	      << std::setw(15) << _T*x+_mu-1
	      << std::setw(15) << _T*y[0]
	      << std::setw(15) << dydx[0]
	      << std::setw(15) << fermi_dirac_dist(x)
	      << std::setw(15) << yscal[0]
	      << std::endl;
    */

    // Take rkqs
    rkqs(y,dydx,NDIM_INTS,x,h,eps,yscal,hdid,hnext);

    // Set next stepsize
    h = hnext;      
  }

  // DEBUGGING
  //std::cout << '\n' << std::endl;

  return _T*y[0];
}
#undef NSTPS_MAX


void FermiDiracIntegrator::derivs(double x, double y[], double dydx[])
{
  double E = _T*x + _mu;
  //        function * Fermi-Dirac dist. * integration vol. elements
  dydx[0] =   (*_f)(E)  * fermi_dirac_dist(x) * E * (E>1.0 ? std::sqrt(E*E-1.0) : 0.0);
}

#define TINY 1.0e-6;
void FermiDiracIntegrator::get_yscal(double h, double x, double y[], double dydx[], double yscal[])
{
  yscal[0] = std::fabs(y[0]) + std::fabs(dydx[0]*h) + TINY;
}
#undef TINY


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
int FermiDiracIntegrator::rkqs(double y[], double dydx[], int n, double& x, double htry, 
			       double eps, double yscal[], double& hdid, double& hnext)
{
  static int i;
  static double errmax,h,htemp,xnew,yerr[NDIM_INTS],ytemp[NDIM_INTS]; 
  
  h=htry;
  for (;;) {
    rkck(y,dydx,n,x,h,ytemp,yerr);
    errmax=0.0;
    for (i=0;i<n;i++){
      if (std::fabs(yerr[i]/yscal[i])>errmax)
	errmax = std::fabs(yerr[i]/yscal[i]);
      if (std::isnan(ytemp[i]))
	errmax = 10.0*eps;
    }
    errmax /= eps;

    if (errmax <= 1.0)
      break;
    htemp= ( isinf(errmax) ? 0.5*h : SAFETY*h*std::pow(errmax,PSHRNK) );
    htemp = (h >= 0.0 ? std::max(htemp,0.1*h) : std::min(htemp,0.1*h));
    xnew=x+htemp;
    if (xnew == x){
      return 1;
    }
    h = htemp;
  }

  if (errmax > ERRCON)
    hnext=SAFETY*h*std::pow(errmax,PGROW);
  else
    hnext=5.0*h;
  x += (hdid=h);

  for (i=0;i<n;i++)
    y[i]=ytemp[i];

  return 0;
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void FermiDiracIntegrator::rkck(double y[], double dydx[], int n, double x, double h, 
				double yout[], double yerr[])
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  static double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  static double ak2[NDIM_INTS],ak3[NDIM_INTS],ak4[NDIM_INTS],
    ak5[NDIM_INTS],ak6[NDIM_INTS],ytemp[NDIM_INTS];

  for (i=0;i<n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  derivs(x+a2*h,ytemp,ak2);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  derivs(x+a3*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  derivs(x+a4*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  derivs(x+a5*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  derivs(x+a6*h,ytemp,ak6);
  for (i=0;i<n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}
#endif
