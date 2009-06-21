#ifndef NUMERICALRECIPES_H
#define NUMERICALRECIPES_H

#include <iostream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include "constants.h"

class NumericalRecipes
{
 public:
  
  class Exception  // Exception class; contains a simple error message
  {
    public:
      const char *text;
      Exception(const char *msg) : text(msg) { }
  };


  //--- Utility routines and variables ---

  static double sqrarg;
  static double SQR(double a);

  static double dsqrarg;
  static double DSQR(double a);

  static double dmaxarg1,dmaxarg2;
  static double DMAX(double a, double b);

  static double dminarg1,dminarg2;
  static double DMIN(double a, double b);

  static double maxarg1,maxarg2;
  static double FMAX(double a, double b);

  static double minarg1,minarg2;
  static double FMIN(double a, double b);

  static long lmaxarg1,lmaxarg2;
  static long LMAX(long a, long b);

  static long lminarg1,lminarg2;
  static long LMIN(long a, long b);

  static int imaxarg1,imaxarg2;
  static int IMAX(int a, int b);

  static int iminarg1,iminarg2;
  static int IMIN(int a, int b);

  static double SIGN(double a, double b);

  static void nrerror(const char error_text[]);
  static float *vector(long nl, long nh);
  static int *ivector(long nl, long nh);
  static unsigned char *cvector(long nl, long nh);
  static unsigned long *lvector(long nl, long nh);
  static double *dvector(long nl, long nh);
  static float **matrix(long nrl, long nrh, long ncl, long nch);
  static double **dmatrix(long nrl, long nrh, long ncl, long nch);
  static int **imatrix(long nrl, long nrh, long ncl, long nch);
  static float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
  static float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
  static float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  static void free_vector(float *v, long nl, long nh);
  static void free_ivector(int *v, long nl, long nh);
  static void free_cvector(unsigned char *v, long nl, long nh);
  static void free_lvector(unsigned long *v, long nl, long nh);
  static void free_dvector(double *v, long nl, long nh);
  static void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
  static void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
  static void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
  static void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
  static void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
  static void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

  //--- Computational Routines ---
  
  static void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2,
		     double **c);
  static void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
		     double x1u, double x2l, double x2u, double x1, double x2, double *ansy,
		     double *ansy1, double *ansy2);
  static double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  template<class F> static double brent(double ax, double bx, double cx, F f, double tol, double *xmin);
  //static void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps, double yscal[],
  //	     double *hdid, double *hnext, void (*derivs)(double, double [], double []), double hmax=0.0);
  static void convlv(double data[], unsigned long n, double respns[], unsigned long m,
		     int isign, double ans[]);
  static void correl(double data1[], double data2[], unsigned long n, double ans[]);
  static void cosft1(double y[], int n);
  static void cosft2(double y[], int n, int isign);
  static double dfridr(double (*func)(double), double x, double h, double *err);
  static double factrl(int n);
  static void fdjac(int n, double x[], double fvec[], double **df,
		    void (*vecfunc)(int, double [], double []));
  static double fmin(double x[]);
  static void four1(double data[], unsigned long nn, int isign);
  static void fourn(double data[], unsigned long nn[], int ndim, int isign);
  static double gammln(double xx);
  static void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		     double *f, double stpmax, int *check, double (*func)(double []));
  static void lubksb(double **a, int n, int *indx, double b[]);
  static void ludcmp(double **a, int n, int *indx, double *d);
  static void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[],
		   void (*derivs)(double, double[], double[]));
  static void mnbrak(double *ax, double *bx, double *cx,
		     double *fa, double *fb, double *fc,
		     double (*func)(double));
  template<class F> static void mnbrak(double *ax, double *bx, double *cx,
				       double *fa, double *fb, double *fc,
				       F func);
  static void newt(double x[], int n, int *check, void (*vecfunc)(int, double [], double []));
  static void odeint(double ystart[], int nvar, double x1, double x2, double eps, double *h1,
		     double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double []),
		     void (*rkqs)(double [], double [], int, double *, double, double, double [],
				  double *, double *, void (*)(double, double [], double []),
				  double),double hmax=0.0);
  static double plgndr(int l, int m, double x);
  static double pythag(double a, double b);
  static void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
  static void realft(double data[], unsigned long n, int isign);
  static void rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,
		    unsigned long nn3, int isign);
  static void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
		   double yerr[], void (*derivs)(double, double [], double []));
  static void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
		   double yscal[], double *hdid, double *hnext,
		   void (*derivs)(double, double [], double []), double hmax=0.0);
  static double rtflsp(double (*func)(double), double x1, double x2, double xacc);
  static double rtsafe(void (*funcd)(double, double *, double *),
		       double x1, double x2, double xacc);
  static void savgol(double c[], int np, int nl, int nr, int ld, int m);
  static void shootf(int n, double v[], double f[]);
  static void sinft(double y[], int n);
  static std::complex<double> spharm(int l, int m, double theta, double phi);
  static double spharm_e(int l, int m, double theta, double phi);
  static double spharm_o(int l, int m, double theta, double phi);
  static void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
  static void svdcmp(double **a, int m, int n, double w[], double **v);
  static void svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma,
		     double **u, double **v, double w[], double *chisq,
		     void (*funcs)(double, double [], int));
  static void tridag(double a[], double b[], double c[], double r[], double u[], unsigned long n);
  template<typename _T> static void tridag(_T a[], _T b[], _T c[], _T r[], _T u[], int n);
  static void twofft(double data1[], double data2[], double fft1[], double fft2[], unsigned long n);
  static int zbrac(double (*func)(double), double *x1, double *x2);
  static double zbrent(double (*func)(double), double x1, double x2, double tol);

  //--- Base class for function objects to calculate derivatives ---
  class derivs_obj
  {
   public:
    virtual ~derivs_obj() { }
    virtual void operator()(double x, const std::valarray<double>& y, std::valarray<double>& dydx) = 0;
  };

  //--- Integrator class ---
  template <class F>
  class Integrator
  {
   public:
    Integrator(unsigned int n, F deriv);
    
    void odeint(std::valarray<double>& ystart, double x1, double x2, double eps,
                double *h1, double hmin, double hmax=0.0);

    void resize(unsigned int n);

    F derivs;

   private:
    void rkqs(std::valarray<double>& y, std::valarray<double>& dydx,
	      double *x, double htry, double eps, std::valarray<double>& yscal,
	      double *hdid, double *hnext, double hmax=0.0);
    void rkck(std::valarray<double>& y, std::valarray<double>& dydx, double x, double h,
	      std::valarray<double>& yout, std::valarray<double>& yerr);

    long odeint_MAXSTP;
    double odeint_TINY;
    double rkqs_SAFETY, rkqs_PGROW, rkqs_PSHRNK, rkqs_ERRCON;

    std::valarray<double> odeint_yscal, odeint_y, odeint_dydx;
    std::valarray<double> rkqs_yerr, rkqs_ytemp;
    std::valarray<double> rkck_ak2, rkck_ak3, rkck_ak4, rkck_ak5, rkck_ak6, rkck_ytemp;
  };

  //--- RootFinder class ---
  template<class F>
  class RootFinder
  {
   public:
    explicit RootFinder(F fun);
    double rtflsp(double x1, double x2, double xacc);
    double rtsafe(double x1, double x2, double xacc);
    double zbrent(double x1, double x2, double tol);

    F func;

   private:
    bool zbrac(double *x1, double *x2);

    int rtflsp_MAXIT;
    double zbrac_FACTOR;
    int zbrac_NTRY;
    int rtsafe_MAXIT;
    int zbrent_ITMAX;
    double zbrent_EPS;
  };

  //--- Chebyshev approximation of functions ---
  class FunctionApproximation
  {
   public:
    template<class _F>
    FunctionApproximation(double xmin, double xmax, double eps, _F func);

    double operator()(double x) const;

    FunctionApproximation deriv() const;

   private:
    double xmin_, xmax_;
    std::valarray<double> coeff_;

    FunctionApproximation(double xmin, double xmax, std::valarray<double> coeff);
  };

  //--- Tridiagonal solver ---
  template<typename _T>
  class TridiagSolver
  {
   public:
    explicit TridiagSolver(int n);
    ~TridiagSolver();

    void operator()(_T a[], _T b[], _T c[], _T r[], _T u[]);

   private:
    int n_;
    _T* gam;
  };


  //--- Variables and parameters ---

  // brent
  static const int brent_ITMAX;
  static const double brent_CGOLD;
  static const double brent_ZEPS;
  // bsstep
  static const int bsstep_KMAXX;
  static const int bsstep_IMAXX;
  static const double bsstep_SAFE1;
  static const double bsstep_SAFE2;
  static const double bsstep_REDMAX;
  static const double bsstep_REDMIN;
  static const double bsstep_TINY;
  static const double bsstep_SCALMX;
  static double **d,*x;
  // dfridr
  static const double dfridr_CON;
  static const double dfridr_CON2;
  static const double dfridr_BIG;
  static const int dfridr_NTAB;
  static const double dfridr_SAFE;
  // fdjac
  static const double fdjac_EPS;
  // lnsrch
  static const double lnsrch_ALF;
  static const double lnsrch_TOLX;
  // ludcmp
  static const double ludcmp_TINY;
  // mnbrak
  static const double mnbrak_GOLD;
  static const double mnbrak_GLIMIT;
  static const double mnbrak_TINY;
  // newt
  static const int newt_MAXITS;
  static const double newt_TOLF;
  static const double newt_TOLMIN;
  static const double newt_TOLX;
  static const double newt_STPMX;
  static int nn;
  static double *fvec;
  static void (*nrfuncv)(int n, double v[], double f[]);
  // odeint
  static const long odeint_MAXSTP;
  static const double odeint_TINY;
  static int kmax,kount;
  static double *xp,**yp,dxsav;
  // rkqs
  static const double rkqs_SAFETY;
  static const double rkqs_PGROW;
  static const double rkqs_PSHRNK;
  static const double rkqs_ERRCON;
  // rtflsp
  static const int rtflsp_MAXIT;
  // rtsafe
  static const int rtsafe_MAXIT;
  // shootf
  static const double shootf_EPS;
  static int nn2,nvar;
  static double x1,x2,xf;
  static void (*derivs)(double x, double y[], double dydx[]);
  static void (*load1)(double x1, double v1[], double y[]);
  static void (*load2)(double x2, double v2[], double y[]);
  static void (*score)(double xf, double y[], double f[]);
  // svdfit
  static const double svdfit_TOL;
  // zbrac
  static const double zbrac_FACTOR;
  static const int zbrac_NTRY;
  // zbrent
  static const int zbrent_ITMAX;
  static const double zbrent_EPS;
};

inline double NumericalRecipes::SQR(double a)  { return ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg); }

inline double NumericalRecipes::DSQR(double a)  { return ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg); }

inline double NumericalRecipes::DMAX(double a, double b)
{
  return (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2));
}

inline double NumericalRecipes::DMIN(double a, double b)
{
  return (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2));
}

inline double NumericalRecipes::FMAX(double a, double b)
{
  return (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2));
}

inline double NumericalRecipes::FMIN(double a, double b)
{
  return (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2));
}

inline long NumericalRecipes::LMAX(long a, long b)
{
  return (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2));
}

inline long NumericalRecipes::LMIN(long a, long b)
{
  return (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2));
}

inline int NumericalRecipes::IMAX(int a, int b)
{
  return (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2));
}

inline int NumericalRecipes::IMIN(int a, int b)
{
  return (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2));
}

inline double NumericalRecipes::SIGN(double a, double b)  { return ((b) >= 0.0 ? std::fabs(a) : -std::fabs(a)); }


//--- Integrator class ---
template<class F>
NumericalRecipes::Integrator<F>::Integrator(unsigned int n, F deriv)
  : derivs(deriv), odeint_MAXSTP(100000000), odeint_TINY(1.0e-26),
  rkqs_SAFETY(0.9), rkqs_PGROW(-0.2), rkqs_PSHRNK(-0.25), rkqs_ERRCON(1.89e-4),
  odeint_yscal(n), odeint_y(n), odeint_dydx(n), rkqs_yerr(n), rkqs_ytemp(n),
  rkck_ak2(n), rkck_ak3(n), rkck_ak4(n), rkck_ak5(n), rkck_ak6(n), rkck_ytemp(n)
{
}

template<class F>
void NumericalRecipes::Integrator<F>::resize(unsigned int n)
{
  odeint_yscal.resize(n);
  odeint_y.resize(n);
  odeint_dydx.resize(n);
  rkqs_yerr.resize(n);
  rkqs_ytemp.resize(n);
  rkck_ak2.resize(n);
  rkck_ak3.resize(n);
  rkck_ak4.resize(n);
  rkck_ak5.resize(n);
  rkck_ak6.resize(n);
  rkck_ytemp.resize(n);
}

template<class F>
void NumericalRecipes::Integrator<F>::odeint(std::valarray<double>& ystart, double x1, double x2,
					     double eps, double *h1, double hmin, double hmax)
{
  int nstp;
  double x,hnext,hdid,h;
  double odeint_MAXSTP_local = odeint_MAXSTP;
    
  if ( (hmax != 0.0) && (std::fabs(*h1) > std::fabs(hmax)) ) 
    (*h1)=hmax; 
  if ( (hmax != 0.0) && (odeint_MAXSTP < std::abs(int((x2-x1)/hmax))) )
    odeint_MAXSTP_local = 2.0*std::fabs((x2-x1)/hmax);
    
  x=x1;
  h=SIGN(*h1,x2-x1);
  odeint_y = ystart;
  for (nstp=1;nstp<=odeint_MAXSTP_local;++nstp)
    {
      derivs(x,odeint_y,odeint_dydx);
      odeint_yscal = std::abs(odeint_y) + std::abs(odeint_dydx*h) + odeint_TINY;
      if ((x+h-x2)*(x+h-x1) > 0.0)
        h = x2-x;
      rkqs(odeint_y,odeint_dydx,&x,h,eps,odeint_yscal,&hdid,&hnext,hmax);
      if ((x-x2)*(x2-x1) >= 0.0)
	  {
	    ystart = odeint_y;
	    *h1 = h;
	    return;
	  }
      if (std::fabs(hnext) <= hmin)
	    throw NumericalRecipes::Exception("Step size too small in Integrator::odeint");
      h=hnext;
    }
  throw NumericalRecipes::Exception("Too many steps in Integrator::odeint");
}
 
template <class F>
void NumericalRecipes::Integrator<F>::rkqs(std::valarray<double>& y, std::valarray<double>& dydx, 
					   double *x, double htry, double eps,
					   std::valarray<double>& yscal, double *hdid,
					   double *hnext, double hmax)
{
  unsigned int i;
  double errmax,h,htemp,xnew;
      
  h=htry;
  for (;;)
  {
    rkck(y,dydx,*x,h,rkqs_ytemp,rkqs_yerr);
    errmax = 0.0;
    for (i=0;i<y.size();++i)
      errmax = std::max(errmax,std::fabs(rkqs_yerr[i]/yscal[i]));
    errmax /= eps;
    if 
      (errmax <= 1.0) break;
    htemp = rkqs_SAFETY*h*std::pow(errmax,rkqs_PSHRNK);
    h = (h >= 0.0 ? std::max(htemp,0.1*h) : std::min(htemp,0.1*h));
    xnew = (*x) + h;
    if (xnew == *x)
      throw NumericalRecipes::Exception("stepsize underflow in Integrator::rkqs");
  }
  if (errmax > rkqs_ERRCON)
    *hnext = rkqs_SAFETY*h*std::pow(errmax,rkqs_PGROW);
  else 
    *hnext = 5.0*h;
  if ( (hmax != 0.0) && (std::fabs(hmax/(*hnext)) < 1.0) )
    (*hnext) = (*hnext)*std::fabs(hmax/(*hnext));
  *x += (*hdid=h);
  y = rkqs_ytemp;
}
  
template <class F>
void NumericalRecipes::Integrator<F>::rkck(std::valarray<double>& y, std::valarray<double>& dydx,
					   double x, double h, std::valarray<double>& yout,
					   std::valarray<double>& yerr)
{
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  static double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
      
  rkck_ytemp = y + b21*h*dydx;
  derivs(x+a2*h,rkck_ytemp,rkck_ak2);
  rkck_ytemp = y + h*(b31*dydx + b32*rkck_ak2);
  derivs(x+a3*h,rkck_ytemp,rkck_ak3);
  rkck_ytemp = y + h*(b41*dydx + b42*rkck_ak2 +b43*rkck_ak3);
  derivs(x+a4*h,rkck_ytemp,rkck_ak4);
  rkck_ytemp = y + h*(b51*dydx + b52*rkck_ak2 + b53*rkck_ak3 + b54*rkck_ak4);
  derivs(x+a5*h,rkck_ytemp,rkck_ak5);
  rkck_ytemp = y + h*(b61*dydx + b62*rkck_ak2 + b63*rkck_ak3 + b64*rkck_ak4 + b65*rkck_ak5);
  derivs(x+a6*h,rkck_ytemp,rkck_ak6);
  yout = y + h*(c1*dydx + c3*rkck_ak3 + c4*rkck_ak4 + c6*rkck_ak6);
  yerr = h*(dc1*dydx + dc3*rkck_ak3 + dc4*rkck_ak4 + dc5*rkck_ak5 + dc6*rkck_ak6);
}


//--- RootFinder class ---
template<class F>
NumericalRecipes::RootFinder<F>::RootFinder(F fun)
  :  func(fun), rtflsp_MAXIT(100000), zbrac_FACTOR(1.6), zbrac_NTRY(2000),
     rtsafe_MAXIT(1000), zbrent_ITMAX(2000), zbrent_EPS(1.0e-8)
{
}

template <class F>
double NumericalRecipes::RootFinder<F>::rtflsp(double x1, double x2, double xacc)
{
  int j;
  double fl,fh,xl,xh,swap,dx,del,f,rtf;

  fl = func(x1);
  fh = func(x2);

  if (fl*fh > 0.0)
  {
    if (zbrac(&x1,&x2))
    {
      fl = func(x1);
      fh = func(x2);
    }
    else {
      std::cerr << "Failed to bracket root in RootFinder::rtflsp\n";
      throw NumericalRecipes::Exception("Failed to bracket root in RootFinder::rtflsp");
    }
  }
  if (fl < 0.0)
  {
    xl = x1;
    xh = x2;
  }
  else
  {
    xl = x2;
    xh = x1;
    swap = fl;
    fl = fh;
    fh = swap;
  }
  dx = xh - xl;
  for (j=1; j<=rtflsp_MAXIT; ++j)
  {
    rtf = xl+dx*fl/(fl-fh);

    f = func(rtf);

    if (f < 0.0)
    {
      del = xl-rtf;
      xl = rtf;
      fl = f;
    }
    else
    {
      del = xh-rtf;
      xh = rtf;
      fh = f;
    }
    dx = xh-xl;
    if (std::fabs(del) < xacc || f == 0.0)
      return rtf;
  }
  std::cerr << "Exceeded maximum iterations in RootFinder::rtflsp\n";
  std::cerr << "Max iterations: " << fl << " @ " << xl << '\t' << fh << " @ " << xh << '\n';
  return ((fl*xh+fh*xl)/(fl+fh));
  //  throw NumericalRecipes::Exception("Exceeded maximum iterations in RootFinder::rtflsp");
  //  return 0.0;
}

template<class F>
bool NumericalRecipes::RootFinder<F>::zbrac(double *x1, double *x2)
{
  int j;
  double f1,f2;

  if (*x1 == *x2)
    throw NumericalRecipes::Exception("Bad initial range in RootFinder::zbrac");
  f1 = func(*x1);
  f2 = func(*x2);
  for (j=1; j<=zbrac_NTRY; ++j)
  {
    if (f1*f2 < 0.0)
      return true;
    if (std::fabs(f1) < std::fabs(f2))
      f1 = func(*x1 += zbrac_FACTOR*(*x1-*x2));
    else
      f2 = func(*x2 += zbrac_FACTOR*(*x2-*x1));
  }
  return false;
}

template<class F>
double NumericalRecipes::RootFinder<F>::rtsafe(double x1, double x2, double xacc)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  func(x1,&fl,&df);
  func(x2,&fh,&df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    throw NumericalRecipes::Exception("Root must be bracketed in Rootfinder::rtsafe");
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=std::fabs(x2-x1);
  dx=dxold;
  func(rts,&f,&df);
  for (j=1;j<=rtsafe_MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
	|| (std::fabs(2.0*f) > std::fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (std::fabs(dx) < xacc) return rts;
    func(rts,&f,&df);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  throw NumericalRecipes::Exception("Maximum number of iterations exceeded in Rootfinder::rtsafe");
  return 0.0;
}

template<class F>
double NumericalRecipes::RootFinder<F>::zbrent(double x1, double x2, double tol)
{
  int iter;
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
  double fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    std::cerr << "Root must be bracketed in Rootfinder::zbrent\n";
    std::cerr << "Bracket: "
	      << std::setw(15) << fa << " @ " << std::setw(15) << a
	      << std::setw(15) << fb << " @ " << std::setw(15) << b
	      << '\n';
    throw NumericalRecipes::Exception("Root must be bracketed in Rootfinder::zbrent");
  }
  fc=fb;
  for (iter=1;iter<=zbrent_ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (std::fabs(fc) < std::fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*zbrent_EPS*std::fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (std::fabs(fb) <= tol || fb == 0.0) return b;
    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=std::fabs(p);
      min1=3.0*xm*q-std::fabs(tol1*q);
      min2=std::fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (std::fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=func(b);
  }
  std::cerr << "Max iterations: "
	    << std::setw(15) << fa << " @ " << std::setw(15) << a
	    << "     "
	    << std::setw(15) << fb << " @ " << std::setw(15) << b << '\n';
  //throw NumericalRecipes::Exception("Maximum number of iterations exceeded in RootFinder::zbrent");
  //return 0.0;
  return ((fb*a+fa*b)/(fa+fb));
}

//--- Chebyshev function approximation ---
template<class _F>
NumericalRecipes::FunctionApproximation::FunctionApproximation(double xmin, double xmax,
							       double eps, _F func)
  : xmin_(xmin), xmax_(xmax), coeff_(50)
{
  size_t n = coeff_.size();
  double err;
  int iter = 1;
  do
  {
    double* f = new double[n];
    double bma = 0.5*(xmax - xmin);
    double bpa = 0.5*(xmax + xmin);
    for (size_t k=0; k<n; ++k)
    {
      double y = std::cos(M_PI*(k+0.5)/n);
      f[k] = func(y*bma+bpa);
    }
    double fac = 2.0/n;
    for (size_t j=0; j<n; ++j)
    {
      double sum=0.0;
      for (size_t k=0; k<n; ++k)
	sum += f[k]*std::cos(M_PI*j*(k+0.5)/n);
      coeff_[j] = fac*sum;
    }
    delete[] f;

    double x;
    err = 0.0;
    for (size_t k=0; k<n; ++k)
    {
      x = xmin + (xmax-xmin)*k/(n-1.0);
      err += std::pow(func(x)-(*this)(x),2);
    }
    err = std::sqrt(err/n);

    if (err > eps)
    {
      n += 10*(1 + iter>>1);
      coeff_.resize(n);
    }

    ++iter;
  } while (err>eps);
}

inline NumericalRecipes::FunctionApproximation::FunctionApproximation(double xmin, double xmax,
							       std::valarray<double> coeff)
  : xmin_(xmin), xmax_(xmax), coeff_(coeff)
{
}

inline double NumericalRecipes::FunctionApproximation::operator()(double x) const
{
  if ((x-xmin_)*(x-xmax_) > 0.0)
    throw NumericalRecipes::Exception("x not in range in FunctionApproximation::operator()");
  double y = (2.0*x-xmin_-xmax_)/(xmax_-xmin_);
  double y2 = 2.0*y;
  double d=0.0, dd=0.0, sv;
  int m = coeff_.size();
  for (int j=m-1; j>=1; --j)
  {
    sv = d;
    d = y2*d-dd+coeff_[j];
    dd = sv;
  }
  return y*d-dd+0.5*coeff_[0];
}

inline NumericalRecipes::FunctionApproximation NumericalRecipes::FunctionApproximation::deriv() const
{
  int n = coeff_.size();
  std::valarray<double> coeff_deriv(n);
  coeff_deriv[n-1] = 0.0;
  coeff_deriv[n-2] = 2*(n-1)*coeff_[n-1];
  for (int j=n-3; j>=0; --j)
    coeff_deriv[j] = coeff_deriv[j+2]+2*(j+1)*coeff_[j+1];
  double con = 2.0/(xmax_-xmin_);
  for (int j=0; j<n; ++j)
    coeff_deriv[j] *= con;
  return FunctionApproximation(xmin_,xmax_,coeff_deriv);
}

//--- Tridiagonal solver ---
template<typename _T>
NumericalRecipes::TridiagSolver<_T>::TridiagSolver(int n)
 : n_(n)
{
  gam = new _T[n];
}

template<typename _T>
NumericalRecipes::TridiagSolver<_T>::~TridiagSolver()
{
  delete[] gam;
}

template<typename _T>
void NumericalRecipes::TridiagSolver<_T>::operator()(_T a[], _T b[], _T c[], _T r[], _T u[])
{
  int j;
  _T bet;
 
  if (b[0] == 0.0)
    nrerror("Error 1 in tridag");
  u[0] = r[0]/(bet=b[0]);
  for (j=1; j<n_; ++j)
  {
    gam[j] = c[j-1]/bet;
    bet = b[j] - a[j]*gam[j];
    if (bet == 0.0)
      nrerror("Error 2 in tridag");
    u[j] = (r[j] - a[j]*u[j-1])/bet;
  }
  for (j=(n_-2); j>=0; --j)
    u[j] -= gam[j+1]*u[j+1];
}


//--- Template versions of routines ---

//--- brent ---
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
template<class F>
double NumericalRecipes::brent(double ax, double bx, double cx, F f, double tol,
			       double *xmin)
{
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=f(x);
  for (iter=1;iter<=brent_ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*std::fabs(x)+brent_ZEPS);
    if (std::fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (std::fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=std::fabs(q);
      etemp=e;
      e=d;
      if (std::fabs(p) >= std::fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=brent_CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=brent_CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(std::fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=f(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }
  nrerror("Too many iterations in brent");
  *xmin=x;
  return fx;
}
#undef SHFT

//--- mnbrak ---
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
template<class F>
void NumericalRecipes::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
			      F func)
{
  double ulim,u,r,q,fu,dum;

  *fa=func(*ax);
  *fb=func(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+mnbrak_GOLD*(*bx-*ax);
  *fc=func(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(std::fabs(q-r),mnbrak_TINY),q-r));
    ulim=(*bx)+mnbrak_GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=func(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+mnbrak_GOLD*(*cx-*bx);
      fu=func(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=func(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+mnbrak_GOLD*(*cx-*bx))
	SHFT(*fb,*fc,fu,func(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=func(u);
    } else {
      u=(*cx)+mnbrak_GOLD*(*cx-*bx);
      fu=func(u);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}
#undef SHFT

// --- tridag ---
template<typename _T>
void NumericalRecipes::tridag(_T a[], _T b[], _T c[], _T r[], _T u[], int n)
{
  int j;
  _T bet;

  _T* gam = new _T[n];
  if (b[0] == 0.0)
    nrerror("Error 1 in tridag");
  u[0] = r[0]/(bet=b[0]);
  for (j=1; j<n; ++j)
  {
    gam[j] = c[j-1]/bet;
    bet = b[j] - a[j]*gam[j];
    if (bet == 0.0)
      nrerror("Error 2 in tridag");
    u[j] = (r[j] - a[j]*u[j-1])/bet;
  }
  for (j=(n-2); j>=0; --j)
    u[j] -= gam[j+1]*u[j+1];
  delete[] gam;
}


#endif
