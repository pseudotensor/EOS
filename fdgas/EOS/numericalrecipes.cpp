#include "numericalrecipes.h"

#define NR_END 1
#define FREE_ARG void*

void NumericalRecipes::nrerror(const char error_text[])
  /* Numerical Recipes standard error handler */
{
  throw Exception(error_text);
}

float *NumericalRecipes::vector(long nl, long nh)
  /* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)std::malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

int *NumericalRecipes::ivector(long nl, long nh)
  /* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)std::malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

unsigned char *NumericalRecipes::cvector(long nl, long nh)
  /* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v=(unsigned char *)std::malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

unsigned long *NumericalRecipes::lvector(long nl, long nh)
  /* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v=(unsigned long *)std::malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

double *NumericalRecipes::dvector(long nl, long nh)
  /* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)std::malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}

float **NumericalRecipes::matrix(long nrl, long nrh, long ncl, long nch)
  /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) std::malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) std::malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **NumericalRecipes::dmatrix(long nrl, long nrh, long ncl, long nch)
  /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) std::malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) std::malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **NumericalRecipes::imatrix(long nrl, long nrh, long ncl, long nch)
  /* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) std::malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) std::malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **NumericalRecipes::submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
				    long newrl, long newcl)
  /* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  /* allocate array of pointers to rows */
  m=(float **) std::malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **NumericalRecipes::convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
  /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) std::malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***NumericalRecipes::f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
  /* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t=(float ***) std::malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) std::malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) std::malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void NumericalRecipes::free_vector(float *v, long nl, long nh)
  /* free a float vector allocated with vector() */
{
  std::free((FREE_ARG) (v+nl-NR_END));
}

void NumericalRecipes::free_ivector(int *v, long nl, long nh)
  /* free an int vector allocated with ivector() */
{
  std::free((FREE_ARG) (v+nl-NR_END));
}

void NumericalRecipes::free_cvector(unsigned char *v, long nl, long nh)
  /* free an unsigned char vector allocated with cvector() */
{
  std::free((FREE_ARG) (v+nl-NR_END));
}

void NumericalRecipes::free_lvector(unsigned long *v, long nl, long nh)
  /* free an unsigned long vector allocated with lvector() */
{
  std::free((FREE_ARG) (v+nl-NR_END));
}

void NumericalRecipes::free_dvector(double *v, long nl, long nh)
  /* free a double vector allocated with dvector() */
{
  std::free((FREE_ARG) (v+nl-NR_END));
}

void NumericalRecipes::free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
  /* free a float matrix allocated by matrix() */
{
  std::free((FREE_ARG) (m[nrl]+ncl-NR_END));
  std::free((FREE_ARG) (m+nrl-NR_END));
}

void NumericalRecipes::free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
  /* free a double matrix allocated by dmatrix() */
{
  std::free((FREE_ARG) (m[nrl]+ncl-NR_END));
  std::free((FREE_ARG) (m+nrl-NR_END));
}

void NumericalRecipes::free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
  /* free an int matrix allocated by imatrix() */
{
  std::free((FREE_ARG) (m[nrl]+ncl-NR_END));
  std::free((FREE_ARG) (m+nrl-NR_END));
}

void NumericalRecipes::free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
  /* free a submatrix allocated by submatrix() */
{
  std::free((FREE_ARG) (b+nrl-NR_END));
}

void NumericalRecipes::free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
  /* free a matrix allocated by convert_matrix() */
{
  std::free((FREE_ARG) (b+nrl-NR_END));
}

void NumericalRecipes::free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
				     long ndl, long ndh)
  /* free a float f3tensor allocated by f3tensor() */
{
  std::free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  std::free((FREE_ARG) (t[nrl]+ncl-NR_END));
  std::free((FREE_ARG) (t+nrl-NR_END));
}

double NumericalRecipes::sqrarg = 0.0;
double NumericalRecipes::dsqrarg = 0.0;
double NumericalRecipes::maxarg1 = 0.0;
double NumericalRecipes::maxarg2 = 0.0;
double NumericalRecipes::minarg1 = 0.0;
double NumericalRecipes::minarg2 = 0.0;
double NumericalRecipes::dmaxarg1 = 0.0;
double NumericalRecipes::dmaxarg2 = 0.0;
double NumericalRecipes::dminarg1 = 0.0;
double NumericalRecipes::dminarg2 = 0.0;
long NumericalRecipes::lmaxarg1 = 0;
long NumericalRecipes::lmaxarg2 = 0;
long NumericalRecipes::lminarg1 = 0;
long NumericalRecipes::lminarg2 = 0;
int NumericalRecipes::imaxarg1 = 0;
int NumericalRecipes::imaxarg2 = 0;
int NumericalRecipes::iminarg1 = 0;
int NumericalRecipes::iminarg2 = 0;

//--- Static initializations ---
// brent
const int NumericalRecipes::brent_ITMAX = 100;
const double NumericalRecipes::brent_CGOLD = 0.3819660;
const double NumericalRecipes::brent_ZEPS = 1.0e-10;
// bsstep
const int NumericalRecipes::bsstep_KMAXX = 8;
const int NumericalRecipes::bsstep_IMAXX = NumericalRecipes::bsstep_KMAXX+1;
const double NumericalRecipes::bsstep_SAFE1 = 0.25;
const double NumericalRecipes::bsstep_SAFE2 = 0.7;
const double NumericalRecipes::bsstep_REDMAX = 1.0e-5;
const double NumericalRecipes::bsstep_REDMIN = 0.7;
const double NumericalRecipes::bsstep_TINY = 1.0e-30;
const double NumericalRecipes::bsstep_SCALMX = 0.1;
// dfridr
const double NumericalRecipes::dfridr_CON = 1.4;
const double NumericalRecipes::dfridr_CON2 = (NumericalRecipes::dfridr_CON*NumericalRecipes::dfridr_CON);
const double NumericalRecipes::dfridr_BIG = 1.0e30;
const int NumericalRecipes::dfridr_NTAB = 10;
const double NumericalRecipes::dfridr_SAFE = 2.0;
// fdjac
const double NumericalRecipes::fdjac_EPS = 1.0e-4;
// lnsrch
const double NumericalRecipes::lnsrch_ALF = 1.0e-4;
const double NumericalRecipes::lnsrch_TOLX = 1.0e-8;
// ludcmp
const double NumericalRecipes::ludcmp_TINY = 1.0e-20;
// mnbrak
const double NumericalRecipes::mnbrak_GOLD = 1.618034;
const double NumericalRecipes::mnbrak_GLIMIT  = 100.0;
const double NumericalRecipes::mnbrak_TINY = 1.0e-20;
// newt
const int NumericalRecipes::newt_MAXITS = 500;
const double NumericalRecipes::newt_TOLF = 1.0e-8;
const double NumericalRecipes::newt_TOLMIN = 1.0e-9;
const double NumericalRecipes::newt_TOLX = 1.0e-10;
const double NumericalRecipes::newt_STPMX = 100.0;
// odeint
const long NumericalRecipes::odeint_MAXSTP = 100000000;
const double NumericalRecipes::odeint_TINY = 1.0e-26;
// rkqs
const double NumericalRecipes::rkqs_SAFETY = 0.9;
const double NumericalRecipes::rkqs_PGROW = -0.2;
const double NumericalRecipes::rkqs_PSHRNK = -0.25;
const double NumericalRecipes::rkqs_ERRCON = 1.89e-4;
// rtflsp
const int NumericalRecipes::rtflsp_MAXIT = 10000;
// rtsafe
const int NumericalRecipes::rtsafe_MAXIT = 10000;
// shootf
const double NumericalRecipes::shootf_EPS = 1.0e-8;
// svdfit
const double NumericalRecipes::svdfit_TOL = 1.0e-8;
// zbrac
const double NumericalRecipes::zbrac_FACTOR = 1.6;
const int NumericalRecipes::zbrac_NTRY = 2000;
// zbrent
const int NumericalRecipes::zbrent_ITMAX = 500;
const double NumericalRecipes::zbrent_EPS = 1.0e-8;

//--- bcucof ---
void NumericalRecipes::bcucof(double y[], double y1[], double y2[], double y12[],
			      double d1, double d2, double **c)
{
  static int wt[16][16]=
  { {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
    {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
    {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
    {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
    {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
    {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
    {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
    {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
    {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1} };
  int l,k,j,i;
  double xx,d1d2,cl[16],x[16];

  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}


//--- bcuint ---
void NumericalRecipes::bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
			      double x1u, double x2l, double x2u, double x1, double x2,
			      double *ansy, double *ansy1, double *ansy2)
{
  int i;
  double t,u,d1,d2,**c;

  c=dmatrix(1,4,1,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  *ansy=(*ansy2)=(*ansy1)=0.0;
  for (i=4;i>=1;i--) {
    *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
  free_dmatrix(c,1,4,1,4);
}

//--- brent ---
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
double NumericalRecipes::brent(double ax, double bx, double cx, double (*f)(double), double tol,
			       double *xmin)
{
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
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
    fu=(*f)(u);
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

//--- bsstep ---
double** NumericalRecipes::d = 0;
double* NumericalRecipes::x = 0;
/*void NumericalRecipes::bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps, double yscal[],
	    double *hdid, double *hnext, void (*derivs)(double, double [], double []), double hmax)
{
  int i,iq,k,kk,km=0;
  static int first=1,kmax,kopt;
  static double epsold = -1.0,xnew;
  double eps1,errmax=0.0,fact,h,red=0.0,scale=0.0,work,wrkmin,xest;
  double *err,*yerr,*ysav,*yseq;
  static double a[bsstep_IMAXX+1];
  static double alf[bsstep_KMAXX+1][bsstep_KMAXX+1];
  static int nseq[bsstep_IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
  int reduct,exitflag=0;

  d=dmatrix(1,nv,1,bsstep_KMAXX);
  err=dvector(1,bsstep_KMAXX);
  x=dvector(1,bsstep_KMAXX);
  yerr=dvector(1,nv);
  ysav=dvector(1,nv);
  yseq=dvector(1,nv);
  if (eps != epsold) {
    *hnext = xnew = -1.0e29;
    eps1=bsstep_SAFE1*eps;
    a[1]=nseq[1]+1;
    for (k=1;k<=bsstep_KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
    for (iq=2;iq<=bsstep_KMAXX;iq++) {
      for (k=1;k<iq;k++)
	alf[k][iq]=std::pow(eps1,(a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.0)*(2*k+1)));
    }
    epsold=eps;
    for (kopt=2;kopt<bsstep_KMAXX;kopt++)
      if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
    kmax=kopt;
  }
  h=htry;
  for (i=1;i<=nv;i++) ysav[i]=y[i];
  if (*xx != xnew || h != (*hnext)) {
    first=1;
    kopt=kmax;
  }
  reduct=0;
  for (;;) {
    for (k=1;k<=kmax;k++) {
      xnew=(*xx)+h;
      if (xnew == (*xx)) nrerror("step size underflow in bsstep");
      mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
      xest=SQR(h/nseq[k]);
      pzextr(k,xest,yseq,y,yerr,nv);
      if (k != 1) {
	errmax=bsstep_TINY;
	for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	km=k-1;
	err[km]=std::pow(errmax/bsstep_SAFE1,1.0/(2*km+1));
      }
      if (k != 1 && (k >= kopt-1 || first)) {
	if (errmax < 1.0) {
	  exitflag=1;
	  break;
	}
	if (k == kmax || k == kopt+1) {
	  red=bsstep_SAFE2/err[km];
	  break;
	}
	else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
	  red=1.0/err[km];
	  break;
	}
	else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
	  red=alf[km][kmax-1]*bsstep_SAFE2/err[km];
	  break;
	}
	else if (alf[km][kopt] < err[km]) {
	  red=alf[km][kopt-1]/err[km];
	  break;
	}
      }
    }
    if (exitflag) break;
    red=FMIN(red,bsstep_REDMIN);
    red=FMAX(red,bsstep_REDMAX);
    h *= red;
    reduct=1;
  }
  *xx=xnew;
  *hdid=h;
  first=0;
  wrkmin=1.0e35;
  for (kk=1;kk<=km;kk++) {
    fact=FMAX(err[kk],bsstep_SCALMX);
    work=fact*a[kk+1];
    if (work < wrkmin) {
      scale=fact;
      wrkmin=work;
      kopt=kk+1;
    }
  }
  *hnext=h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact=FMAX(scale/alf[kopt-1][kopt],bsstep_SCALMX);
    if (a[kopt+1]*fact <= wrkmin) {
      *hnext=h/fact;
      kopt++;
    }
  }
  free_dvector(yseq,1,nv);
  free_dvector(ysav,1,nv);
  free_dvector(yerr,1,nv);
  free_dvector(x,1,bsstep_KMAXX);
  free_dvector(err,1,bsstep_KMAXX);
  free_dmatrix(d,1,nv,1,bsstep_KMAXX);
}*/

//--- convlv ---
void NumericalRecipes::convlv(double data[], unsigned long n, double respns[],
			      unsigned long m, int isign, double ans[])
{
  unsigned long i,no2;
  double dum,mag2,*fft;

  fft=dvector(1,n<<1);
  for (i=1;i<=(m-1)/2;i++)
    respns[n+1-i]=respns[m+1-i];
  for (i=(m+3)/2;i<=n-(m-1)/2;i++)
    respns[i]=0.0;
  twofft(data,respns,fft,ans,n);
  no2=n>>1;
  for (i=2;i<=n+2;i+=2) {
    if (isign == 1) {
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
      ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
    } else if (isign == -1) {
      if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
	nrerror("Deconvolving at response zero in convlv");
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
      ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
    } else nrerror("No meaning for isign in convlv");
  }
  ans[2]=ans[n+1];
  realft(ans,n,-1);
  free_dvector(fft,1,n<<1);
}

//--- correl ---
void NumericalRecipes::correl(double data1[], double data2[], unsigned long n, double ans[])
{
  unsigned long no2,i;
  double dum,*fft;

  fft=dvector(1,n<<1);
  twofft(data1,data2,fft,ans,n);
  no2=n>>1;
  for (i=2;i<=n+2;i+=2) {
    ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
    ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  realft(ans,n,-1);
  free_dvector(fft,1,n<<1);
}

//--- cosft1 ---
#define PI 3.141592653589793
void NumericalRecipes::cosft1(double y[], int n)
{
  int j,n2;
  double sum,y1,y2;
  double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;

  theta=PI/n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  sum=0.5*(y[1]-y[n+1]);
  y[1]=0.5*(y[1]+y[n+1]);
  n2=n+2;
  for (j=2;j<=(n>>1);j++) {
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
    y1=0.5*(y[j]+y[n2-j]);
    y2=(y[j]-y[n2-j]);
    y[j]=y1-wi*y2;
    y[n2-j]=y1+wi*y2;
    sum += wr*y2;
  }
  realft(y,n,1);
  y[n+1]=y[2];
  y[2]=sum;
  for (j=4;j<=n;j+=2) {
    sum += y[j];
    y[j]=sum;
  }
}
#undef PI

//--- cosft2 ---
#define PI 3.141592653589793
void NumericalRecipes::cosft2(double y[], int n, int isign)
{
  int i;
  double sum,sum1,y1,y2,ytemp;
  double theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;

  theta=0.5*PI/n;
  wr1=cos(theta);
  wi1=sin(theta);
  wpr = -2.0*wi1*wi1;
  wpi=sin(2.0*theta);
  if (isign == 1) {
    for (i=1;i<=n/2;i++) {
      y1=0.5*(y[i]+y[n-i+1]);
      y2=wi1*(y[i]-y[n-i+1]);
      y[i]=y1+y2;
      y[n-i+1]=y1-y2;
      wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
      wi1=wi1*wpr+wtemp*wpi+wi1;
    }
    realft(y,n,1);
    for (i=3;i<=n;i+=2) {
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=y[i]*wr-y[i+1]*wi;
      y2=y[i+1]*wr+y[i]*wi;
      y[i]=y1;
      y[i+1]=y2;
    }
    sum=0.5*y[2];
    for (i=n;i>=2;i-=2) {
      sum1=sum;
      sum += y[i];
      y[i]=sum1;
    }
  } else if (isign == -1) {
    ytemp=y[n];
    for (i=n;i>=4;i-=2) y[i]=y[i-2]-y[i];
    y[2]=2.0*ytemp;
    for (i=3;i<=n;i+=2) {
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=y[i]*wr+y[i+1]*wi;
      y2=y[i+1]*wr-y[i]*wi;
      y[i]=y1;
      y[i+1]=y2;
    }
    realft(y,n,-1);
    for (i=1;i<=n/2;i++) {
      y1=y[i]+y[n-i+1];
      y2=(0.5/wi1)*(y[i]-y[n-i+1]);
      y[i]=0.5*(y1+y2);
      y[n-i+1]=0.5*(y1-y2);
      wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
      wi1=wi1*wpr+wtemp*wpi+wi1;
    }
  }
}
#undef PI

//--- dfridf ---
double NumericalRecipes::dfridr(double (*func)(double), double x, double h, double *err)
{
    int i,j;
    double errt,fac,hh,**a,ans=0.0;

    if (h == 0.0) nrerror("h must be nonzero in dfridr");
    a=dmatrix(1,dfridr_NTAB,1,dfridr_NTAB);
    hh=h;
    a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
    *err=dfridr_BIG;
    for (i=2;i<=dfridr_NTAB;i++) {
	hh /= dfridr_CON;
	a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
	fac=dfridr_CON2;
	for (j=2;j<=i;j++) {
	    a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	    fac=dfridr_CON2*fac;
	    errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
	    if (errt <= *err) {
		*err=errt;
		ans=a[j][i];
	    }
	}
	if (fabs(a[i][i]-a[i-1][i-1]) >= dfridr_SAFE*(*err)) break;
    }
    free_dmatrix(a,1,dfridr_NTAB,1,dfridr_NTAB);
    return ans;
}

//--- factrl ---
double NumericalRecipes::factrl(int n)
{
  static int ntop=4;
  static double a[33]={1.0,1.0,2.0,6.0,24.0};
  int j;

  if (n < 0) nrerror("Negative factorial in routine factrl");
  if (n > 32) return exp(gammln(n+1.0));
  while (ntop<n) {
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}

//--- fdjac ---
void NumericalRecipes::fdjac(int n, double x[], double fvec[], double **df,
			     void (*vecfunc)(int, double [], double []))
{
  int i,j;
  double h,temp,*f;

  f=dvector(1,n);
  for (j=1;j<=n;j++) {
    temp=x[j];
    h=fdjac_EPS*fabs(temp);
    if (h == 0.0) h=fdjac_EPS;
    x[j]=temp+h;
    h=x[j]-temp;
    (*vecfunc)(n,x,f);
    x[j]=temp;
    for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
  }
  free_dvector(f,1,n);
}

//--- fmin ---
double NumericalRecipes::fmin(double x[])
{
  int i;
  double sum;

  (*nrfuncv)(nn,x,fvec);
  for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
  return 0.5*sum;
}

//--- four1 ---
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void NumericalRecipes::four1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP

//--- fourn ---
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void NumericalRecipes::fourn(double data[], unsigned long nn[], int ndim, int isign)
{
  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  double tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;

  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	  for (i3=i1;i3<=ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    SWAP(data[i3],data[i3rev]);
	    SWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
	for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	  for (i2=i1;i2<=ip3;i2+=ifp2) {
	    k1=i2;
	    k2=k1+ifp1;
	    tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
	    tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
#undef SWAP

//--- gammln ---
double NumericalRecipes::gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

//--- lnsrch ---
void NumericalRecipes::lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
                              double *f, double stpmax, int *check, double (*func)(double []))
{
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0,fold2=0.0,rhs1,rhs2,slope,sum,temp,test,tmplam=0.0;

  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=std::sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=lnsrch_TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
      return;
    } else if (*f <= fold+lnsrch_ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) nrerror("Roundoff problem in lnsrch.");
	  else tmplam=(-b+std::sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}

//--- lubksb ---
void NumericalRecipes::lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

//--- ludcmp ---
void NumericalRecipes::ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv=dvector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=ludcmp_TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_dvector(vv,1,n);
}

//--- mmid ---
void NumericalRecipes::mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
	  double yout[], void (*derivs)(double, double[], double[]))
{
  int n,i;
  double x,swap,h2,h,*ym,*yn;

  ym=dvector(1,nvar);
  yn=dvector(1,nvar);
  h=htot/nstep;
  for (i=1;i<=nvar;i++) {
    ym[i]=y[i];
    yn[i]=y[i]+h*dydx[i];
  }
  x=xs+h;
  (*derivs)(x,yn,yout);
  h2=2.0*h;
  for (n=2;n<=nstep;n++) {
    for (i=1;i<=nvar;i++) {
      swap=ym[i]+h2*yout[i];
      ym[i]=yn[i];
      yn[i]=swap;
    }
    x += h;
    (*derivs)(x,yn,yout);
  }
  for (i=1;i<=nvar;i++)
    yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
  free_dvector(yn,1,nvar);
  free_dvector(ym,1,nvar);
}

//--- mnbrak ---
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void NumericalRecipes::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
			      double (*func)(double))
{
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+mnbrak_GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(std::fabs(q-r),mnbrak_TINY),q-r));
    ulim=(*bx)+mnbrak_GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
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
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+mnbrak_GOLD*(*cx-*bx))
	SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+mnbrak_GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}
#undef SHFT

//--- newt ---
int NumericalRecipes::nn = 0;
double * NumericalRecipes::fvec = 0;
void (*NumericalRecipes::nrfuncv)(int n, double v[], double f[]) = 0;

#define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);\
	free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}
void NumericalRecipes::newt(double x[], int n, int *check, void (*vecfunc)(int, double [], double []))
{
  int i,its,j,*indx;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

  indx=ivector(1,n);
  fjac=dmatrix(1,n,1,n);
  g=dvector(1,n);
  p=dvector(1,n);
  xold=dvector(1,n);
  fvec=dvector(1,n);
  nn=n;
  nrfuncv=vecfunc;
  f=fmin(x);
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test < 0.01*newt_TOLF) {
    *check=0;
    FREERETURN
      }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=newt_STPMX*FMAX(std::sqrt(sum),(double)n);
  for (its=1;its<=newt_MAXITS;its++) {
    fdjac(n,x,fvec,fjac,vecfunc);
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;
    for (i=1;i<=n;i++) p[i] = -fvec[i];
    ludcmp(fjac,n,indx,&d);
    lubksb(fjac,n,indx,p);
    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < newt_TOLF) {
      *check=0;
      FREERETURN
	}
    if (*check) {
      test=0.0;
      den=FMAX(f,0.5*n);
      for (i=1;i<=n;i++) {
	temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      *check=(test < newt_TOLMIN ? 1 : 0);
      FREERETURN
	}
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < newt_TOLX) FREERETURN
			    }
  nrerror("MAXITS exceeded in newt");
}
#undef FREERETURN

//--- odeint ---
int NumericalRecipes::kmax = 0;
int NumericalRecipes::kount = 0;
double * NumericalRecipes::xp = 0;
double ** NumericalRecipes::yp = 0;
double NumericalRecipes::dxsav = 0.0;

void NumericalRecipes::odeint(double ystart[], int nvar, double x1, double x2, double eps,
			      double *h1,double hmin, int *nok, int *nbad, 
			      void (*derivs)(double, double [], double []),
			      void (*rkqs)(double [], double [], int, double *, double, double, 
					   double [],double *, double *,
					   void (*)(double, double [], double []), double), 
			      double hmax)
{
  int nstp,i;
  double xsav=0.0,x,hnext,hdid,h;
  double *yscal,*y,*dydx;
  double odeint_MAXSTP_local=odeint_MAXSTP;
  
  if (hmax != 0.0 && fabs(*h1)>fabs(hmax)) 
    (*h1)=hmax; 
  if (hmax != 0.0 && odeint_MAXSTP < (int) fabs((x2-x1)/hmax))
    odeint_MAXSTP_local = 2.0*fabs((x2-x1)/hmax);

  yscal=dvector(1,nvar);
  y=dvector(1,nvar);
  dydx=dvector(1,nvar);
  x=x1;
  h=SIGN(*h1,x2-x1);
  *nok = (*nbad) = kount = 0;
  for (i=1;i<=nvar;i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=1;nstp<=odeint_MAXSTP_local;nstp++) {
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+odeint_TINY;
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp[++kount]=x;
      for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,hmax);
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++) ystart[i]=y[i];
      if (kmax) {
	xp[++kount]=x;
	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
      }
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      *h1 = h;
      return;
    }
    if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
    h=hnext;
  }
  nrerror("Too many steps in routine odeint");
}

//--- plgndr ---
double NumericalRecipes::plgndr(int l, int m, double x)
{
  double fact,pll=0.0,pmm,pmmp1,somx2;
  int i,ll;

  if (m < 0 || m > l || fabs(x) > 1.0)
    nrerror("Bad arguments in routine plgndr");
  pmm=1.0;
  if (m > 0) {
    somx2=std::sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
	pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}

//--- pythag ---
double NumericalRecipes::pythag(double a, double b)
{
  double absa,absb;
  absa=std::fabs(a);
  absb=std::fabs(b);
  if (absa > absb) return absa*std::sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*std::sqrt(1.0+SQR(absa/absb)));
}

//--- pzextr ---
void NumericalRecipes::pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
  int k1,j;
  double q,f2,f1,delta,*c;

  c=dvector(1,nv);
  x[iest]=xest;
  for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
  if (iest == 1) {
    for (j=1;j<=nv;j++) d[j][1]=yest[j];
  } else {
    for (j=1;j<=nv;j++) c[j]=yest[j];
    for (k1=1;k1<iest;k1++) {
      delta=1.0/(x[iest-k1]-xest);
      f1=xest*delta;
      f2=x[iest-k1]*delta;
      for (j=1;j<=nv;j++) {
	q=d[j][k1];
	d[j][k1]=dy[j];
	delta=c[j]-q;
	dy[j]=f1*delta;
	c[j]=f2*delta;
	yz[j] += dy[j];
      }
    }
    for (j=1;j<=nv;j++) d[j][iest]=dy[j];
  }
  free_dvector(c,1,nv);
}

//--- realft ---
void NumericalRecipes::realft(double data[], unsigned long n, int isign)
{
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;

  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}

//--- rlft3 ---
void NumericalRecipes::rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,
			     unsigned long nn3, int isign)
{
  unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
  double theta,wi,wpi,wpr,wr,wtemp;
  double c1,c2,h1r,h1i,h2r,h2i;

  if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != long(nn1*nn2*nn3))
    nrerror("rlft3: problem with dimensions or contiguity of data array\n");
  c1=0.5;
  c2 = -0.5*isign;
  theta=isign*(6.28318530717959/nn3);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  nn[1]=nn1;
  nn[2]=nn2;
  nn[3]=nn3 >> 1;
  if (isign == 1) {
    fourn(&data[1][1][1]-1,nn,3,isign);
    for (i1=1;i1<=nn1;i1++)
      for (i2=1,j2=0;i2<=nn2;i2++) {
	speq[i1][++j2]=data[i1][i2][1];
	speq[i1][++j2]=data[i1][i2][2];
      }
  }
  for (i1=1;i1<=nn1;i1++) {
    j1=(i1 != 1 ? nn1-i1+2 : 1);
    wr=1.0;
    wi=0.0;
    for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
      for (i2=1;i2<=nn2;i2++) {
	if (i3 == 1) {
	  j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
	  h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
	  h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
	  h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
	  h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
	  data[i1][i2][1]=h1r+h2r;
	  data[i1][i2][2]=h1i+h2i;
	  speq[j1][j2]=h1r-h2r;
	  speq[j1][j2+1]=h2i-h1i;
	} else {
	  j2=(i2 != 1 ? nn2-i2+2 : 1);
	  j3=nn3+3-(i3<<1);
	  h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
	  h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
	  h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
	  h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
	  data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
	  data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
	  data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
	  data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
	}
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  }
  if (isign == -1)
    fourn(&data[1][1][1]-1,nn,3,isign);
}

//--- rkck ---
void NumericalRecipes::rkck(double y[], double dydx[], int n, double x, double h, double yout[],
			    double yerr[], void (*derivs)(double, double [], double []))
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

  ak2=dvector(1,n);
  ak3=dvector(1,n);
  ak4=dvector(1,n);
  ak5=dvector(1,n);
  ak6=dvector(1,n);
  ytemp=dvector(1,n);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=1;i<=n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  free_dvector(ytemp,1,n);
  free_dvector(ak6,1,n);
  free_dvector(ak5,1,n);
  free_dvector(ak4,1,n);
  free_dvector(ak3,1,n);
  free_dvector(ak2,1,n);
}

//--- rkqs ---
void NumericalRecipes::rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
			    double yscal[], double *hdid, double *hnext, 
			    void (*derivs)(double, double [], double []), double hmax)
{
  int i;
  double errmax,h,htemp,xnew,*yerr,*ytemp;

  yerr=dvector(1,n);
  ytemp=dvector(1,n);
  h=htry;
  for (;;) {
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    errmax=0.0;
    for (i=1;i<=n;i++) 
      errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if 
      (errmax <= 1.0) break;
    htemp=rkqs_SAFETY*h*std::pow(errmax,rkqs_PSHRNK);
    h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
    xnew=(*x)+h;
    if 
      (xnew == *x) nrerror("stepsize underflow in rkqs");
  }
  if 
    (errmax > rkqs_ERRCON) *hnext=rkqs_SAFETY*h*std::pow(errmax,rkqs_PGROW);
  else 
    *hnext=5.0*h;
  if 
    (hmax != 0.0 && std::fabs(hmax/(*hnext)) < 1.0 ) (*hnext)=(*hnext)*std::fabs(hmax/(*hnext));
  *x += (*hdid=h);
  for (i=1;i<=n;i++) 
    y[i]=ytemp[i];
  free_dvector(ytemp,1,n);
  free_dvector(yerr,1,n);
}

//--- rtflsp ---
double NumericalRecipes::rtflsp(double (*func)(double), double x1, double x2, double xacc)
{
  int j;
  double fl,fh,xl,xh,swap,dx,del,f,rtf;

  fl=(*func)(x1);
  fh=(*func)(x2);
  if (fl*fh > 0.0) nrerror("Root must be bracketed in rtflsp");
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xl=x2;
    xh=x1;
    swap=fl;
    fl=fh;
    fh=swap;
  }
  dx=xh-xl;
  for (j=1;j<=rtflsp_MAXIT;j++) {
    rtf=xl+dx*fl/(fl-fh);
    f=(*func)(rtf);
    if (f < 0.0) {
      del=xl-rtf;
      xl=rtf;
      fl=f;
    } else {
      del=xh-rtf;
      xh=rtf;
      fh=f;
    }
    dx=xh-xl;
    if (fabs(del) < xacc || f == 0.0) return rtf;
  }
  nrerror("Maximum number of iterations exceeded in rtflsp");
  return 0.0;
}

//--- rtsafe ---
double NumericalRecipes::rtsafe(void (*funcd)(double, double *, double *),
				double x1, double x2, double xacc)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  (*funcd)(x1,&fl,&df);
  (*funcd)(x2,&fh,&df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    nrerror("Root must be bracketed in rtsafe");
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
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df);
  for (j=1;j<=rtsafe_MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
	|| (fabs(2.0*f) > fabs(dxold*df))) {
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
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,&f,&df);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  nrerror("Maximum number of iterations exceeded in rtsafe");
  return 0.0;
}

//--- savgol ---
void NumericalRecipes::savgol(double c[], int np, int nl, int nr, int ld, int m)
{
  int imj,ipj,j,k,kk,mm,*indx;
  double d,fac,sum,**a,*b;

  if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m)
    nrerror("bad args in savgol");
  indx=ivector(1,m+1);
  a=dmatrix(1,m+1,1,m+1);
  b=dvector(1,m+1);
  for (ipj=0;ipj<=(m << 1);ipj++) {
    sum=(ipj ? 0.0 : 1.0);
    for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
    for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
    mm=IMIN(ipj,2*m-ipj);
    for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
  }
  ludcmp(a,m+1,indx,&d);
  for (j=1;j<=m+1;j++) b[j]=0.0;
  b[ld+1]=1.0;
  lubksb(a,m+1,indx,b);
  for (kk=1;kk<=np;kk++) c[kk]=0.0;
  for (k = -nl;k<=nr;k++) {
    sum=b[1];
    fac=1.0;
    for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
    kk=((np-k) % np)+1;
    c[kk]=sum;
  }
  free_dvector(b,1,m+1);
  free_dmatrix(a,1,m+1,1,m+1);
  free_ivector(indx,1,m+1);
}

//--- shootf ---
int NumericalRecipes::nn2 = 0;
int NumericalRecipes::nvar = 0;
double NumericalRecipes::x1 = 0.0;
double NumericalRecipes::x2 = 0.0;
double NumericalRecipes::xf = 0.0;
void (*NumericalRecipes::derivs)(double x, double y[], double dydx[]) = 0;
void (*NumericalRecipes::load1)(double x1, double v1[], double y[]) = 0;
void (*NumericalRecipes::load2)(double x2, double v2[], double y[]) = 0;
void (*NumericalRecipes::score)(double xf, double y[], double f[]) = 0;

void NumericalRecipes::shootf(int n, double v[], double f[])
{
  int i,nbad,nok;
  double h1,hmin=0.0,*f1,*f2,*y;

  f1=dvector(1,nvar);
  f2=dvector(1,nvar);
  y=dvector(1,nvar);
  kmax=0;
  h1=(x2-x1)/1000.0;
  (*load1)(x1,v,y);
  odeint(y,nvar,x1,xf,shootf_EPS,&h1,hmin,&nok,&nbad,derivs,&rkqs);
  (*score)(xf,y,f1);
  h1=(x2-x1)/1000.0;
  (*load2)(x2,&v[nn2],y);
  odeint(y,nvar,x2,xf,shootf_EPS,&h1,hmin,&nok,&nbad,derivs,&rkqs);
  (*score)(xf,y,f2);
  for (i=1;i<=n;i++) f[i]=f1[i]-f2[i];
  free_dvector(y,1,nvar);
  free_dvector(f2,1,nvar);
  free_dvector(f1,1,nvar);
}

//--- sinft ---
void NumericalRecipes::sinft(double y[], int n)
{
  int j,n2=n+2;
  double sum,y1,y2;
  double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

  theta=3.14159265358979/(double) n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  y[1]=0.0;
  for (j=2;j<=(n>>1)+1;j++) {
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
    y1=wi*(y[j]+y[n2-j]);
    y2=0.5*(y[j]-y[n2-j]);
    y[j]=y1+y2;
    y[n2-j]=y1-y2;
  }
  realft(y,n,1);
  y[1]*=0.5;
  sum=y[2]=0.0;
  for (j=1;j<=n-1;j+=2) {
    sum += y[j];
    y[j]=y[j+1];
    y[j+1]=sum;
  }
}

//--- spharm ---
std::complex<double> NumericalRecipes::spharm(int l, int m, double theta, double phi)
{
  if ((m>0 && m>l) || (m<0 && m < -l)) return 0;
  return (m<0 ? std::pow(-1.0,m) : 1.0)*std::sqrt((2.0*l+1.0)*factrl(l-std::abs(m))/(4.0*M_PI*factrl(l+std::abs(m))))*
    plgndr(l,abs(m),std::cos(theta))*std::complex<double>(std::cos(m*phi),std::sin(m*phi));
}

//--- spharm_e ---
double NumericalRecipes::spharm_e(int l, int m, double theta, double phi)
{
  if ((m>0 && m>l) || (m<0)) return 0.0;
  return std::sqrt((2.0*l+1.0)*factrl(l-m)/(2.0*M_PI*factrl(l+m)))*
    plgndr(l,m,std::cos(theta))*std::cos(m*phi);
}

//--- spharm_o ---
double NumericalRecipes::spharm_o(int l, int m, double theta, double phi)
{
  if ((m>0 && m>l) || (m<0)) return 0.0;
  return std::sqrt((2.0*l+1.0)*factrl(l-m)/(2.0*M_PI*factrl(l+m)))*
    plgndr(l,m,std::cos(theta))*std::sin(m*phi);
}

//--- svbksb ---
void NumericalRecipes::svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]){
  int jj,j,i;
  double s,*tmp;
 
  tmp=dvector(1,n);
  for (j=1;j<=n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=1;i<=m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=1;j<=n;j++) {
    s=0.0;
    for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free_dvector(tmp,1,n);
}

//--- svdcmp ---
void NumericalRecipes::svdcmp(double **a, int m, int n, double w[], double **v)
{
  int flag,i,its,j,jj,k,l=0,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1=dvector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += std::fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(std::sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += std::fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(std::sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(std::fabs(w[i])+std::fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((double)(std::fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(std::fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(std::fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_dvector(rv1,1,n);
}

//--- svdfit ---
void NumericalRecipes::svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma,
			      double **u, double **v, double w[], double *chisq,
			      void (*funcs)(double, double [], int))
{
  int j,i;
  double wmax,tmp,thresh,sum,*b,*afunc;

  b=dvector(1,ndata);
  afunc=dvector(1,ma);
  for (i=1;i<=ndata;i++) {
    (*funcs)(x[i],afunc,ma);
    tmp=1.0/sig[i];
    for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=svdfit_TOL*wmax;
  for (j=1;j<=ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  svbksb(u,w,v,ndata,ma,b,a);
  *chisq=0.0;
  for (i=1;i<=ndata;i++) {
    (*funcs)(x[i],afunc,ma);
    for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
    *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  free_dvector(afunc,1,ma);
  free_dvector(b,1,ndata);
}

//--- tridag ---
void NumericalRecipes::tridag(double a[], double b[], double c[], double r[], double u[],
			      unsigned long n)
{
  unsigned long j;
  double bet,*gam;
 
  gam = dvector(1,n);
  if (b[1] == 0.0)
    nrerror("Error 1 in tridag");
  u[1] = r[1]/(bet=b[1]);
  for (j=2; j<=n; ++j)
  {
    gam[j] = c[j-1]/bet;
    bet = b[j] - a[j]*gam[j];
    if (bet == 0.0)
      nrerror("Error 2 in tridag");
    u[j] = (r[j] - a[j]*u[j-1])/bet;
  }
  for (j=(n-1); j>=1; --j)
    u[j] -= gam[j+1]*u[j+1];
  free_dvector(gam,1,n);
}


//--- twofft ---
void NumericalRecipes::twofft(double data1[], double data2[], double fft1[], double fft2[],
			      unsigned long n)
{
  unsigned long nn3,nn2,jj,j;
  double rep,rem,aip,aim;

  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2)
  {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  four1(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2)
  {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}

//--- zbrac ---
int NumericalRecipes::zbrac(double (*func)(double), double *x1, double *x2)
{
  int j;
  double f1,f2;

  if (*x1 == *x2)
    nrerror("Bad initial range in zbrac");
  f1=(*func)(*x1);
  f2=(*func)(*x2);
  for (j=1;j<=zbrac_NTRY;j++)
  {
    if (f1*f2 < 0.0)
      return 1;
    if (std::fabs(f1) < std::fabs(f2))
      f1=(*func)(*x1 += zbrac_FACTOR*(*x1-*x2));
    else
      f2=(*func)(*x2 += zbrac_FACTOR*(*x2-*x1));
  }
  return 0;
}

//--- zbrent ---
double NumericalRecipes::zbrent(double (*func)(double), double x1, double x2, double tol)
{
  int iter;
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    nrerror("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=zbrent_ITMAX;iter++)
  {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
    {
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
    fb=(*func)(b);
  }
  nrerror("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
