#include<config.h>
#include<math.h>

#ifndef GMP
#include<alg/alg_remez.h>
#else

CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// alg_remez.C
//
// AlgRemez is relevant to the Rational Hybrid Molecular Dynamics
// algorithm.  The Remez algorithm is used for generating the nth root
// rational approximation.
//
// Note this class can only be used if
// the gnu multiprecision (GNU MP) library is present.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <string.h>
#include <alg/common_arg.h>
#include <alg/cg_arg.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_remez.h>
CPS_START_NAMESPACE

// Constructor
AlgRemez::AlgRemez(Float lower, Float upper, long precision) 
{
  cname = "AlgRemez";
  char *fname = "AlgRemez()";
  VRB.Func(cname,fname);

  prec = precision;
  bigfloat::setDefaultPrecision(prec);

  apstrt = lower;
  apend = upper;
  apwidt = apend - apstrt;
  VRB.Flow(cname,fname,"Approximation bounds are [%e,%e]\n", (Float)apstrt,(Float)apend);

  alloc = 0;
  n = 0;
  d = 0;

  // Only require the approximation spread to be less than tolerance 
  if (sizeof(Float)==sizeof(double)) {
    tolerance = 1e-15;
  } else if (sizeof(Float)==sizeof(float)) {
    tolerance = 1e-8;
  }

}

// Destructor
AlgRemez::~AlgRemez()
{
  char *fname = "~AlgRemez()";
  VRB.Func(cname,fname);

  if (alloc) {
    VRB.Sfree(cname,fname, "param",param);
    delete [] param;
    VRB.Sfree(cname,fname, "roots",roots);
    delete [] roots;
    VRB.Sfree(cname,fname, "poles",poles);
    delete [] poles;
    VRB.Sfree(cname,fname, "xx",xx);
    delete [] xx;
    VRB.Sfree(cname,fname, "mm",mm);
    delete [] mm;
  }
}

// Free memory and reallocate as necessary
void AlgRemez::allocate(int degree)
{
  char *fname = "allocate(int)";
  VRB.Func(cname,fname);

  // Arrays have previously been allocated, deallocate first, then allocate
  if (alloc) {
    VRB.Sfree(cname,fname, "param",param);
    delete [] param;
    VRB.Sfree(cname,fname, "roots",roots);
    delete [] roots;
    VRB.Sfree(cname,fname, "poles",poles);
    delete [] poles;
    VRB.Sfree(cname,fname, "xx",xx);
    delete [] xx;
    VRB.Sfree(cname,fname, "mm",mm);
    delete [] mm;
  }

  int size = 2*degree + 1;

  // Note use of new and delete in memory allocation - cannot run on qcdsp
  param = new bigfloat[size+1];
  if(param == 0) ERR.Pointer(cname,fname,"param");
  VRB.Smalloc(cname,fname,"param",param,size+1 * sizeof(bigfloat));

  roots = new bigfloat[size];
  if(roots == 0) ERR.Pointer(cname,fname,"roots");
  VRB.Smalloc(cname,fname,"roots",roots,degree * sizeof(bigfloat));

  poles = new bigfloat[size];
  if(poles == 0) ERR.Pointer(cname,fname,"poles");
  VRB.Smalloc(cname,fname,"poles",poles,degree * sizeof(bigfloat));

  xx = new bigfloat[size+2];
  if(xx== 0) ERR.Pointer(cname,fname,"xx");
  VRB.Smalloc(cname,fname,"xx",xx,size+2 * sizeof(bigfloat));

  mm = new bigfloat[size+2];
  if(mm == 0) ERR.Pointer(cname,fname,"mm");
  VRB.Smalloc(cname,fname,"mm",mm,size+2 * sizeof(bigfloat));

  alloc = 1;
}

// Reset the bounds of the approximation
void AlgRemez::setBounds(Float lower, Float upper)
{
  char *fname = "setbounds(Float, Float)";
  VRB.Func(cname,fname);

  apstrt = lower;
  apend = upper;
  apwidt = apend - apstrt;
}

// Generate the rational approximation x^(pnum/pden)
Float AlgRemez::generateApprox(int degree, unsigned long pnum, unsigned long pden)
{
  char *fname = "generateApprox(int, unsigned long, unsigned long)";
  VRB.Func(cname,fname);

  // Reallocate arrays, since degree has changed
  if (degree != n) allocate(degree);

  step = new bigfloat[2*degree+2];
  if(step == 0) ERR.Pointer(cname,fname,"step");
  VRB.Smalloc(cname,fname,"step",step,(2*degree+2) * sizeof(bigfloat));

  power_num = pnum;
  power_den = pden;
  spread = 1.0e37;
  iter = 0;

  d = degree;
  n = degree;
  neq = n + d + 1;

  initialGuess();
  stpini(step);

  while (spread > tolerance) { //iterate until convergance

    if (iter++%100==0) 
      VRB.Flow(cname,fname,"Iteration %d, spread %e delta %e\n", iter-1,(Float)spread,(Float)delta);
    equations();

    if (delta < tolerance)
      ERR.General(cname, fname,"Delta too small, try increasing precision\n");

    search(step);

  }

  int sign;
  Float error = (Float)getErr(mm[0],&sign);
  VRB.Result(cname,fname,"Converged at %d iterations, error = %e\n",iter,error);

  // Once the approximation has been generated, calculate the roots
  if(!root()) ERR.General(cname,fname,"Root finding failed\n");
  
  VRB.Sfree(cname,fname, "step",step);
  delete [] step;

  // Return the maximum error in the approximation
  return error;
}

// Return the partial fraction expansion of the approximation x^(pnum/pden)
int AlgRemez::getPFE(Float *Res, Float *Pole, Float *Norm) {
  char *fname = "getPFE(Float*, Float*, Float*)";
  VRB.Func(cname,fname);

  if (!alloc) ERR.General(cname,fname,"Approximation not yet generated\n");

  bigfloat *r = new bigfloat[n];
  if(r == 0) ERR.Pointer(cname,fname,"r");
  VRB.Smalloc(cname,fname,"r",r,n * sizeof(bigfloat));
  bigfloat *p = new bigfloat[d];
  if(p == 0) ERR.Pointer(cname,fname,"p");
  VRB.Smalloc(cname,fname,"p",p,d * sizeof(bigfloat));
  
  for (int i=0; i<n; i++) {
    r[i] = roots[i];
    p[i] = poles[i];
  }
  
  // Perform a partial fraction expansion
  pfe(r, p, norm);

  // Convert to Float and return
  *Norm = (Float)norm;
  for (int i=0; i<n; i++) {
    Res[i] = (Float)r[i];
    Pole[i] = (Float)p[i];
  }

  VRB.Sfree(cname,fname, "r",r);
  delete [] r;
  VRB.Sfree(cname,fname, "p",p);
  delete [] p;

  // Where the smallest shift is located
  return 0;
}

// Return the partial fraction expansion of the approximation x^(-pnum/pden)
int AlgRemez::getIPFE(Float *Res, Float *Pole, Float *Norm) {
  char *fname = "getIPFE(Float*, Float*, Float*)";
  VRB.Func(cname,fname);

  if (!alloc) ERR.General(cname,fname,"Approximation not yet generated\n");

  bigfloat *r = new bigfloat[d];
  if(r == 0) ERR.Pointer(cname,fname,"r");
  VRB.Smalloc(cname,fname,"r",r,n * sizeof(bigfloat));
  bigfloat *p = new bigfloat[n];
  if(p == 0) ERR.Pointer(cname,fname,"p");
  VRB.Smalloc(cname,fname,"p",p,d * sizeof(bigfloat));
  
  // Want the inverse function
  for (int i=0; i<n; i++) {
    r[i] = poles[i];
    p[i] = roots[i];
  }

  // Perform a partial fraction expansion
  pfe(r, p, (bigfloat)1l/norm);

  // Convert to Float and return
  *Norm = (Float)((bigfloat)1l/(norm));
  for (int i=0; i<n; i++) {
    Res[i] = (Float)r[i];
    Pole[i] = (Float)p[i];
  }

  VRB.Sfree(cname,fname, "r",r);
  delete [] r;
  VRB.Sfree(cname,fname, "p",p);
  delete [] p;

  // Where the smallest shift is located
  return 0;
}

// Initial values of maximal and minimal errors
void AlgRemez::initialGuess() {
  char *fname = "initialGuess()";
  VRB.Func(cname,fname);
  // Supply initial guesses for solution points
  long ncheb = neq;			// Degree of Chebyshev error estimate
  bigfloat a, r;

  // Find ncheb+1 extrema of Chebyshev polynomial

  a = ncheb;
  mm[0] = apstrt;
  for (long i = 1; i < ncheb; i++) {
    r = 0.5 * (1 - cos((M_PI * i)/(double) a));
    //r *= sqrt_bf(r);
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    mm[i] = apstrt + r * apwidt;
  }
  mm[ncheb] = apend;

  a = 2.0 * ncheb;
  for (long i = 0; i <= ncheb; i++) {
    r = 0.5 * (1 - cos(M_PI * (2*i+1)/(double) a));
    //r *= sqrt_bf(r); // Squeeze to low end of interval
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    xx[i] = apstrt + r * apwidt;
  }
}

// Initialise step sizes
void AlgRemez::stpini(bigfloat *step) {
  char *fname = "stpini()";
  VRB.Func(cname,fname);

  xx[neq+1] = apend;
  delta = 0.25;
  step[0] = xx[0] - apstrt;
  for (int i = 1; i < neq; i++) step[i] = xx[i] - xx[i-1];
  step[neq] = step[neq-1];
}

// Search for error maxima and minima
void AlgRemez::search(bigfloat *step) {
  char *fname = "search()";
  VRB.Func(cname,fname);
  bigfloat a, q, xm, ym, xn, yn, xx0, xx1;
  int i, j, meq, emsign, ensign, steps;

  meq = neq + 1;
  bigfloat *yy = new bigfloat[meq];
  if(yy == 0) ERR.Pointer(cname,fname,"yy");
  VRB.Smalloc(cname,fname,"yy",yy,meq * sizeof(bigfloat));

  bigfloat eclose = 1.0e30;
  bigfloat farther = 0l;

  j = 1;
  xx0 = apstrt;

  for (i = 0; i < meq; i++) {
    steps = 0;
    xx1 = xx[i]; // Next zero
    if (i == meq-1) xx1 = apend;
    xm = mm[i];
    ym = getErr(xm,&emsign);
    q = step[i];
    xn = xm + q;
    if (xn < xx0 || xn >= xx1) {	// Cannot skip over adjacent boundaries
      q = -q;
      xn = xm;
      yn = ym;
      ensign = emsign;
    } else {
      yn = getErr(xn,&ensign);
      if (yn < ym) {
	q = -q;
	xn = xm;
	yn = ym;
	ensign = emsign;
      }
    }
  
    while(yn >= ym) {		// March until error becomes smaller.
      if (++steps > 10) break;
      ym = yn;
      xm = xn;
      emsign = ensign;
      a = xm + q;
      if (a == xm || a <= xx0 || a >= xx1) break;// Must not skip over the zeros either side.
      xn = a;
      yn = getErr(xn,&ensign);
    }

    mm[i] = xm;			// Position of maximum
    yy[i] = ym;			// Value of maximum

    if (eclose > ym) eclose = ym;
    if (farther < ym) farther = ym;

    xx0 = xx1; // Walk to next zero.
  } // end of search loop

  q = (farther - eclose);	// Decrease step size if error spread increased
  if (eclose != 0.0) q /= eclose; // Relative error spread
  if (q >= spread) delta *= 0.5; // Spread is increasing; decrease step size
  spread = q;

  for (i = 0; i < neq; i++) {
    q = yy[i+1];
    if (q != 0.0) q = yy[i] / q  - (bigfloat)1l;
    else q = 0.0625;
    if (q > (bigfloat)0.25) q = 0.25;
    q *= mm[i+1] - mm[i];
    step[i] = q * delta;
  }
  step[neq] = step[neq-1];
  
  for (i = 0; i < neq; i++) {	// Insert new locations for the zeros.
    xm = xx[i] - step[i];
    if (xm <= apstrt) continue;
    if (xm >= apend) continue;
    if (xm <= mm[i]) xm = (bigfloat)0.5 * (mm[i] + xx[i]);
    if (xm >= mm[i+1]) xm = (bigfloat)0.5 * (mm[i+1] + xx[i]);
    xx[i] = xm;
  }

  VRB.Sfree(cname,fname, "yy",yy);
  delete [] yy;
}

// Solve the equations
void AlgRemez::equations(void) {
  char *fname = "equations()";
  VRB.Func(cname,fname);
  bigfloat x, y, z;
  int i, j, ip;
  bigfloat *aa;

  bigfloat *AA = new bigfloat[(neq)*(neq)];
  if(AA == 0) ERR.Pointer(cname,fname,"AA");
  VRB.Smalloc(cname,fname,"AA",AA,(neq)*(neq)* sizeof(bigfloat));
  bigfloat *BB = new bigfloat[neq];
  if(BB == 0) ERR.Pointer(cname,fname,"BB");
  VRB.Smalloc(cname,fname,"BB",BB,(neq) * sizeof(bigfloat));
  
  for (i = 0; i < neq; i++) {	// set up the equations for solution by simq()
    ip = neq * i;		// offset to 1st element of this row of matrix
    x = xx[i];			// the guess for this row
    y = func(x);		// right-hand-side vector

    z = (bigfloat)1l;
    aa = AA+ip;
    for (j = 0; j <= n; j++) {
      *aa++ = z;
      z *= x;
    }

    z = (bigfloat)1l;
    for (j = 0; j < d; j++) {
      *aa++ = -y * z;
      z *= x;
    }
    BB[i] = y * z;		// Right hand side vector
  }

  // Solve the simultaneous linear equations.
  if (simq(AA, BB, param, neq)) ERR.General(cname,fname,"simq failed\n");;

  VRB.Sfree(cname,fname, "AA",AA);
  delete [] AA;
  VRB.Sfree(cname,fname, "BB",BB);
  delete [] BB;

}

// Evaluate the rational form P(x)/Q(x) using coefficients
// from the solution vector param
bigfloat AlgRemez::approx(const bigfloat x) {
  char *fname = "approx(bigfloat)";
  VRB.Func(cname,fname);
  bigfloat yn, yd;
  int i;

  // Work backwards toward the constant term.
  yn = param[n];		// Highest order numerator coefficient
  for (i = n-1; i >= 0; i--) yn = x * yn  +  param[i]; 
  yd = x + param[n+d];	// Highest degree coefficient = 1.0
  for (i = n+d-1; i > n; i--) yd = x * yd  +  param[i];

  return(yn/yd);
}

// Compute size and sign of the approximation error at x
bigfloat AlgRemez::getErr(bigfloat x, int *sign) {
  char *fname = "getErr(bigfloat)";
  VRB.Func(cname,fname);
  bigfloat e, f;

  f = func(x);
  e = approx(x) - f;
  if (f != 0) e /= f;
  if (e < (bigfloat)0.0) {
    *sign = -1;
    e = -e;
  }
  else *sign = 1;
  
  return(e);
}

// Calculate function required for the approximation
bigfloat AlgRemez::func(const bigfloat x) {
  char *fname = "func(bigfloat)";
  VRB.Func(cname,fname);

  bigfloat y,dy,f=1l,df;

  // initial guess to accelerate convergance
  y = (bigfloat)pow((double)x,(double)((bigfloat)1l/(bigfloat)power_den));
  while (abs_bf(f)>(bigfloat)1l/pow_bf((bigfloat)10,prec)) { // approx good to 10^(-prec)
    f = pow_bf(y,power_den) - x;
    df = (bigfloat)power_den*pow_bf(y,power_den-1);// need power_den-1 because of diff
    dy = f/df;
    y -= dy;
  }
  return pow_bf(y,power_num);
}

// Solve the system AX=B
int AlgRemez::simq(bigfloat A[], bigfloat B[], bigfloat X[], int n) {

  char *fname = "simq(bigfloat*, bigfloat*, bigfloat*, int, int)";
  VRB.Func(cname,fname);

  int i, j, ij, ip, ipj, ipk, ipn;
  int idxpiv, iback;
  int k, kp, kp1, kpk, kpn;
  int nip, nkp, nm1;
  bigfloat em, q, rownrm, big, size, pivot, sum;
  bigfloat *aa;

  int *IPS = (int*)smalloc((neq) * sizeof(int));		// simq() work vector
  if(IPS == 0) ERR.Pointer(cname,fname,"IPS");
  VRB.Smalloc(cname,fname,"IPS",IPS,(neq) * sizeof(int));

  nm1 = n - 1;
  // Initialize IPS and X
  
  ij = 0;
  for (i = 0; i < n; i++) {
    IPS[i] = i;
    rownrm = 0.0;
    for(j = 0; j < n; j++) {
      q = abs_bf(A[ij]);
      if(rownrm < q) rownrm = q;
      ++ij;
    }
    if (rownrm == (bigfloat)0l) {
      VRB.Result(cname,fname,"simq rownrm=0\n");
      VRB.Sfree(cname,fname, "IPS",IPS);
      delete [] IPS;
      return(1);
    }
    X[i] = (bigfloat)1.0 / rownrm;
  }
  
  for (k = 0; k < nm1; k++) {
    big = 0.0;
    idxpiv = 0;
    
    for (i = k; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      size = abs_bf(A[ipk]) * X[ip];
      if (size > big) {
	big = size;
	idxpiv = i;
      }
    }
    
    if (big == (bigfloat)0l) {
      VRB.Result(cname,fname,"simq big=0\n");
      VRB.Sfree(cname,fname, "IPS",IPS);
      delete [] IPS;
      return(2);
    }
    if (idxpiv != k) {
      j = IPS[k];
      IPS[k] = IPS[idxpiv];
      IPS[idxpiv] = j;
    }
    kp = IPS[k];
    kpk = n*kp + k;
    pivot = A[kpk];
    kp1 = k+1;
    for (i = kp1; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      em = -A[ipk] / pivot;
      A[ipk] = -em;
      nip = n*ip;
      nkp = n*kp;
      aa = A+nkp+kp1;
      for (j = kp1; j < n; j++) {
	ipj = nip + j;
	A[ipj] = A[ipj] + em * *aa++;
      }
    }
  }
  kpn = n * IPS[n-1] + n - 1;	// last element of IPS[n] th row
  if (A[kpn] == (bigfloat)0l) {
    VRB.Result(cname,fname,"simq A[kpn]=0\n");
    VRB.Sfree(cname,fname, "IPS",IPS);
    delete [] IPS;
    return(3);
  }

  
  ip = IPS[0];
  X[0] = B[ip];
  for (i = 1; i < n; i++) {
    ip = IPS[i];
    ipj = n * ip;
    sum = 0.0;
    for (j = 0; j < i; j++) {
      sum += A[ipj] * X[j];
      ++ipj;
    }
    X[i] = B[ip] - sum;
  }
  
  ipn = n * IPS[n-1] + n - 1;
  X[n-1] = X[n-1] / A[ipn];
  
  for (iback = 1; iback < n; iback++) {
    //i goes (n-1),...,1
    i = nm1 - iback;
    ip = IPS[i];
    nip = n*ip;
    sum = 0.0;
    aa = A+nip+i+1;
    for (j= i + 1; j < n; j++) 
      sum += *aa++ * X[j];
    X[i] = (X[i] - sum) / A[nip+i];
  }
  
  VRB.Sfree(cname,fname, "IPS",IPS);
  delete [] IPS;
  return(0);
}

// Calculate the roots of the approximation
int AlgRemez::root() {
  char *fname = "root()";
  VRB.Func(cname,fname);
  long i,j;
  bigfloat x,dx=0.05;
  bigfloat upper=1, lower=-10000;
  bigfloat tol = 1e-20;

  bigfloat *poly = new bigfloat[neq+1];
  if(poly == 0) ERR.Pointer(cname,fname,"poly");
  VRB.Smalloc(cname,fname,"poly",poly,(2*n+2) * sizeof(bigfloat));

  // First find the numerator roots
  for (i=0; i<=n; i++) poly[i] = param[i];
  for (i=n-1; i>=0; i--) {
    roots[i] = rtnewt(poly,i+1,lower,upper,tol);
    if (roots[i] == 0.0) {
      VRB.Warn(cname,fname,"Failure to converge on root\n");
      return 0;
    }
    poly[0] = -poly[0]/roots[i];
    for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/roots[i];
  }

  // Now find the denominator roots
  poly[d] = 1l;
  for (i=0; i<d; i++) poly[i] = param[n+1+i];
  for (i=d-1; i>=0; i--) {
    poles[i] = rtnewt(poly,i+1,lower,upper,tol);
    if (poles[i] == 0.0) {
      VRB.Warn(cname,fname,"Failure to converge on root");
      return 0;
    }
    poly[0] = -poly[0]/poles[i];
    for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/poles[i];
  }

  norm = param[n];
  VRB.Result(cname,fname,"Normalisation constant is %e\n",(double)norm);
  for (i=0; i<n; i++)
    VRB.Result(cname,fname,"%ld root = %e pole = %e\n",i,(Float)roots[i],(Float)poles[i]);

  VRB.Sfree(cname,fname, "poly",poly);
  delete [] poly;

  return 1;
}

// Evaluate the polynomial
bigfloat AlgRemez::polyEval(bigfloat x, bigfloat *poly, long size) {
  char *fname = "polyEval(bigfloat, bigfloat *, long)";
  VRB.Func(cname,fname);
  bigfloat f = poly[size];
  for (int i=size-1; i>=0; i--) f = f*x + poly[i];
  return f;
}

// Evaluate the differential of the polynomial
bigfloat AlgRemez::polyDiff(bigfloat x, bigfloat *poly, long size) {
  char *fname = "polyDiff(bigfloat, bigfloat *, long)";
  VRB.Func(cname,fname);
  bigfloat df = (bigfloat)size*poly[size];
  for (int i=size-1; i>0; i--) df = df*x + (bigfloat)i*poly[i];
  return df;
}

// Newton's method to calculate roots
bigfloat AlgRemez::rtnewt(bigfloat *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc) {
  char *fname = "rtnnewt(bigfloat *, long, bigfloatm bigfloat, bigfloat)";
  VRB.Func(cname,fname);
  int j;
  bigfloat df, dx, f, rtn;
  
  rtn=(bigfloat)0.5*(x1+x2);
  for (j=1; j<=JMAX;j++) {
    f = polyEval(rtn, poly, i);
    df = polyDiff(rtn, poly, i);
    dx = f/df;
    rtn -= dx;
    if ((x1-rtn)*(rtn-x2) < (bigfloat)0.0)
      VRB.Warn(cname,fname,"Jumped out of brackets in rtnewt\n");
    if (abs_bf(dx) < xacc) return rtn;
  }
  VRB.Warn(cname,fname,"Maximum number of iterations exceeded in rtnewt\n");
  return 0.0;
}

// Evaluate the partial fraction expansion of the rational function
// with res roots and poles poles.  Result is overwritten on input
// arrays.
void AlgRemez::pfe(bigfloat *res, bigfloat *poles, bigfloat norm) {
  char *fname = "pfe(bigfloat *res, bigfloat *poles)";
  VRB.Func(cname,fname);

  int i,j,small;
  bigfloat temp;
  bigfloat *numerator = new bigfloat[n];
  bigfloat *denominator = new bigfloat[d];

  // Construct the polynomials explicitly 
  for (i=1; i<n; i++) {
    numerator[i] = 0l;
    denominator[i] = 0l;
  }
  numerator[0]=1l;
  denominator[0]=1l;

  for (j=0; j<n; j++) {
    for (i=n-1; i>=0; i--) {
      numerator[i] *= -res[j];
      denominator[i] *= -poles[j];
      if (i>0) {
	numerator[i] += numerator[i-1];
	denominator[i] += denominator[i-1];
      }
    }
  }

  // Convert to proper fraction form.
  // Fraction is now in the form 1 + n/d, where O(n)+1=O(d)
  for (i=0; i<n; i++) numerator[i] -= denominator[i];

  // Find the residues of the partial fraction expansion and absorb the
  // coefficients.
  for (i=0; i<n; i++) {
    res[i] = 0l;
    for (j=n-1; j>=0; j--) {
      res[i] = poles[i]*res[i]+numerator[j];
    }
    for (j=n-1; j>=0; j--) {
      if (i!=j) res[i] /= poles[i]-poles[j];
    }
    res[i] *= norm;
  }  

  // res now holds the residues
  j = 0;
  for (i=0; i<n; i++) poles[i] = -poles[i];

  // Move the ordering of the poles from smallest to largest
  for (j=0; j<n; j++) {
    small = j;
    for (i=j+1; i<n; i++) {
      if (poles[i] < poles[small]) small = i;
    }
    if (small != j) {
      temp = poles[small];
      poles[small] = poles[j];
      poles[j] = temp;
      temp = res[small];
      res[small] = res[j];
      res[j] = temp;
    }
    VRB.Flow(cname,fname,"Residue = %e, Pole = %e\n", (double)res[j], (double)poles[j]);
  }

  delete [] numerator;
  delete [] denominator;
}

CPS_END_NAMESPACE

#endif