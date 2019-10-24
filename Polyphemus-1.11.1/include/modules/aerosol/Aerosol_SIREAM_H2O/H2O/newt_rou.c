#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"

#define NRANSI

#define ALF 1.0e-4 /* 1.0e-10 */
#define EPS 1.0e-4 /* 1.0e-4 */
#define MAXITS 100 /*400*/
#define STPMX 100.0
#define TINY 1.0e-20;
#define TOLF 1.0e-3 /* 1.0e-4 */
#define TOLMIN 1.0e-6   /* 1.0e-9 */
#define TOLX 1.0e-3    /*1.0e-7*/

#define NR_END 1
#define FREE_ARG char*


#define FREEDRETURN {free_dvector(fvec_double,1,n);			\
    free_dvector(xold_double,1,n); free_dvector(p_double,1,n);		\
    free_dvector(g_double,1,n);free_dmatrix(fjac_double,1,n,1,n);	\
    free_ivector(indx,1,n);return;}
#define FREEDRETURN1 {free_dvector1(fvec1_double,1,n);			\
    free_dvector1(xold_double,1,n); free_dvector1(p_double,1,n);        \
    free_dvector1(g_double,1,n);free_dmatrix1(fjac_double,1,n,1,n);	\
    free_ivector1(indx,1,n);return;}
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);        \
    free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);	\
    free_ivector(indx,1,n);return;}
#define FREERETURN1 {free_vector1(fvec1,1,n);free_vector1(xold,1,n);	\
    free_vector1(p,1,n);free_vector1(g,1,n);free_matrix1(fjac,1,n,1,n);	\
    free_ivector1(indx,1,n);return;}

int nn, nn1;
float *fvec, *fvec1;
void (*nrfuncv)(int n, float v[], float f[]);
void (*nrfuncv1)(int n, float v[], float f[]);

double *fvec_double, *fvec1_double;
void (*nrfuncv_double)(int n, double v[], double f[]);
void (*nrfuncv1_double)(int n, double v[], double f[]);

/* global variables */
extern int anrerrflag, bnrerrflag;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(nn, nn1, fvec, fvec1, nrfuncv, nrfuncv1)
#pragma omp threadprivate(fvec_double, fvec1_double, nrfuncv_double,	\
                          nrfuncv1_double)
#pragma omp threadprivate(anrerrflag, bnrerrflag)
#endif

/*****************NEWT FUNCTIONS THAT USE FLOAT VARS *******************/


void fdjac(int n, float x[], float fvec[], float **df,
           void (*vecfunc)(int, float [], float []))
{
  int i, j;
  float h, temp, *f;

  f = vector(1, n);
  for (j = 1; j <= n; j++)
    {
      temp = x[j];
      h = EPS * fabs(temp);
      if (h == 0.0) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      (*vecfunc)(n, x, f);
      x[j] = temp;
      for (i = 1; i <= n; i++) df[i][j] = (f[i] - fvec[i]) / h;
    }
  free_vector(f, 1, n);
}

void fdjac1(int n, float x[], float fvec[], float **df,
            void (*vecfunc)(int, float [], float []))
{
  int i, j;
  float h, temp, *f;

  f = vector1(1, n);
  for (j = 1; j <= n; j++)
    {
      temp = x[j];
      h = EPS * fabs(temp);
      if (h == 0.0) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      (*vecfunc)(n, x, f);
      x[j] = temp;
      for (i = 1; i <= n; i++) df[i][j] = (f[i] - fvec[i]) / h;
    }
  free_vector1(f, 1, n);
}

float fmin0(float x[])
{
  int i;
  float sum;
  (*nrfuncv)(nn, x, fvec);
  for (sum = 0.0, i = 1; i <= nn; i++) sum += SQR(fvec[i]);
  return 0.5 * sum;
}

float fmin1(float x[])
{
  int i;
  float sum;

  (*nrfuncv1)(nn1, x, fvec1);
  for (sum = 0.0, i = 1; i <= nn1; i++) sum += SQR(fvec1[i]);
  return 0.5 * sum;
}

void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
            float *f, float stpmax, int *check, float(*func)(float []))
{
  int i;
  float a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp,
    test, tmplam;

  /* These initializations have been added to remove warning messages
     from compilers. They don't change anything to the method as
     f2 and alam2 are initialized later in the loop before being used. */
  f2 = 0.;
  alam2 = 0.;

  *check = 0;
  for (sum = 0.0, i = 1; i <= n; i++) sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++) p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.", 'b');
  test = 0.0;
  for (i = 1; i <= n; i++)
    {
      temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
      if (temp > test) test = temp;
    }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;)
    {
      for (i = 1; i <= n; i++) x[i] = xold[i] + alam * p[i];
      *f = (*func)(x);
      if (alam < alamin)
        {
          for (i = 1; i <= n; i++) x[i] = xold[i];
          *check = 1;
          return;
        }
      else if (*f <= fold + ALF * alam * slope) return;
      else
        {
          if (alam == 1.0)
            tmplam = -slope / (2.0 * (*f - fold - slope));
          else
            {
              rhs1 = *f - fold - alam * slope;
              rhs2 = f2 - fold - alam2 * slope;
              a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
              b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
              if (a == 0.0) tmplam = -slope / (2.0 * b);
              else
                {
                  disc = b * b - 3.0 * a * slope;
                  if (disc < 0.0) tmplam = 0.5 * alam;
                  else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0 * a);
                  else tmplam = -slope / (b + sqrt(disc));
                }
              if (tmplam > 0.5 * alam)
                tmplam = 0.5 * alam;
            }
        }
      alam2 = alam;
      f2 = *f;
      alam = FMAX(tmplam, 0.1 * alam);
    }
}

void lnsrch1(int n, float xold[], float fold, float g[], float p[], float x[],
             float *f, float stpmax, int *check, float(*func)(float []))
{
  int i;
  float a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp,
    test, tmplam;

  /* These initializations have been added to remove warning messages
     from compilers. They don't change anything to the method as
     f2 and alam2 are initialized later in the loop before being used. */
  f2 = 0.;
  alam2 = 0.;

  *check = 0;
  for (sum = 0.0, i = 1; i <= n; i++) sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++) p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch1.", 'b');
  test = 0.0;
  for (i = 1; i <= n; i++)
    {
      temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
      if (temp > test) test = temp;
    }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;)
    {
      for (i = 1; i <= n; i++) x[i] = xold[i] + alam * p[i];
      *f = (*func)(x);
      if (alam < alamin)
        {
          for (i = 1; i <= n; i++) x[i] = xold[i];
          *check = 1;
          return;
        }
      else if (*f <= fold + ALF * alam * slope) return;
      else
        {
          if (alam == 1.0)
            tmplam = -slope / (2.0 * (*f - fold - slope));
          else
            {
              rhs1 = *f - fold - alam * slope;
              rhs2 = f2 - fold - alam2 * slope;
              a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
              b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
              if (a == 0.0) tmplam = -slope / (2.0 * b);
              else
                {
                  disc = b * b - 3.0 * a * slope;
                  if (disc < 0.0) tmplam = 0.5 * alam;
                  else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0 * a);
                  else tmplam = -slope / (b + sqrt(disc));
                }
              if (tmplam > 0.5 * alam)
                tmplam = 0.5 * alam;
            }
        }
      alam2 = alam;
      f2 = *f;
      alam = FMAX(tmplam, 0.1 * alam);
    }
}

void lubksb(float **a, int n, int *indx, float b[])
{
  int i, ii = 0, ip, j;
  float sum;

  for (i = 1; i <= n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii)
        for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
      else if (sum) ii = i;
      b[i] = sum;
    }
  for (i = n; i >= 1; i--)
    {
      sum = b[i];
      for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}
void lubksb1(float **a, int n, int *indx, float b[])
{
  int i, ii = 0, ip, j;
  float sum;

  for (i = 1; i <= n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii)
        for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
      else if (sum) ii = i;
      b[i] = sum;
    }
  for (i = n; i >= 1; i--)
    {
      sum = b[i];
      for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}


void ludcmp(float **a, int n, int *indx, float *d)
{
  int i, imax, j, k;
  float big, dum, sum, temp;
  float *vv;

  /* Initializes imax to remove warning message from the compiler.
     This initialization shall not change the method at all. */
  imax = 1;

  vv = vector(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++)
    {
      big = 0.0;
      for (j = 1; j <= n; j++)
        {
          if ((temp = fabs(a[i][j])) > big) big = temp;
        }
      if (big == 0.0)
        {
          nrerror("Singular matrix in routine ludcmp", 'b');
          /* BKP 0900 exit(1);*/
          return;

        }
      vv[i] = 1.0 / big;
    }
  for (j = 1; j <= n; j++)
    {
      for (i = 1; i < j; i++)
        {
          sum = a[i][j];
          for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0.0;
      for (i = j; i <= n; i++)
        {
          sum = a[i][j];
          for (k = 1; k < j; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((dum = vv[i] * fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (k = 1; k <= n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            }
          *d = -(*d);
          vv[imax] = vv[j];
        }
      indx[j] = imax;
      if (a[j][j] == 0.0) a[j][j] = TINY;
      if (j != n)
        {
          dum = 1.0 / (a[j][j]);
          for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
  free_vector(vv, 1, n);
}

void ludcmp1(float **a, int n, int *indx, float *d)
{
  int i, imax, j, k;
  float big, dum, sum, temp;
  float *vv;

  /* Initializes imax to remove warning message from the compiler.
     This initialization shall not change the method at all. */
  imax = 1;

  vv = vector1(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++)
    {
      big = 0.0;
      for (j = 1; j <= n; j++)
        {
          if ((temp = fabs(a[i][j])) > big) big = temp;
        }
      if (big == 0.0)
        {
          nrerror("Singular matrix in routine ludcmp1", 'b');
          /* BKP 0900 exit(1); */
          return;
        }
      vv[i] = 1.0 / big;
    }
  for (j = 1; j <= n; j++)
    {
      for (i = 1; i < j; i++)
        {
          sum = a[i][j];
          for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0.0;
      for (i = j; i <= n; i++)
        {
          sum = a[i][j];
          for (k = 1; k < j; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((dum = vv[i] * fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (k = 1; k <= n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            }
          *d = -(*d);
          vv[imax] = vv[j];
        }
      indx[j] = imax;
      if (a[j][j] == 0.0) a[j][j] = TINY;
      if (j != n)
        {
          dum = 1.0 / (a[j][j]);
          for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
  free_vector1(vv, 1, n);
}

void newt(float x[], int n, int *check,
          void (*vecfunc)(int, float [], float []))
{
  void fdjac(int n, float x[], float fvec[], float **df,
             void (*vecfunc)(int, float [], float []));
  float fmin0(float x[]);
  void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
              float * f, float stpmax, int * check, float(*func)(float []));
  void lubksb(float **a, int n, int * indx, float b[]);
  void ludcmp(float **a, int n, int * indx, float * d);
  int i, its, j, *indx;
  float d, den, f, fold, stpmax, sum, temp, test, **fjac, *g, *p, *xold;

  indx = ivector(1, n);
  fjac = matrix(1, n, 1, n);
  g = vector(1, n);
  p = vector(1, n);
  xold = vector(1, n);
  fvec = vector(1, n);
  nn = n;
  nrfuncv = vecfunc;
  f = fmin0(x); /* fvec[i] calculated */
  test = 0.0;
  for (i = 1; i <= n; i++)
    {
      if (fabs(fvec[i]) > test) test = fabs(fvec[i]);
    }
  if (test < 0.01 * TOLF)
    {
      *check = 0;
      FREERETURN
        }
  for (sum = 0.0, i = 1; i <= n; i++) sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt(sum), (float)n);
  for (its = 1; its <= MAXITS; its++)
    {
      fdjac(n, x, fvec, fjac, vecfunc);
      for (i = 1; i <= n; i++)
        {
          for (sum = 0.0, j = 1; j <= n; j++) sum += fjac[j][i] * fvec[j];
          g[i] = sum;
        }
      for (i = 1; i <= n; i++) xold[i] = x[i];
      fold = f;
      for (i = 1; i <= n; i++) p[i] = -fvec[i];
      ludcmp(fjac, n, indx, &d);
      if (bnrerrflag == 1) FREERETURN
                             lubksb(fjac, n, indx, p);
      lnsrch(n, xold, fold, g, p, x, &f, stpmax, check, fmin0);
      test = 0.0;
      for (i = 1; i <= n; i++)
        if (fabs(fvec[i]) > test) test = fabs(fvec[i]);
      if (test < TOLF)
        {
          *check = 0;
          FREERETURN
            }
      if (*check)
        {
          test = 0.0;
          den = FMAX(f, 0.5 * n);
          for (i = 1; i <= n; i++)
            {
              temp = fabs(g[i]) * FMAX(fabs(x[i]), 1.0) / den;
              if (temp > test) test = temp;
            }
          *check = (test < TOLMIN ? 1 : 0);
          FREERETURN
            }
      test = 0.0;
      for (i = 1; i <= n; i++)
        {
          temp = (fabs(x[i] - xold[i])) / FMAX(fabs(x[i]), 1.0);
          if (temp > test) test = temp;
        }
      if (test < TOLX)
        {
          FREERETURN
            }
    }
  nrerror("MAXITS exceeded in newt", 'b');
  FREERETURN
    }

void newt1(float x[], int n, int *check,
           void (*vecfunc)(int, float [], float []))
{
  void fdjac1(int n, float x[], float fvec1[], float **df,
              void (*vecfunc)(int, float [], float []));
  float fmin1(float x[]);
  void lnsrch1(int n, float xold[], float fold, float g[], float p[], float x[],
               float * f, float stpmax, int * check, float(*func)(float []));
  void lubksb1(float **a, int n, int * indx, float b[]);
  void ludcmp1(float **a, int n, int * indx, float * d);
  int i, its, j, *indx;
  float d, den, f, fold, stpmax, sum, temp, test, **fjac, *g, *p, *xold;

  indx = ivector1(1, n);
  fjac = matrix1(1, n, 1, n);
  g = vector1(1, n);
  p = vector1(1, n);
  xold = vector1(1, n);
  fvec1 = vector1(1, n);
  nn1 = n;
  nrfuncv1 = vecfunc;
  f = fmin1(x);
  test = 0.0;
  for (i = 1; i <= n; i++)
    if (fabs(fvec1[i]) > test) test = fabs(fvec1[i]);
  if (test < 0.01 * TOLF)
    {
      *check = 0;
      FREERETURN1
        }
  for (sum = 0.0, i = 1; i <= n; i++) sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt(sum), (float)n);
  for (its = 1; its <= MAXITS; its++)
    {
      fdjac1(n, x, fvec1, fjac, vecfunc);
      for (i = 1; i <= n; i++)
        {
          for (sum = 0.0, j = 1; j <= n; j++) sum += fjac[j][i] * fvec1[j];
          g[i] = sum;
        }
      for (i = 1; i <= n; i++) xold[i] = x[i];
      fold = f;
      for (i = 1; i <= n; i++) p[i] = -fvec1[i];
      ludcmp1(fjac, n, indx, &d);
      if (bnrerrflag == 1) FREERETURN1
                             lubksb1(fjac, n, indx, p);
      lnsrch1(n, xold, fold, g, p, x, &f, stpmax, check, fmin1);
      test = 0.0;
      for (i = 1; i <= n; i++)
        if (fabs(fvec1[i]) > test) test = fabs(fvec1[i]);
      if (test < TOLF)
        {
          *check = 0;
          FREERETURN1
            }
      if (*check)
        {
          test = 0.0;
          den = FMAX(f, 0.5 * n);
          for (i = 1; i <= n; i++)
            {
              temp = fabs(g[i]) * FMAX(fabs(x[i]), 1.0) / den;
              if (temp > test) test = temp;
            }
          *check = (test < TOLMIN ? 1 : 0);
          FREERETURN1
            }
      test = 0.0;
      for (i = 1; i <= n; i++)
        {
          temp = (fabs(x[i] - xold[i])) / FMAX(fabs(x[i]), 1.0);
          if (temp > test) test = temp;
        }
      if (test < TOLX)
        {
          FREERETURN1
            }
    }
  nrerror("MAXITS exceeded in newt1", 'b');
  FREERETURN1
    }

/***************NUMERICAL RECIPES UTILITY FUNCTIONS ******************/

void nrerror(char error_text[], char err_source)
/* Numerical Recipes standard error handler,
   modified to get error source BKP 08/01/00*/
{
  if (err_source == 'b')  bnrerrflag = 1;
  if (err_source == 'a')  anrerrflag = 1;
  if (err_source == 'g')
    {
      fprintf(stderr, "allocation error...now exiting to system...\n");
      exit(1);
    }
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v) nrerror("allocation failure in vector()", 'g');
  return v - nl + NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()", 'g');
  return v - nl + NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v = (unsigned char *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()", 'g');
  return v - nl + NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v = (unsigned long *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()", 'g');
  return v - nl + NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()", 'g');
  return v - nl + NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()", 'g');
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t)((nrow + NR_END) * sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()", 'g');
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t)((nrow + NR_END) * sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()", 'g');
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
                  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
  float **m;

  /* allocate array of pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()", 'g');
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
   and ncol=nch-ncl+1. The routine should be called with the address
   &a[0][0] as the first argument. */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure in convert_matrix()", 'g');
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl] = a - ncl;
  for (i = 1, j = nrl + 1; i < nrow; i++, j++) m[j] = m[j - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t = (float ***) malloc((size_t)((nrow + NR_END) * sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()", 'g');
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (float **) malloc((size_t)((nrow * ncol + NR_END) * sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()", 'g');
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (float *) malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()", 'g');
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
    {
      t[i] = t[i - 1] + ncol;
      t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
      for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
    }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG)(b + nrl - NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG)(b + nrl - NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

float *vector1(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v) nrerror("allocation failure in vector1()", 'g');
  return v - nl + NR_END;
}

int *ivector1(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v) nrerror("allocation failure in ivector1()", 'g');
  return v - nl + NR_END;
}

unsigned char *cvector1(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v = (unsigned char *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector1()", 'g');
  return v - nl + NR_END;
}

unsigned long *lvector1(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v = (unsigned long *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
  if (!v) nrerror("allocation failure in lvector1()", 'g');
  return v - nl + NR_END;
}

double *dvector1(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v) nrerror("allocation failure in dvector1()", 'g');
  return v - nl + NR_END;
}

float **matrix1(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix1()", 'g');
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix1()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix1(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t)((nrow + NR_END) * sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix1()", 'g');
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix1()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix1(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t)((nrow + NR_END) * sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix1()", 'g');
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix1()", 'g');
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **submatrix1(float **a, long oldrl, long oldrh, long oldcl, long oldch,
                   long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
  float **m;

  /* allocate array of pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix1()", 'g');
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix1(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
   and ncol=nch-ncl+1. The routine should be called with the address
   &a[0][0] as the first argument. */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t)((nrow + NR_END) * sizeof(float*)));
  if (!m) nrerror("allocation failure in convert_matrix1()", 'g');
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl] = a - ncl;
  for (i = 1, j = nrl + 1; i < nrow; i++, j++) m[j] = m[j - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor1(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t = (float ***) malloc((size_t)((nrow + NR_END) * sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor1()", 'g');
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (float **) malloc((size_t)((nrow * ncol + NR_END) * sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor1()", 'g');
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (float *) malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor1()", 'g');
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
    {
      t[i] = t[i - 1] + ncol;
      t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
      for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
    }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector1(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_ivector1(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_cvector1(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_lvector1(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_dvector1(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_matrix1(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_dmatrix1(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_imatrix1(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_submatrix1(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG)(b + nrl - NR_END));
}

void free_convert_matrix1(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG)(b + nrl - NR_END));
}

void free_f3tensor1(float ***t, long nrl, long nrh, long ncl, long nch,
                    long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}



/**************** NWET FUNCTIONS THAT USE DOUBLE VARS*******************/

void fdjac_double(int n, double x[], double fvec_double[], double **df,
                  void (*vecfunc)(int, double [], double []))
{
  int i, j;
  double h, temp, *f;

  f = dvector(1, n);
  for (j = 1; j <= n; j++)
    {
      temp = x[j];
      h = EPS * fabs(temp);
      if (h == 0.0) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      (*vecfunc)(n, x, f);
      x[j] = temp;
      for (i = 1; i <= n; i++) df[i][j] = (f[i] - fvec_double[i]) / h;
    }
  free_dvector(f, 1, n);
}

void fdjac1_double(int n, double x[], double fvec_double[], double **df,
                   void (*vecfunc)(int, double [], double []))
{
  int i, j;
  double h, temp, *f;

  f = dvector1(1, n);
  for (j = 1; j <= n; j++)
    {
      temp = x[j];
      h = EPS * fabs(temp);
      if (h == 0.0) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      (*vecfunc)(n, x, f);
      x[j] = temp;
      for (i = 1; i <= n; i++) df[i][j] = (f[i] - fvec_double[i]) / h;
    }
  free_dvector1(f, 1, n);
}

double fmin_double(double x[])
{
  int i;
  double sum;
  (*nrfuncv_double)(nn, x, fvec_double);
  for (sum = 0.0, i = 1; i <= nn; i++) sum += SQR(fvec_double[i]);
  return 0.5 * sum;
}

double fmin1_double(double x[])
{
  int i;
  double sum;

  (*nrfuncv1_double)(nn1, x, fvec1_double);
  for (sum = 0.0, i = 1; i <= nn1; i++) sum += SQR(fvec1_double[i]);
  return 0.5 * sum;
}

void lnsrch_double(int n, double xold[], double fold, double g[],
                   double p[], double x[], double *f, double stpmax, int *check,
                   double(*func)(double []))
{
  int i;
  double a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp,
    test, tmplam;

  /* These initializations have been added to remove warning messages
     from compilers. They don't change anything to the method as
     f2 and alam2 are initialized later in the loop before being used. */
  f2 = 0.;
  alam2 = 0.;

  *check = 0;
  for (sum = 0.0, i = 1; i <= n; i++) sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++) p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch_double.", 'a');
  test = 0.0;
  for (i = 1; i <= n; i++)
    {
      temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
      if (temp > test) test = temp;
    }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;)
    {
      for (i = 1; i <= n; i++) x[i] = xold[i] + alam * p[i];
      *f = (*func)(x);
      if (alam < alamin)
        {
          for (i = 1; i <= n; i++) x[i] = xold[i];
          *check = 1;
          return;
        }
      else if (*f <= fold + ALF * alam * slope) return;
      else
        {
          if (alam == 1.0)
            tmplam = -slope / (2.0 * (*f - fold - slope));
          else
            {
              rhs1 = *f - fold - alam * slope;
              rhs2 = f2 - fold - alam2 * slope;
              a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
              b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
              if (a == 0.0) tmplam = -slope / (2.0 * b);
              else
                {
                  disc = b * b - 3.0 * a * slope;
                  if (disc < 0.0) tmplam = 0.5 * alam;
                  else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0 * a);
                  else tmplam = -slope / (b + sqrt(disc));
                }
              if (tmplam > 0.5 * alam)
                tmplam = 0.5 * alam;
            }
        }
      alam2 = alam;
      f2 = *f;
      alam = FMAX(tmplam, 0.1 * alam);
    }
}

void lnsrch1_double(int n, double xold[], double fold, double g[],
                    double p[], double x[], double *f, double stpmax,
                    int *check1, double(*func)(double []))
{
  int i;
  double a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp, test, tmplam;

  /* These initializations have been added to remove warning messages
     from compilers. They don't change anything to the method as
     f2 and alam2 are initialized later in the loop before being used. */
  f2 = 0.;
  alam2 = 0.;

  *check1 = 0;
  for (sum = 0.0, i = 1; i <= n; i++) sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 1; i <= n; i++) p[i] *= stpmax / sum;
  for (slope = 0.0, i = 1; i <= n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch_double.", 'a');
  test = 0.0;

  for (i = 1; i <= n; i++)
    {
      temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
      if (temp > test) test = temp;
    }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;)
    {

      for (i = 1; i <= n; i++) x[i] = xold[i] + alam * p[i];

      /* BKP somehow *check1 is changed after the for statement */
      /* what if I reset it to 0  in each loop ? */
      (*check1) = 0;

      *f = (*func)(x);

      if (alam < alamin)
        {
          for (i = 1; i <= n; i++) x[i] = xold[i];
          *check1 = 1;
          return;
        }
      else if (*f <= fold + ALF * alam * slope)
        {
          return;
        }
      else
        {
          if (alam == 1.0)
            tmplam = -slope / (2.0 * (*f - fold - slope));
          else
            {
              rhs1 = *f - fold - alam * slope;
              rhs2 = f2 - fold - alam2 * slope;
              a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
              b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
              if (a == 0.0) tmplam = -slope / (2.0 * b);
              else
                {
                  disc = b * b - 3.0 * a * slope;
                  if (disc < 0.0) tmplam = 0.5 * alam;
                  else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0 * a);
                  else tmplam = -slope / (b + sqrt(disc));
                }
              if (tmplam > 0.5 * alam)
                tmplam = 0.5 * alam;
            }
        }
      alam2 = alam;
      f2 = *f;
      alam = FMAX(tmplam, 0.1 * alam);
    }
}

void lubksb_double(double **a, int n, int *indx, double b[])
{
  int i, ii = 0, ip, j;
  double sum;

  for (i = 1; i <= n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii)
        for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
      else if (sum) ii = i;
      b[i] = sum;
    }
  for (i = n; i >= 1; i--)
    {
      sum = b[i];
      for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}
void lubksb1_double(double **a, int n, int *indx, double b[])
{
  int i, ii = 0, ip, j;
  double sum;

  for (i = 1; i <= n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii)
        for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
      else if (sum) ii = i;
      b[i] = sum;
    }
  for (i = n; i >= 1; i--)
    {
      sum = b[i];
      for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}

void ludcmp_double(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  /* Initializes imax to remove warning message from the compiler.
     This initialization shall not change the method at all. */
  imax = 1;

  vv = dvector(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++)
    {
      big = 0.0;
      for (j = 1; j <= n; j++)
        {
          if ((temp = fabs(a[i][j])) > big) big = temp;

        }
      if (big == 0.0)
        {
          nrerror("Singular matrix in routine ludcmp_double", 'a');
          /* BKP 0900 exit(1) ; */
          return;
        }
      vv[i] = 1.0 / big;
    }
  for (j = 1; j <= n; j++)
    {
      for (i = 1; i < j; i++)
        {
          sum = a[i][j];
          for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0.0;
      for (i = j; i <= n; i++)
        {
          sum = a[i][j];
          for (k = 1; k < j; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((dum = vv[i] * fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (k = 1; k <= n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            }
          *d = -(*d);
          vv[imax] = vv[j];
        }
      indx[j] = imax;
      if (a[j][j] == 0.0) a[j][j] = TINY;
      if (j != n)
        {
          dum = 1.0 / (a[j][j]);
          for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
  free_dvector(vv, 1, n);
}

void ludcmp1_double(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  /* Initializes imax to remove warning message from the compiler.
     This initialization shall not change the method at all. */
  imax = 1;

  vv = dvector1(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++)
    {
      big = 0.0;
      for (j = 1; j <= n; j++)
        {
          if ((temp = fabs(a[i][j])) > big) big = temp;

        }
      if (big == 0.0)
        {
          nrerror("Singular matrix in routine ludcmp_double", 'a');
          /* BKP 0900 exit(1); */
          return;
        }
      vv[i] = 1.0 / big;
    }
  for (j = 1; j <= n; j++)
    {
      for (i = 1; i < j; i++)
        {
          sum = a[i][j];
          for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0.0;
      for (i = j; i <= n; i++)
        {
          sum = a[i][j];
          for (k = 1; k < j; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((dum = vv[i] * fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (k = 1; k <= n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            }
          *d = -(*d);
          vv[imax] = vv[j];
        }
      indx[j] = imax;
      if (a[j][j] == 0.0) a[j][j] = TINY;
      if (j != n)
        {
          dum = 1.0 / (a[j][j]);
          for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
  free_dvector1(vv, 1, n);
}

void newt_double(double x[], int n, int *check,
                 void (*vecfunc)(int, double [], double []))
{
  void fdjac_double(int n, double x[], double fvec_double[], double **df, void (*vecfunc)(int, double [], double []));
  double fmin_double(double x[]);
  void lnsrch_double(int n, double xold_double[], double fold, double g_double[], double p_double[], double x[], double * f, double stpmax, int * check, double(*func)(double []));
  void lubksb_double(double **a, int n, int * indx, double b[]);
  void ludcmp_double(double **a, int n, int * indx, double * d);
  int i, its, j, *indx;
  double d, den, f, fold, stpmax, sum, temp, test, **fjac_double, *g_double, *p_double, *xold_double;

  indx = ivector(1, n);
  fjac_double = dmatrix(1, n, 1, n);
  g_double = dvector(1, n);
  p_double = dvector(1, n);
  xold_double = dvector(1, n);
  fvec_double = dvector(1, n);
  nn = n;
  nrfuncv_double = vecfunc;
  f = fmin_double(x); /* fvec_double[i] calculated */
  test = 0.0;
  for (i = 1; i <= n; i++)
    {
      if (fabs(fvec_double[i]) > test) test = fabs(fvec_double[i]);
    }
  if (test < 0.01 * TOLF)
    {
      *check = 0;
      FREEDRETURN
        }
  for (sum = 0.0, i = 1; i <= n; i++) sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt(sum), (double)n);
  for (its = 1; its <= MAXITS; its++)
    {
      fdjac_double(n, x, fvec_double, fjac_double, vecfunc);
      for (i = 1; i <= n; i++)
        {
          for (sum = 0.0, j = 1; j <= n; j++) sum += fjac_double[j][i] * fvec_double[j];
          g_double[i] = sum;
        }
      for (i = 1; i <= n; i++) xold_double[i] = x[i];
      fold = f;
      for (i = 1; i <= n; i++) p_double[i] = -fvec_double[i];
      ludcmp_double(fjac_double, n, indx, &d);
      if (anrerrflag == 1) FREEDRETURN
                             lubksb_double(fjac_double, n, indx, p_double);
      lnsrch_double(n, xold_double, fold, g_double, p_double, x, &f, stpmax, check, fmin_double);
      test = 0.0;
      for (i = 1; i <= n; i++)
        if (fabs(fvec_double[i]) > test) test = fabs(fvec_double[i]);
      if (test < TOLF)
        {
          *check = 0;
          FREEDRETURN
            }
      if (*check)
        {
          test = 0.0;
          den = FMAX(f, 0.5 * n);
          for (i = 1; i <= n; i++)
            {
              temp = fabs(g_double[i]) * FMAX(fabs(x[i]), 1.0) / den;
              if (temp > test) test = temp;
            }
          *check = (test < TOLMIN ? 1 : 0);

          FREEDRETURN
            }
      test = 0.0;
      for (i = 1; i <= n; i++)
        {
          temp = (fabs(x[i] - xold_double[i])) / FMAX(fabs(x[i]), 1.0);
          if (temp > test) test = temp;
        }
      if (test < TOLX)
        {

          FREEDRETURN
            }
    }
  nrerror("MAXITS exceeded in newt_double", 'a');
  FREEDRETURN
    }

void newt1_double(double x[], int n, int *check1, void (*vecfunc)(int, double [], double []))
{
  void fdjac1_double(int n, double x[], double fvec1_double[], double **df, void (*vecfunc)(int, double [], double []));
  double fmin1_double(double x[]);
  void lnsrch1_double(int n, double xold_double[], double fold, double g_double[], double p_double[], double x[], double * f, double stpmax, int * check1, double(*func)(double []));
  void lubksb1_double(double **a, int n, int * indx, double b[]);
  void ludcmp1_double(double **a, int n, int * indx, double * d);
  int i, its, j, *indx;
  double d, den, f, fold, stpmax, sum, temp, test, **fjac_double, *g_double, *p_double, *xold_double;

  indx = ivector1(1, n);
  fjac_double = dmatrix1(1, n, 1, n);
  g_double = dvector1(1, n);
  p_double = dvector1(1, n);
  xold_double = dvector1(1, n);
  fvec1_double = dvector1(1, n);
  nn1 = n;
  nrfuncv1_double = vecfunc;
  f = fmin1_double(x);
  test = 0.0;
  for (i = 1; i <= n; i++)
    if (fabs(fvec1_double[i]) > test) test = fabs(fvec1_double[i]);
  if (test < 0.01 * TOLF)
    {
      *check1 = 0;
      FREEDRETURN1
        }
  for (sum = 0.0, i = 1; i <= n; i++) sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt(sum), (double)n);

  for (its = 1; its <= MAXITS; its++) /* begin iteration loop */
    {
      fdjac1_double(n, x, fvec1_double, fjac_double, vecfunc);
      for (i = 1; i <= n; i++)
        {
          for (sum = 0.0, j = 1; j <= n; j++) sum += fjac_double[j][i] * fvec1_double[j];
          g_double[i] = sum;
        }
      for (i = 1; i <= n; i++) xold_double[i] = x[i];
      fold = f;
      for (i = 1; i <= n; i++) p_double[i] = -fvec1_double[i];
      ludcmp1_double(fjac_double, n, indx, &d);
      if (anrerrflag == 1) FREEDRETURN1
                             lubksb1_double(fjac_double, n, indx, p_double);
      lnsrch1_double(n, xold_double, fold, g_double, p_double, x, &f, stpmax, check1, fmin1_double);

      test = 0.0;
      for (i = 1; i <= n; i++)
        if (fabs(fvec1_double[i]) > test) test = fabs(fvec1_double[i]);
      if (test < TOLF)
        {
          *check1 = 0;                   /* check1 = false, good return */
          FREEDRETURN1
            }
      if (*check1)
        {
          test = 0.0;
          den = FMAX(f, 0.5 * n);
          for (i = 1; i <= n; i++)
            {
              temp = fabs(g_double[i]) * FMAX(fabs(x[i]), 1.0) / den;
              if (temp > test) test = temp;
            }
          *check1 = (test < TOLMIN ? 1 : 0);
          FREEDRETURN1
            }
      test = 0.0;
      for (i = 1; i <= n; i++)
        {
          temp = (fabs(x[i] - xold_double[i])) / FMAX(fabs(x[i]), 1.0);
          if (temp > test) test = temp;
        }
      if (test < TOLX)
        {
          FREEDRETURN1
            }
    }
  nrerror("MAXITS exceeded in newt1_double", 'a');
  FREEDRETURN1
    }

#undef NRANSI

#undef ALF
#undef EPS
#undef MAXITS
#undef STPMX
#undef TOLX
#undef TINY
#undef TOLF
#undef TOLMIN
#undef TOLX


#undef FREERETURN
#undef FREERETURN1
#undef FREEDRETURN
#undef FREEDRETURN1
