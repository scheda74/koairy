#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "nr.h"
#include "glodef.h"


/* global variables */
extern double MW[NAAERO + 1];
extern double KA[NAMOL + 1];
extern double VPAtorr[NAMOL + 1];
extern double totA[NAMOL + 1];
extern double HVAPA[NAMOL + 1];
extern double K0A[NAMOL + 1];
extern int aidx[NAAERO + 1], aidxmol[NAMOL + 1];
extern int NK[NAMOL + 1];
extern int Newtflag;
extern int oligoflag;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(MW, KA, VPAtorr, totA, HVAPA, aidx, aidxmol,	\
                          NK, Newtflag)
#endif

/* global functions */
extern void newt1_double(double x[], int n, int * chk, void (*vecfunc)(int, double[], double[]));
extern void ThermoA(double xx[], double ac[], int n);
extern void partiaabs(int n, double xx[], double ff[]);
extern void partiabsapprox(int n, double aerod[]);
extern void Kpart(int n, double ac[], double MWom, double KB[], double VP[], double HVAPA[], double K0A[], int oligoflag);

/***************************************************************************
Purpose: Partition of Type A condensables by absorption
         Ai orig used in Unifac --> gamma,
         gamma --> Ki; Ki --> the distribution of OC
         Particle phase OC (Ai final) calculated
         need f = Ai orig - Ai final = 0, f = delta aerosol is output

         1. the outer loop routine for newt_double, the globally convergent
            multi-dimensional Newton's method, used to solve the
            non-linear simultaneous equations.
   2. when newt is not used (Newtflag = 0) gamma and K calculated
            based on Ai orig and partition is calculated without iterations.

Preconditions: called by fmin_double (a subroutine of Newt_double)
               or by amainabs (aabsorp.c) if newt not used.

Key calls: unidriver, Kpart, partiaabs, partiabsapprox

Notes: 1. defined (hardcoded in glo.h or glodef.h) values used:

Revision history:  1. Developed by Betty Pun, AER, Nov 00 Under CARB
                   2. Added direct solution procedure (Newtflag = 0) Feb 01
       3. Added extern variable for Hvap to pass to Kpart
          July 2005
***************************************************************************/
void TypeAabs(int n, double aero[], double f[])
{
  /*
    n = no. of partitioning species
    aero = guess aero conc. from newt (main) or input aero conc.
    given x (aero) calculate ac and K, calculate x final
    (aerod) from partition.  f = aerod - aero = 0 when they agree
  */

  double * aerod;                /* output concentration */
  double * aerotmp;  /* work array for unifac and Kpart with NAMOL species */
  /* aerotmp has space for all species, but aerod has only non-zero ones */
  double tmol, tmolinv;         /* total mole, inverse total mole */
  double * x;
  /* mole fraction input into unidriver, index 0 .. n-1 */
  double * ac, * acshift;
  /* ac = activity coefficient calculated by unifac, index 0 .. n-1 */
  /* acshift = activity coefficients, index 1 .. n */
  int i, j;
  int idx;
  int chk1d;
  double MWom;
  int iteratenum;   /* set to iterate if Newtflag = 0 and no initial PM */

  aerotmp = calloc(NAMOL + 1, sizeof(double));
  aerod = calloc(n + 1, sizeof(double));
  x = calloc(NAMOL + 1, sizeof(double));
  ac = calloc(NAMOL + 1, sizeof(double));
  acshift = calloc(NAMOL + 1, sizeof(double));

  j = 1;
  for (i = 1; i <= n; i ++)
    {
      if (aidxmol[i] == j) aerotmp[j] = aero[i];
      j++;
      aerod[i] = aero[i];
    }

  iteratenum = 1;

  /* if Newtflag == 1, tmol > 0. because of initial guess
     tmol == 0 may happen when Newtflag = 0.
     In this case, we get an estimate of tmol and aerod by iteration
     using sequential substitution */

  if (Newtflag == 0)
    {
      tmol = 0.;
      for (i = 1; i <= n; i ++)
        {
          tmol += aerod[i] / MW[aidx[i]];
        }
      if (tmol == 0.0)
        {
          for (i = 1; i <= NAMOL; i++)
            {
              aerod[i] = 0.1 * totA[aidxmol[i]];
              aerotmp[i] = aerod[i];
              /* printf("aerod[%d] = %lf\n", i, aerod[i]); */
            }
          iteratenum = 3;
        }
    }

  while (iteratenum > 0)
    {

      tmol = 0.;
      for (i = 1; i <= n; i ++)
        {
          tmol += aerod[i] / MW[aidx[i]];
        }

      tmolinv = 1. / tmol;
      /*   printf("typeaabs: tmol = %le \n", tmol); */

      idx = 1;
      MWom = 0.;

      for (i = 0; i < NAMOL; i ++)
        {
          if (aerotmp[i + 1] >= 0.0)
            {
              x[i] = aerotmp[i + 1] / MW[aidx[idx]] / tmol;
              MWom += x[i] * MW[aidx[idx]];
              idx++;
            }
          else x[i] = 0.0;
        }

      /* Call c Driver of fortran program */
      ThermoA(x, ac, NAMOL + 1);

      /* calculate partition coefficient */
      for (i = 0; i < NAMOL; i++)
        {
          acshift[i + 1] = ac[i];
        }
      Kpart(NAMOL, acshift, MWom, KA, VPAtorr, HVAPA, K0A, oligoflag);

      if (Newtflag == 1)
        {
          /* calculate aerod[i] simulataneously */
          newt1_double(aerod, n, &chk1d, partiaabs);

          if (chk1d != 0)
            printf("check not zero! Global minimum not found in Aabs.\n");

          for (i = 1; i <= n; i ++)
            {
              f[i] = aerod[i] - aero[i];
            }
        }
      else
        {
          partiabsapprox(NAMOL, aerod);
          for (i = 1; i <= NAMOL; i++)
            {
              aerotmp[i] = aerod[i];
              aero[i] = aerod[i];
            }
        }
      iteratenum--;
    }

  free(x);
  free(aerod);
  free(ac);
  free(acshift);
  free(aerotmp);
}

