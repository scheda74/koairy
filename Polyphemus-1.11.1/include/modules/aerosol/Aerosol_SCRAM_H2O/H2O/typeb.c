#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "nr.h"
#include "glodef.h"

/* global functions */
extern void newt1 (float x[], int n, int * chk, void (*vecfunc)(int, float[], float[]));
extern void thermob (double xx[], double ac[], int n);
extern void bpartitionfunc (int n, float xx[], float ff[]);
extern void Kpart (int n, double ac[], double MWom, double KB[], double VPB[], double HVAPB[], double K0B[], int oligoflag);
extern void bpartfuncapprox (int n, float aero[]);

/* global variables */
extern double PAOM;
extern double MWB[NBSP+1];
extern double KB[NBSP+1];
extern double VPB[NBSP+1];
extern double HVAPB[NBSP+1];
extern int Newtflag;
extern int aidxb[NBSP+1];
extern int oligoflag;
extern double K0B[NBSP+1];

/* Properties of primary organic absorbing medium, defined in glo.h */
extern double xaom[(NBSPAOM+1)]; /* mole fraction of non-volatile organics */
/* extern double fom; fraction of TSP that forms POA */
extern double MWaom;     /* mean MW or POA */
extern double MWaom_mix;
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(PAOM, MWB, KB, VPB, HVAPB, Newtflag, aidxb, xaom, \
                           MWaom, MWaom_mix)
#endif


/***************************************************************************
Purpose: Partition of Type B condensables, 
         Ai orig used in Unifac --> gamma, 
         gamma --> Ki; Ki --> the distribution of OC
         Particle phase OC (Ai final) calculated
         need f = Ai orig - Ai final = 0, f = delta aerosol is output

         1. the outer loop routine for newt, the globally convergent 
            multi-dimensional Newton's method, used to solve the 
            non-linear simultaneous equations.
         2. when newt is not used, gamma and K calculated based on Ai orig
            and partition is calculated for A without iterating

Preconditions: called by fmin (a subroutine of Newt) 
               or Main (when newt not used)

Key calls: unidriver, Kpart, partitionfunc, partfuncapprox	    

Notes: 1. defined (hardcoded in glo.h or glodef.h) values used: 
          fom = 0.1 (mass fraction of TSP forming AOM)
	  not used because PAOM is passed in
          MWaom (non-volatile) = 280
          NBSPAOM = 5,
          xaom = {0., 0.4, 0.05, 0.15, 0.12, 0.28};
	  default breakdown of AOM
	  compound			mass frac.	mole frac.
	  C24 alkanoic acid		0.50		0.40
	  C18 alkenoic acid		0.05		0.05
	  actonyl syringol		0.11		0.15
	  C20 alkane			0.17		0.12
	  arom. dicarboxylic acid	0.17		0.28
  
Revision history:  Developed by Betty Pun, AER, Jan 99 Under EPRI
                   Modified by Betty Pun, AER, Nov 99 Under CARB
                   1. increase the number of partitioning compound
                   2. conform to models-3 coding standard
		   3. allow the selection of equations to solve when 
		      using Newt.
		   4. added extern variable HVAPB for call to Kpart
***************************************************************************/
void TypeB (int n, float aero[], float f[])
{
  /* 
     n = no. of partitioning species
     aero = guess aero conc. from newt (main) or input aero conc.
     given x (aero) calculate ac and K, calculate x final
     (aerof) from partition.  f = aerof - aero = 0 when they agree 
  */

  double orgp[(NBSPAOM+1)];
  /* No. of moles of non-volatile organic absorbing compounds */
  float * aerof;                /* output concentration */
  float * aerotmp;  /* working array for unifac and Kpart with NBSP species */
  double tmolaom;               /* total mole absorbing organic medium */
  double tmol, tmolinv;         /* total mole, inverse total mole */ 
  double MWom;                  /* molecular weight of all organics */
  double * x;	
  /* mole fraction input into unidriver, index 0 .. n-1 */
  double * ac;
  /* activity coefficient calculated by unifac, index 0 .. n-1 */

  int i;
  int idx;
  int chk;

  aerotmp = calloc (NBSP+1, sizeof(float)); 
  aerof = calloc (n+1, sizeof(float));
  /* aerotmp has space for all species, but aerof has only non-zero ones */
  x = calloc (NBSPAOM+NBSP, sizeof(double));	
  ac = calloc (NBSPAOM+NBSP, sizeof(double));


  for (i = 1; i <= n; i ++){
    aerotmp[aidxb[i]] = aero[i];
    aerof[i] = aero[i];
  }
  
  /* Calculate mole fraction and MWom */
  /* tmolaom = PAOM / MWaom; */
  tmolaom = PAOM / MWaom_mix;

  for (i = 1; i <= NBSPAOM; i ++){
    orgp[i] = tmolaom * xaom[i]; /* umoles if input ug */
  }

  tmol = tmolaom;
  for (i = 1; i <= n; i ++){
    tmol += aerof[i] / MWB[aidxb[i]];
  }
  tmolinv = 1./ tmol;
  
  for (i = 1; i <= NBSPAOM; i ++)
    x[i-1] = orgp[i] * tmolinv;

  for (i = 1; i <= NBSP; i ++)
    x[NBSPAOM+i-1] = aerotmp[i] / MWB[i] * tmolinv;

  MWom = 0.;
  for (i = 0; i < NBSPAOM; i ++){
    MWom += x[i] * MWaom;
  }
  for (i = NBSPAOM; i < (NBSP + NBSPAOM); i ++){
    MWom += x[i] * MWB[i - NBSPAOM +1];
  }

  /* Call c Driver of fortran program */ 
  thermob (x, ac, (NBSPAOM + NBSP));
  
  /* calculate partition coefficient */
  Kpart(NBSP, &(ac[NBSPAOM-1]), MWom, KB, VPB, HVAPB, K0B, oligoflag);
  
  if (Newtflag == 1){
    /* calculate ai's simulataneously */
    newt1 (aerof, n, &chk, bpartitionfunc);

    if (chk != 0) 
      printf ("check not zero! Global minimum not found.\n");   

    for (i = 1; i <= n; i ++){
      f[i] = aerof[i] - aero[i];
    
    }
  } 
  else { 
    /* if Newt is not used, modify aero directly */
    bpartfuncapprox (n, aero);

  }

  
  free (x);
  free (aerof);
  free (aerotmp);
  free (ac);
  
}

