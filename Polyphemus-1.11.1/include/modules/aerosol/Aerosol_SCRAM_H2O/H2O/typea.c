
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "nr_double.h"
#include "glodef.h"
#define TINY 1.0e-8

/* global variables */
extern int naero;
extern int aidx[NAAERO+1];
extern double RH, acHP, totA[NAMOL+1], LWC; 
extern double MW[NAAERO+1];
extern int NK[NAMOL+1];
extern double negcharge;
extern double GAMMAinf[NAMOL];

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(naero, aidx, RH, acHP, totA, LWC, MW, NK, \
                          negcharge) 
#endif

/* global functions */
extern void thermoa (double xx[], double gamma[], int n);
extern void partfunc2 (int n, double aero[], double ac[], double ff[]);
/* extern void partitionfunc(int n, double aero[], double ac[], double f[]); */

/***************************************************************************
Purpose: Type A module: takes input of particle concentrations and returns
         deviations from equilibrium

Preconditions: this subroutine is called from newt, the globally 
               convergent multi-dimensional Newton's method used to 
               solve the non-linear simultaneous equations.

Subroutines called: unidriver and partition
                    (unidriver is always called with all 7 molecules that may
		    be present) 
		     
Revisions: 1. Developed by Betty Pun, AER, Jan 99, under EPRI funding for 
              prototype SOA module with 2 condensable compounds (malic acid,
	      glyoxalic acid) and water

	   2. Modified November 99 to accept 6 model compounds + water under
	      CARB funding (for the list of condensables see main)

	   3. Modified treatment of negative test concentrations so that
	      program does not exit prematurely.  B. Pun Jan 2000.
**************************************************************************/
void TypeA (int n, double aero[], double f[])
{
  FILE * fout1;

  double tmolaom, tmol, tmolinv, MWom;
  double tmpx;
  double  x[NAMOL+1];	
  /* mole fraction input into unidriver, index 0 .. n */
  double gammar[NAMOL+1] ;
  /* activity coefficient calculated by unifac, index 0 .. n */
  double gamma[NAAERO+1];
  /* index 1 .. n and used in partition calculation 
     gamma[1 ..n] = activity coefficients Henry's Law
     gamma contains 1 for activity coefficients of ions
  */

  int i,j, ieq, jeq;

 /* x = calloc ((NAMOL+1), sizeof(double));	
    gamma = calloc ((naero+1), sizeof(double));
    gammar = calloc ((NAMOL+1), sizeof(double));  */

  for (i = 0; i <= NAMOL; i++) {
    x[i] =0.;
    gammar[i] = 0.;
  }
  for (i = 0; i <= NAAERO; i++) gamma[i] = 0.;


  /* Calculate mole fraction (molecules only) */
  tmol = LWC / MW[0];
  for (i = 1; i <= naero; i ++){
    /* deal with negative aero conc */
    if (aero[i] >= 0.0)
    {
      tmol += aero[i] / MW[(aidx[i])];
    }
  }


  if (tmol > 0.) tmolinv = 1./ tmol;
  else {
    printf("Total Mole <= 0 in Type A! Skip.\n");
    negcharge = 0.0;
    exit (1);
  }
  /* 
     particulate water associated with dissolved species  
     ions treated as molecules.  Therefore, for the purpose  
     of calculating gamma, use seven molecules: 
     (0-5) 6 solute species, (6) water 
  */

  x[NAMOL] = LWC / MW[0] / tmol;/*mole fraction of H2O*/
  jeq = 1;
  ieq = 1;

  for (i = 0; i < NAMOL; i ++)
  {
    if (aidx[ieq] == jeq)
    {
      for (j = 0; j < NK[i+1]; j ++)
      {
	if (aero[ieq] >= 0.0)
	{
	  x[i] += aero[ieq]/MW[aidx[ieq]]/tmol;
	  /*if(aero[ieq]>HUGE)
	    printf("x[%d](%e)=aero[%d](%e)/MW[%d](%e)/tmol(%e)\n",i,x[i],
		   ieq,aero[ieq],aidx[ieq],MW[aidx[ieq]],tmol);*/
	}
	ieq ++;
      }
    }
    else x[i]= TINY;
    jeq += NK[i+1];
  }
  
/*  if(aero[ieq]>1e+4)
    printf("TypeA: middle\n\n");*/
  
  for ( i = 0; i <= NAMOL; i ++) {
    if (x[i] >= 0.0){
    }
    else {
      /* if (x[i] < 0.0) { */
      printf("TypeA: x[%d] < 0 (= %e).  Skip...\n", i, x[i]);
      negcharge = 0.0;
      exit (1);
    }
  }

  /* Call c Driver of fortran program to calculate solute a.c. */ 
  thermoa (x, gammar, NAMOL+1);/*x should be unchanged, why check if after*/
  for ( i = 0; i < NAMOL; i ++) gammar[i]=gammar[i]/GAMMAinf[i];
    
  for ( i = 0; i <= NAMOL; i ++) {
    if (x[i] >= 0.0){
    }
    else {
      /* if (x[i] < 0.0) { */
      printf("TypeA: x[%d] < 0 (= %e).  Exiting...\n", i, x[i]);
      negcharge = 0.0;
      exit (1);
    }
  }
  jeq = 1;
  ieq = 1;
  for (i = 0; i < NAMOL; i ++) {
    if (jeq == aidx[ieq]){
      gamma[ieq] = gammar[i];      /* first assign Raoult's law gamma */
      ieq += NK[i+1];
    }
    jeq += NK[i+1];
  }

  /* calculate solute activity coefficient at infinite dilution 
     according to the standard state of Raoult's Law 
  */
  jeq = 1;
  ieq = 1;
  for (i = 0; i < NAMOL; i ++) {
	/* used when 0.0 was used  if (x[i] > TINY) { */
    if (jeq == aidx[ieq]){
	  /* try setting anion a.c. the same as molecule 6-1-99 */
      for (j = 1; j < NK[i+1]; j ++) {
		gamma[ieq+j] = gamma[ieq];  /* gamma[2] = 1.0; */
		/* gamma[3] = 1.0;*/
      }
      ieq += NK[i+1];
    }
    jeq += NK[i+1];
  }

  /* calculate deviation from equilibrium (solutes) */
  partfunc2(naero, aero, gamma, f);


  /* calculate negcharge */

  negcharge = 0.0;
  ieq = 1;
  jeq = 1;     
  for (i = 1; i <= NAMOL; i ++){
    if (aidx[ieq] == jeq) {
      for (j = 0; j < NK[i]; j ++) {
	negcharge += aero[ieq]*j/MW[aidx[ieq]];
	ieq ++;
      }
    }
    jeq += NK[i];
  }
  
  /*printf("negcharge=%e \n",negcharge);*/
}


#undef TINY
