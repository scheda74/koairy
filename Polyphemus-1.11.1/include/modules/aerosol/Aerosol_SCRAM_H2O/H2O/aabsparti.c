#include <stdio.h>
#include <stdlib.h>
#include "glodef.h"

/* global variables */
extern double VPCrit;
extern int aidx[NAAERO+1], aidxmol[NAMOL+1];

extern double VPAtorr[NAMOL+1];

extern double totA[NAMOL+1];
/* input total concentration gas + particle phase */

extern double KA[NAMOL+1];
/* array of partition coefficients calculated by Kpart */

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(aidx, aidxmol, VPAtorr, totA, KA)
#endif

/***********************************************************************
Purpose: Equilibrium calculation for dry Type A absorptive partition 
         Calculate deviation from equilibrium given set of 
	 concentrations and partition constants
	 
Preconditions: called by fmin1_double, subroutine of newt1_double 

Revision History: Developed by Betty Pun, AER, Nov 00 under CARB funding
**************************************************************************/
void partiaabs (int n, double aerod[], double ff[]) {

  /* 
     aerod[1 .. (n-1)] are the PM phase concentration
     given aerod, calculate ff (deviation from equilibrium)
  */

  /* follow NR convension use index 1 .. n */
  double * a;
  double * g;
  double sumPM;  	
  /* sumPM = condensed OC (a[]) total amount absorbing phase */
  int i;

  a = calloc (n+1, sizeof(double));
  g = calloc (n+1, sizeof(double));

  sumPM = 0.;	
  
  for (i = 1; i <= n; i ++)    {
    a[i] = aerod[i];
    /* calculated sumPM = (fom * TSP + sum a[i]) */
    sumPM += a[i];
  }

  /* to calculate g and */	
  /* satisfy Kpart[i] = a[i] / sumPM / g[i] or mass balance */
  for (i = 1; i <= n; i ++){
    if (VPAtorr[aidx[i]] > VPCrit) {
      g[i] = totA[aidxmol[i]] - a[i];
      ff[i] = a[i] / sumPM / g[i] - KA[aidxmol[i]]; 
      /* add code to deal with negative gas conc when low TSP but
	 lots of condensable */
      if (g[i] <= 0) {
	g[i] = a[i]/sumPM/KA[aidxmol[i]];
	ff[i] = totA[aidxmol[i]] - g[i] - a[i];
      }
    }
    else {
      g[i] = a[i]/sumPM/KA[aidxmol[i]];
      ff[i] = totA[aidxmol[i]] - g[i] - a[i];
    } 
  }
  free (a);
  free (g);
	  
}

/***********************************************************************
Purpose: Approximated equilibrium calculation for Type B partition 
         Fixed absorbing medium (primary + secondary) amount and composition
	 used to calculate partition constants
	 
Preconditions: called by TypeAabs 

Notes: 1. all NAMOL species are passed in, therefore aidx/aidxmol not used.
 
Revision History: Developed by Betty Pun, AER, Feb 02 99 under CARB funding
		  as an alternative to solving simulataneous equations 
		  using newt. 

**************************************************************************/
void partiabsapprox (int n, double aerod[]) {
  /* aerod[1 .. n] are the PM phase concentration */

  double sumPM;  	
  /* sumPM = Minit + condensed OC (a[]) total amount absorbing phase */

  int i;

  sumPM = 0.;
	
  for (i = 1; i <= n; i ++)
      sumPM += aerod[i];
  
  /* to satisfy Kpart[i] = a[i] / sumPM / g[i]
     K * sumPM = a/(c-a)
     K * sumPM * c = a * K * sumPM + a
     a = K * sumPM * c / (1 + K * sumPM) 
   */
  
  for (i = 1; i <= n; i ++){
    aerod[i] = KA[i] * sumPM * totA[i] / (1 + KA[i] * sumPM);
  }
  
}

