#include <stdio.h>
#include <stdlib.h>
#include "nr_double.h"
#include "glodef.h"
#include "binsolu.h"

double *aerozsr;  /* aero array used in zsr calculations only in this file */

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(aerozsr)
#endif

extern double totA[NAMOL+1], RH;
extern double MW[NAAERO+1];
extern int NK[NAMOL+1];
extern int naero;
extern int aidx[NAAERO+1];
extern int zsrflag;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(totA, RH, MW, NK, naero, aidx, zsrflag)
#endif
  
/**************************************************************************
Purpose: ZSR is used to calculate the amount of water associated with Type  
         A organic molecules only when zsrflag = 1; when zsrflag = 0 use 
	 Newt and UNIFAC to solve for Aw = RH

Preconditinos: called by TypeA routine

Subroutines called: Newt1_double (calls soalwcfunc)
                Newt1_double solves one implicit equation for a.c.(water) = RH

Revision history: 1. Developed by Betty Pun, AER, Feb 99 under EPRI funding
                     for the prototype SOA module.
		  2. Modified November 99 by Betty Pun, AER, under CARB
		     funding to adhere to models-3 coding standards
		  3. Modified to use ZSR or unifac to calculate water 
		     associated with organics, as specified by zsrflag
		     by Betty Pun, Nov, 99.  A file with xi at given Aw
		     is included binsolu.h
		  4. Fixed declaration of dLWC to have dimension of 2; 
		     since dLWC[1] is used.
***************************************************************************/
double ZSR (double *aerog)
{

  /* aerog contains the guesses of the aero species */

  void lwcorgfunc(int n, double x[], double f[]);
  extern void newt1_double(double x[], int n, int *check,void (*vecfunc)(int, double[], double[]));
 
  double molesolute;            /* total organic solute (umole)*/
  double dLWC[2];               /* water associated with organics ug/m3 */
  int i, j, chk, ieq, jeq;
  double aeromolec[NAMOL+1];     /* amt of organic solute (ions as solutes) */
                                /* in umole/m3 air, used in zsr calc */
  double molal1, molal2;        /* molalbin values bracketing RH of interest */
  double molalw;                /* linear interpolation of molal1 and molal2 */
  double binLWC[NAMOL+1];        /* LWC ug/m3 for binary solution */ 
  double totbinLWC;             /* sum ZSR LWC for binary solution */
  
  aerozsr = calloc( (naero+1), sizeof(double));
  /* 12/4/00 initialize chk */
  chk = 0;     
  for (i = 1; i <= naero; i ++)  {
    aerozsr[i] = aerog[i];
    
  }

  /* particulate water associated with dissolved species  
     ions treated as molecules.  Therefore, for the purpose  
     of calculating gamma, use 7 molecules 
  */

  if (zsrflag == 0) {

  /* guess Xw = RH */
  molesolute = 0.;
  for (i = 1; i <= naero; i ++){
    molesolute += aerozsr[i]/ MW[(aidx[i])] ;
    
  }


  dLWC[1] = molesolute * RH / (1.-RH)*MW[0];

  /* Call c Driver of fortran program to calculate solvent a.c. */ 
  /* activity of water solvent = RH */
  newt1_double (dLWC, 1, &chk, lwcorgfunc);


  /* if (chk != 0)
     printf("spurious convergence: gradient = 0\n"); */

  }
  else {                        /* zsrflag = 1 */
    ieq = 1;
    jeq = 1;

    /* find the umole/m3 air of the solutes including ions */
    for (i = 1; i <= NAMOL; i ++){
      aeromolec[i] = 0.0;
      if (aidx[ieq] == jeq) {
	for (j = 0; j < NK[i]; j ++) {
	  
	  aeromolec[i] += aerozsr[ieq]/MW[aidx[ieq]];
	  ieq ++;
	}
      }
      else aeromolec[i]= 0.0;
      jeq += NK[i];
    }                                    

    for (j = 1; j <= NAMOL; j ++){
      
      i = j;         
                     
      if (aeromolec[j] > 0.0) {

	/* look up umol j / ug water corresponding to RH */
	molal1 = molalbin[((int)(RH/RHgrad))][j-1];
	molal2 = molalbin[((int)(RH/RHgrad))+1][j-1];

	/* interpolate */
	
	molalw=molal1+(molal2-molal1)/RHgrad*(RH - RHgrad*((int)(RH/RHgrad))); 

	binLWC[j] =   aeromolec[j] / molalw;
	
      } /* end if */
      else binLWC[j] = 0.0;
    } /* end for */

    /* zsr equation */
    
    totbinLWC = 0.0;
    for (j = 1; j <= NAMOL; j ++) {
      if ((j < 1) || (j > NAMOL)) printf("j = %d\n", j);
      else {
	
	totbinLWC += binLWC[j];
      }
    }

    dLWC[1] = totbinLWC;
    
  }   /* end else zsrflag = 1 */
  free (aerozsr);
  return (dLWC[1]);
}


/**************************************************************************
Purpose: lwcorgfunc is the function used by NEWT1_double to determine 
         if activity of water equals RH in a given solution.

Preconditions: called by Newt1_double (which is called by ZSR)

Notes: the input water content associated with organics is given in dLWC[1]
       n is the no. of organic molecules
       f[1] is 100 x the deviation from RH = a.c.(water)

Revision History: 1. Developed by Betty Pun, AER, Feb. 99 under EPRI funding
                     for prototype SOA module with 2 compounds, malic acid 
		     and glyoxalic acid.
		  2. Modified by Betty Pun, AER, Nov. 99 under CARB funding
                     to partition as many as 6 organic compounds
		  3. Modified to comply with Models-3 coding standard, 
                     Nov. 99 

***************************************************************************/
void lwcorgfunc(int n, double dLWC[], double f[])
{
  extern void thermoa (double xx[], double gamma[], int n);
  /* unifac always called with 7 molecules, order as shown in main. */  

  double x[NAMOL+1], gamma[NAMOL+1], tmol, tmolinv;
  int i, j, ieq, jeq;

  tmol = dLWC[1] / MW[0];
  if (tmol < 0.0) tmol = 0.0;

  for (i = 1; i <= naero; i ++){
    tmol += aerozsr[i] / MW[(aidx[i])];
  }

  if (tmol != 0.) tmolinv = 1./ tmol; 
  else { 
    printf("Total Mole = 0 in ZSRsub! Exit.\n"); 
    exit (1); 
  }
  
   x[NAMOL] = dLWC[1] / MW[0] / tmol; 
   if (x[NAMOL] < 0.0) x[NAMOL] = 0.0;

  jeq = 1;
  ieq = 1;

  for (i = 1; i <= NAMOL; i ++){
    x[i-1]=0.0;
    if (aidx[ieq] == jeq) {
      for (j = 0; j < NK[i]; j ++) {
        x[i-1] += aerozsr[ieq]/MW[aidx[ieq]]/tmol;
        ieq ++;
      }
    }
    else x[i-1]= 0.0;
    jeq += NK[i];
  }                                                                            

  thermoa (x, gamma, NAMOL+1); 

  f[1] = 100. * (gamma[NAMOL] * x[NAMOL] / RH - 1.0);
  
}


