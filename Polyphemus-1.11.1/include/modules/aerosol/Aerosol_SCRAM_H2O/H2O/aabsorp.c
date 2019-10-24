#include <stdio.h>
#include "nr.h"       /* NR convension index 1 .. n;
                         most arrays run from 1 to NASP+1 */
#include "glodef.h"

/* global variables */
extern double totA[NAMOL+1], g[NAMOL+1];
extern int anrerrflag, Newtflag;
extern int NK[NAMOL+1], aidx[NAAERO+1], aidxmol[NAMOL+1];
extern double VPAtorr[NAMOL+1];
extern double VPCrit; 

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(totA, g, anrerrflag, Newtflag, NK, aidx, aidxmol, \
                          VPAtorr)
#endif

/**************************************************************************
Purpose:  Absorption for partitioning Type A compounds with no existing water 
          If Newtflag = 1, use Newt_double, the globally convergent 
	  multi-dimensional Newton's method to solve the non-linear 
          simultaneous equations.
          Inner loop: partitionfunc (given Ci, Ki, solve Ai, Gi)
          Outer loop: iterate between Ai used in Unifac and Ai
                      calculated in partitionfunc (in this file)
          If Newtflag = 0, partition calculated based on input aerosol mass
          and composition (fixed composition and mass of absorbing medium)

Called by: amain 

Key calls:  Newt_double, TypeAabs
                                                          
Notes: 1. Parameters hard coded in glo.h (global variable include file)
          or glodef.h:
          used in bmain:
            NASP, VPi of condensing species
       2. Parameters include file unifacparam.h used in unidriver
       3. aabaero contains initial PM concentrations (used when newtflag = 0),
          and then gets changed to new concentrations.
Revision History:  1. Developed by Betty Pun, AER, DEC 00 Under CARB
                   2. Add code to calculate absorption by assuming fixed 
		      concentration of existing absorbing medium (not 
		      iterating) for 3-D application, Feb 01.
*************************************************************************/
void amainabs (double aabaero[] )
{
  /* global functions */
  extern void newt_double (double x[], int n, int * check, 
		    void (*vecfunc)(int, double[], double[]));
  extern void TypeAabs (int n, double x[], double f[]);    

  int neq;                         /* no. of equations to be solved by Newt */
  int check;                       /* flag used by newt */
  /* need to allocate space for check to avoid overwriting x[1] */     
  int idx;                         /* counter for aidx */
  int i,j;                         /* dummy counters */
  check = 0;

  if (Newtflag == 1) {
    for (i = 1; i <= NAMOL; i++) aabaero[i] = 0.0;
    /* from this point on in newt, aero contains only non-zero elements */
    idx = 1;
    j = 1;
    for (i = 1; i <= NAMOL; i ++) {
      /* assign indices initial guess for PM phase conc */
      if (totA[i] > PRAC_ZERO){
	aidxmol[idx] = i;
	aidx[idx] = j;
	if (VPAtorr[i] > VPCrit) aabaero[idx] = 0.3 * totA[i];
	else aabaero[idx] = 0.9 * totA[i];
	idx ++;
      }
      j += NK[i];
    }                                                                    
    neq = idx-1;
    
    newt_double (aabaero, neq, &check, TypeAabs);

    if (check != 0) {
      printf ("Type A absorp. check not zero! Global minimum not found.\n");
      anrerrflag = 1;
    }
    idx = 1;
    for (i = 1; i <= NAMOL; i ++){
      if (i == aidxmol[idx]){
        
	/* error trap for negative conc */
	if (aabaero[idx] < 0.0) {
	  aabaero[idx] = 0.0;
	}
	else if ( aabaero[idx] > totA[i]) {
	  aabaero[idx] = totA[i];    
	}
	
	g[i] = totA[i]- aabaero[idx]; 
	idx ++;
      }
      else {
	/* with the use of prac_zero, gas[i] not necessarily 0.0
	   when equation not solved gas[i] = 0.0; */
	g[i] = totA[i];
      }
    }  /* end for */
    /* reassign aabaero to contain all NASP (including totA = 0) 
     elements for output in main */
    for (i = 1; i <= NAMOL; i++) 
      aabaero[i] = totA[i] - g[i];
  } /* end if newtflag = 1 */
  else {
    /* newtflag = 0; all concentrations passed to Typeaabs */
    idx = 1;
    j = 1;
    for (i = 1; i <= NAMOL; i++) {
      aidxmol[idx] = i;
      aidx[idx] = j;
      idx ++;
      j += NK[i];
    }
    TypeAabs(NAMOL, aabaero, NULL);
    for (i = 1; i <= NAMOL; i++) g[i] = totA[i] - aabaero[i];
  }

  return;
}
