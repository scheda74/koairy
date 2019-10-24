#include <stdio.h>
#include <stdlib.h>
#define DIMFUN 10
/* max dimension of unifac parameter arrays for functional groups
   = NFUNC in unifacparam.h */

#define DIMMOL 10
/* max dimension of unifac parameter arrays for functional groups
   = NMOL in unifacparam.h */

#define TINY 1e-8
                                          

/* global variables and functions */

extern double temperature;
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(temperature)
#endif

/* for sun implementation */
extern void unifac_(int * NMOL, int * NFUNC, int NU[][DIMMOL], double X[], double A[][DIMFUN], double * RG, double * QG, double * Z, double * temperature, double * GAMA);

/* for ibm implmentation */
/*                                                                             
extern void unifac(int * NMOL, int * NFUNC, int NU[][DIMMOL], double X[], double A[][DIMFUN], double * RG, double * QG, double * Z, double * temperature, double * GAMA); 
*/

/* for c++ formulation: */
/* extern "C" {void unifac(int NMOL,int NFUNC, int NU[][50], double X[], \
   double A[][50], double * RG, double * QG, double Z, double temperature, \
   double * GAMA);}; 
*/


/**************************************************************************** 
Purpose: C (CPP) Driver for Unifac.for routine
         Translation of Pradeep Saxena's original driver 
         from Fortran to CPP

Preconditions: called by TypeA

Key Calls: Unifac

Notes: 1. Current formulation is for IBM machine
          Sun machine requires different "decorated name" to call 
	  fortran subroutine from C

	  extern void unifac_(int * NMOL, int * NFUNC, int NU[][50], \
	  double X[], double A[][50], double * RG, double * QG, \
	  double * Z, double * temperature, double * GAMA);    

	  unifac_(&NMOL,&NFUNC,nut,X,at,RG,QG,&Z,&temperature,GAMMA);   

       2. using ln linker (UNIX), fortran code needs address of variables and 
	  arrays.

       3. this routine uses indices of 0 .. n-1 and is the only one that 
          does not follow the NR tradition of indexing arrays using 1 .. n

       4. only one set of concentrations (X) are inputted 
          at each call.  RMOL is a one-dimension array.
	  Therefore, Nsolutions = No. of UNIFAC calls = 1.

       5. input X(i) = mole fraction of each component, 
          n = NSP, output GAMMA = ac.  File input UNIFAC parameters.  
	  No file output.

       6. takes 6 compounds

       7. xpass is used to receive x inputs, x is local and modified 
          (normalized) in unidriver.  while xpass, pointing to x, in typea
          is not changed because it continues to be used by typea.

       8. since C and fortran stores array in different order, we use original
          input order for A and NU and have to transpose them in C before
	  sending them into the fortran unifac routine
        

Revision History: Unidriver code translated by Betty Pun, AER, Apr 98 
                  under EPRI following P. Saxena's UNIFAC code. 
		  IMPLICIT REAL*8 (A-H,O-Z)
		  REAL*8 = double
		  variables whose names start with I, J, K, L, M, N
		  are INTEGER type.  Default INTEGER is INTERGER*4, 
		  INTEGER*4 = int.
		  error file will be opened in unifac subroutine
                  Amount of water is read from file, not fixed at 55.5 moles.

		  Incorporated into Type B calculations by Betty Pun, AER
		  Nov 98.  
		  
		  Modified by Betty Pun, Jan 99.  
		  Added renormalization if sum(xi) not 1

		  Modified by Betty Pun, December 99 to use an include file
		  instead of input read for unifac parameters

   
****************************************************************************/  
void unidriver (double XPASS[], double GAMMA[], int n)
{

#include "unifacparam.h"  /* uses DIMMOL and DIMFUN */

  FILE * fout1;
  int i,j;
  double X[DIMMOL];             /* local variable for x (normalized)*/	

  int nut[DIMFUN][DIMMOL];      /* transpose of NU */
  double at[DIMFUN][DIMFUN];    /* transpose of A */

  double sumx;                  /* sum of mole or mole fraction */

  if (NMOL != n){
    printf ("NMOL in parameter file does not match no. of solutes + water\n");
    exit (1);
  }

  if (NMOL != DIMMOL) {
    printf ("NMOL in parameter file does not match DIMMOL \n");
    exit (1);
  }

  if (NFUNC != DIMFUN) {
    printf ("NFUNC in parameter file does not match DIMFUNC\n");
    exit (1);
  }

  /*  transpose NU, A for passing to Fortran */

  for (i = 0; i < DIMMOL; i ++) {
    for (j = 0; j < DIMFUN; j ++) 
      
      nut[j][i] = NU[i][j];
  }

  for (i = 0; i < DIMFUN; i ++) {
    for (j = 0; j < DIMFUN; j ++) 
      at[j][i] = A[i][j];
   }  

  for (i = 0; i < n; i ++) X[i] = XPASS[i];

  /* check that x sums to 1 or renormalize */
  sumx = 0.0;
  for (i = 0; i < NMOL; i++){
    sumx += X[i];
    /* if (i == 3) printf("%le\t", X[i]);*/
  }

  if (((sumx - 1.0) > TINY) || ((sumx - 1.0) < -(TINY))){
    
    for (i = 0; i < NMOL; i ++){
      X[i] = X[i]  / sumx;
    }
  }
  
  unifac_(&NMOL,&NFUNC,nut,X,at,RG,QG,&Z,&temperature,GAMMA);

  /* for ibm implmentation */
  /* unifac(&NMOL,&NFUNC,nut,X,at,RG,QG,&Z,&temperature,GAMMA); */         

  return;
}

#undef TINY
#undef DIMFUN
#undef DIMMOL


