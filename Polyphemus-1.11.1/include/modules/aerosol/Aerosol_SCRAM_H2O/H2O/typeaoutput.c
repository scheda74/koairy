/* preprocessor DEFINES and INCLUDES */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "glodef.h"   /* define statements for C preprocessor */
#include "glo.h"      /* external or global variables for the OA module */


void typeaoutput (int *tboaflag, double *aero, double *aeros, double deltaLWC)
{
  FILE * fout;
  /* global variables: totA[NAMOL+1], g[NAMOL+1], temperature, acHP, RH, LWC,
     negcharge declared in glo.h */
  int ieq, i, j;

 
  fout = fopen ("a.out" , "a");
  fprintf(fout, "ITS = %d\tRH = %lf\nLWC = %lf\t", (*tboaflag), RH, LWC);  
  fprintf(fout, "acHP = %lf\ntemperature = %lf\n", acHP, temperature); 
  fprintf(fout, "Total (microgram/m3)\tGas\tParticulate\tIonic Form(s)\n");

  if (LWC > LWC_ZERO) {
    /* counting molecules and ions in aero array */
    ieq = 1;   
    for (i = 1; i <= NAMOL; i ++){
      fprintf (fout, "%6.3le\t%6.3le\t", totA[i], g[i]); 
      printf ("%6.3le\t", totA[i]); /* if totA < PRACZERO gas can >0 */
      if (totA[i] > PRAC_ZERO){
	if (totA[i] != g[i]) {
	  if (aeros[i] != 0.0) { 
	    fprintf(fout, "%6.3le (dry)", aeros[i]);
	    printf("%6.3le (dry)", aeros[i]);
	  }	
	  else {
	    for (j = 0; j < NK[i]; j++) {
	      printf ("%6.3le\t", *(aero+ieq)); 
	      fprintf (fout, "%6.3le\t", aero[ieq]);
	      ieq ++; 
	    }                             /* end for */
	  }
	}     /* if totA != g */
      }                                 /* end if totA > 0 */
      fprintf(fout, "\n");  
    }                                   /* end for */                     
  }                                     /* if LWC > 0 */
  else {
    fprintf (fout, "no initial water, LWC = DLWC\n"); 
    printf ( "no initial water, LWC = DLWC\n"); 
    for (ieq = 1; ieq <= NAAERO; ieq ++) 
      printf("oamain debug: aero[%d]=%le\n", ieq, aero[ieq]); 
    ieq = 1;
    for (i = 1; i <= NAMOL; i ++) {
      fprintf (fout, "%le\t %le\t", totA[i], g[i]);  
      printf ("%le\t", totA[i]); /* if totA < PRACZERO gas > 0 */
      if (totA[i] > PRAC_ZERO){     
	if (totA[i] > g[i]) {     /* LWC = 0, g = totA or VP or absorp */
	  if (aeros[i] != 0.0){     /* dry particle */
	    fprintf (fout, "%le (solid)", aeros[i]); 
	    printf ("%le (solid)", aeros[i]); 
	  }
	  else {                    /* new aq phase */ 
	    fprintf (fout, "%le\t", aero[ieq]); 
	     printf ("%le\t", aero[ieq]); 
	    for (j = 1; j <= NK[i]; j ++) ieq ++;                    
	  }
	}
      }    /* if totA > 0 */
      fprintf(fout, "\n" ); 
      printf("\n" ); 
    } 
  }

  printf("total negative charges (micromole/m3): ");
  printf("%9.6le\n", negcharge);        

  fprintf (fout,"Water associated with organics (aq)\n");
  fprintf (fout, "%9.6le\n", (deltaLWC));

  printf ("Water associated with organics (aq): ");
  printf ("%9.6le\n", (deltaLWC));
  
  fclose (fout);           
}
