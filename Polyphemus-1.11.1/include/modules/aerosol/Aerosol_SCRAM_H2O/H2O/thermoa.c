#include <stdio.h>
#include <stdlib.h>
#include "glodef.h"

#define TINY2 1e-6
                                          

/* global variables and functions */
/* for sun implementation */
extern int thermoflag;
extern void unidriver(double XPASS[], double GAMMA[], int n);
extern double GAMMAinf[NAMOL];

/**************************************************************************** 
Purpose: Call the chosen thermodynamic model 
   
****************************************************************************/  

void thermoa (double XPASS[], double GAMMA[], int n)
{
  int i;
  
  if (thermoflag==1) /* unifac */
	  unidriver(XPASS,GAMMA,n);
	
  if (thermoflag==0) /*Gamma constants */
	{
	  for (i=0; i<n; i++)
		GAMMA[i]=GAMMAinf[i];
	}
  if (thermoflag==2) /*ideal*/
	{
	  for (i=0; i<n; i++)
		GAMMA[i]=GAMMAinf[i];
	}
  
  
  return;
}

#undef TINY2



