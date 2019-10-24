#include <stdio.h>
#include <stdlib.h>
#include "glodef.h"

#define TINY2 1e-6
                                          

/* global variables and functions */
/* for sun implementation */
extern int thermoflag;
extern void bunidriver(double XPASS[], double GAMMA[], int n);

/**************************************************************************** 
Purpose: Call the chosen thermodynamic model 
   
****************************************************************************/  

void thermob (double XPASS[], double GAMMA[], int n)
{
  int i;
  double GAMMAinf[NBSP+NBSPAOM] = {1.0, 1.0, 1.0, 1.0, 1.0, 5.66, 2.29, 1.16, 1.16, 5.99, 3.81};
  double GAMMAf[NBSP+NBSPAOM] = {1.0, 1.0, 1.0, 1.0, 1.0, 5.66, 2.29, 1.16, 1.16, 5.99, 3.81};
  
  if (thermoflag==1)
	{
	  bunidriver(XPASS,GAMMA,n);
	}
  if (thermoflag==0)
	{
	  for (i=0; i<n; i++)
		{
		  if (XPASS[i] <= TINY2) /* if infinite dilution */
			{
			  GAMMA[i]=GAMMAinf[i];
				}
		  else
			{
			  GAMMA[i]=GAMMAf[i];
				}
		}
	}
  if (thermoflag==2) /*ideal*/
	{
	  for (i=0; i<n; i++)
		GAMMA[i]=1.0;
	}
  
  return;
}

#undef TINY2



