/* preprocessor DEFINES and INCLUDES */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "glodef.h"   /* define statements for C preprocessor */
#include "glo.h"      /* external or global variables for the OA module */

void typeboutput(double PAOM, float gasb[], float aerob[])
{
  FILE *fout;
  int i;

  fout = fopen("b.out" , "a");

  fprintf(fout, "PAOM = %le\n", PAOM);
  fprintf(fout, " \t  gas  \t  aero \n");
  for (i = 1; i <= NBSP; i ++)
    {
      if (aerob[i] > 0.0)
        fprintf(fout, "%d\t%le\t%le\n", i, gasb[i], aerob[i]);
      else fprintf(fout, "%d\t%le\n", i, gasb[i]);
    }

  fclose(fout);
}
