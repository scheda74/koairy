#include <stdio.h>
#include <stdlib.h>
#include "glodef.h"
#define maxtmpn 13
#define Critsol 100.0


/* global variables */
extern double totA[NAMOL + 1], RH, LWC, acHP;
extern double Keff[NAAERO + 1], MW[NAAERO + 1];
extern int NK[NAMOL + 1];
extern int aidx[NAAERO + 1];
extern double g[NAMOL + 1];
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(totA, RH, LWC, acHP, Keff, MW, NK, aidx, g)
#endif

/*********************************************************************

Purpose: Type A solute partition calculations
         given aero, calculate f (deviation from equilibrium)

Preconditions: Called by TypeA.

Revisions: 1. Developed by Betty Pun, AER, January, 99 under EPRI funding for
              prototype SOA module for 3 compounds, 5 equations.

     2. Modified by Betty Pun, June 99, to use aidex and n to
        select equations

     3. Modified November, 99 under CARB funding to accomodate 6
              compounds and select equations based on aindex and H. Highly
        soluble dissociative compounds: gas phase calculated by H
        and use mass balance equation.  Not so soluble compounds, gas
        phase calculated by mass balance and use H equation.

     4. Added code to fix problem with with negative gas
        concentrations, B. Pun January 2000.
        When compound not too soluble by Critsol but large LWC
        causes partition mostly into aq phase, switch back to
        partition equation for soluble species.
     5. bkp 6/00 change criteria to deal with low LWC
          old criterion was:
         if ((NK[i] >= 2) && (Keff[aidx[ieq]]/acHP > Critsol))
**********************************************************************/
void partfunc2(int n, double aero[], double gamma[], double f[])
{
  /*
    aero[1 .. n] are the PM phase concentration
    gamma[1 .. n] are the activity coefficients (Henry's Law)
  */

  int i, j, jj, ieq, jeq;          /* counters */
  double sumion;                   /* sum of all ions from one molecule */


  jeq = 1;
  i = 1;
  ieq = 1;

  while (ieq <= n)
    {

      if (aidx[ieq] == jeq)
        {

          for (j = 0; j < NK[i]; j++)
            {

              if (j == 0)                              /* solve Henry's equation */
                {
                  if (NK[i] >= 2 && Keff[aidx[ieq]] * LWC > Critsol)
                    {
                      /* most in aq phase */
                      /* gas-phase concentration determined by Henry's law */
                      g[i] = aero[ieq] * gamma[ieq] / LWC / Keff[aidx[ieq]];
                      sumion = 0.;
                      for (jj = 1; jj < NK[i]; jj ++)
                        sumion += aero[jj + ieq] * MW[aidx[ieq]] / MW[jj + aidx[ieq]];
                      f[ieq] = totA[i] - g[i] - aero[ieq] - sumion;
                    }
                  else
                    {
                      sumion = 0.;
                      for (jj = 1; jj < NK[i]; jj++)
                        sumion += aero[jj + ieq] * MW[aidx[ieq]] / MW[jj + aidx[ieq]];
                      /* gas phase by mass balance */
                      g[i] = totA[i] - aero[ieq] - sumion;
                      f[ieq] = aero[ieq] * gamma[ieq] / g[i] / LWC - Keff[aidx[ieq]];

                      /* 1/00 BP add the following to deal with -ve gas conc */
                      if (g[i] <= 0.0)
                        {
                          g[i] = aero[ieq] * gamma[ieq] / LWC / Keff[aidx[ieq]];
                          f[ieq] = totA[i] - g[i] - aero[ieq] - sumion;
                        }  /* end if */
                    } /* end else */
                  ieq ++;
                }   /* end henry's equation j = 0 */
              else
                {
                  /* to satisfy acid dissociation constants */
                  f[ieq] = acHP * aero[ieq] * gamma[ieq] / aero[ieq - 1] / gamma[ieq - 1] / Keff[aidx[ieq]] - 1.;
                  ieq ++;
                }   /* else acid dissociation */
            }     /* for j matched molecule */
        }       /* if idx[ieq] = jeq , condenable present*/
      jeq += NK[i];
      i++;
    }         /* end while ieq */
  return;
}





