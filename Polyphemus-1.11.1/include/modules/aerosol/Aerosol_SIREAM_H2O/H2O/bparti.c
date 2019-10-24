#include <stdio.h>
#include <stdlib.h>
#include "glodef.h"


/* global variables */
extern double VPCrit;
extern int aidxb[NBSP + 1];
extern double PAOM;
/* non-volatile primary PM for purpose of organic partition */

extern double VPB[NBSP + 1];

extern double cb[NBSP + 1];
/* input total concentration gas + particle phase */

extern double KB[NBSP + 1];
/* array of partition coefficients calculated by Kpart */

extern double fom;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(aidxb, PAOM, VPB, cb, KB, fom)
#endif


/***********************************************************************
Purpose: Equilibrium calculation for Type B partition
         Calculate deviation from equilibrium given set of
   concentrations and partition constants

Preconditions: called by fmin1, subroutine of newt1

Notes: 1. fom * TSP = Minit;  fom = 0.1 (hardcoded)
          In this module, TSP has no organic condensables
          Minit is amount non-volatile primary OC that serves
          as part of the absorbing organic phase.

       2. input PAOM = fom * TSP (pre-empts 1)

Revision History: Developed by Betty Pun, AER, Jan 99 under EPRI funding
                  Modified by Betty Pun, AER, Nov 99 under CARB funding
      1. to adhere to models-3 coding standard
      2. to solve selected equation based on aidx
      Added code to deal with negative gas concentration.
      B. Pun Jan, 2000.  This condition is encountered when
      a compound with high vapor pressure partitions
      significantly into the particle phase due to the presence
      of other SOA (rather than TSP).  To facilitate solution,
      equation to solve switched from Henry's law to mass
      balance.  Gas phase determined by Henry's law, since
      aerosol conc may exceed total input.
**************************************************************************/
void bpartitionfunc(int n, float aerof[], float ff[])
{

  /*
    aerof[1 .. (n-1)] are the PM phase concentration
    given aerof, calculate ff (deviation from equilibrium)
  */

  /* follow NR convension use index 1 .. n */
  double * a;
  double * g;
  double sumPM;
  /* sumPM = Minit + condensed OC (a[]) total amount absorbing phase */
  int i;

  a = calloc(n + 1, sizeof(double));
  g = calloc(n + 1, sizeof(double));

  sumPM = PAOM;

  for (i = 1; i <= n; i ++)
    {
      a[i] = aerof[i];

      /* calculated sumPM = PAOM + sum a[i]) */
      sumPM += a[i];
    }

  /* to calculate g and */
  /* satisfy Kpart[i] = a[i] / sumPM / g[i] or mass balance */
  for (i = 1; i <= n; i ++)
    {
      if (PAOM == 0.)
        {
          if (VPB[aidxb[i]] > VPCrit)
            {
              g[i] = cb[aidxb[i]] - a[i];
              ff[i] = a[i] / sumPM / g[i] - KB[aidxb[i]];
              /* add code to deal with negative gas conc when low TSP but
                 lots of condensable */
              if (g[i] <= 0)
                {
                  g[i] = a[i] / sumPM / KB[aidxb[i]];
                  ff[i] = cb[aidxb[i]] - g[i] - a[i];
                }
            }
          else
            {
              g[i] = a[i] / sumPM / KB[aidxb[i]];
              ff[i] = cb[aidxb[i]] - g[i] - a[i];
            }

        }
      else
        {
          if (VPB[aidxb[i]] / PAOM > VPCrit)
            {
              g[i] = cb[aidxb[i]] - a[i];
              ff[i] = a[i] / sumPM / g[i] - KB[aidxb[i]];
            }
          else
            {
              g[i] = a[i] / sumPM / KB[aidxb[i]];
              ff[i] = cb[aidxb[i]] - g[i] - a[i];
            }
        }

    }
  free(a);
  free(g);

}



/***********************************************************************
Purpose: Approximated equilibrium calculation for Type B partition
         Fixed absorbing medium (primary + secondary) amount and composition
   used to calculate partition constants

Preconditions: called by type B

Notes: 1. In this module, PAOM is primary absorbing organic material
          (non volatile) that is part of the absorbing organic phase.

Revision History: Developed by Betty Pun, AER, Nov 99 under CARB funding
      as an alternative to solving simulataneous equations
      using newt. Equations to be solved as selected based on
      aidx

**************************************************************************/
void bpartfuncapprox(int n, float aerof[])
{
  /* aerof[1 .. (n-1)] are the PM phase concentration
   */

  double sumPM;
  /* sumPM = PAOM + condensed OC (a[]) total amount absorbing phase */

  int i;

  sumPM = PAOM;

  for (i = 1; i <= n; i ++)
    sumPM += aerof[i];

  /* calculated sumPM = PAOM + sum a[i]) */

  /* to satisfy Kpart[i] = a[i] / sumPM / g[i]
     K * sumPM = a/(c-a)
     K * sumPM * c = a * K * sumPM + a
     a = K * sumPM * c / (1 + K * sumPM)
  */

  for (i = 1; i <= n; i ++)
    {
      aerof[i] = KB[aidxb[i]] * sumPM * cb[aidxb[i]] / (1 + KB[aidxb[i]] * sumPM);
    }

}

