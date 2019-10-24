#include <stdio.h>
#include "nr.h"       /* NR convension index 1 .. n;
                         most arrays run from 1 to NBSP+1 */
#include "glodef.h"

/* global variables */
extern double PAOM;
extern double cb[NBSP + 1];
extern int Newtflag, bnrerrflag;
extern int aidxb[NBSP + 1];
extern double VPB[NBSP + 1];
extern double VPCrit;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(PAOM, cb, Newtflag, bnrerrflag, aidxb, VPB)
#endif

/**************************************************************************
Purpose:  Main for Type B module
          Calculates total condensables Ci
          Determines if set of simultaneous equations is solved.
          If so, calls Newt, the globally convergent multi-dimensional
          Newton's method to solve the non-linear simultaneous equations.
          Inner loop: partitionfunc (given Ci, Ki, solve Ai, Gi)
          Outer loop: iterate between Ai used in Unifac and Ai
                      calculated in partitionfunc (in this file)
          If not, partition calculated based on input aerosol mass
          and composition (fixed composition and mass of absorbing medium)

Key calls:  Newt, TypeB

Notes: 1. Parameters hard coded in glo.h (global variable include file)
          or glodef.h:
          used in bmain:
            NBSP, VPBi, MWBi of condensing species
          used in typeb:
            NBSPAOM, fom, MWBom, xaom[i]
       2. Parameters include file unifacparam.h used in unidriver
       3. Partitioning molecules:
          Anthropogenic 1 methyl nitro benzoic acid
    Anthropogenic 2 methyl hydroxy benzoic acid
    Biogenic 1 C15 diene aldehyde with OH and NO3 groups
    Biogenic 2 humulone aldehyde
    Biogenic 3 nopinone

Revision History:  Developed by Betty Pun, AER, Jan 99 Under EPRI
                   Modified by Betty Pun, AER, Nov 99 Under CARB
                   1. Increase the no. of partitioning compounds
                   2. Conform to models-3 coding standards
                   3. Bypass NEWT: so that partition is calculated based on
                      fixed absorbing organic phase amount and composition
                      (primary organic carbon + input SOA).  Newtflag = 0
                   4. Allow the selection of a subset of equations to solve
                      when using Newt.
                   5. Betty Pun, April 2000.  Combine Type A and Type B
                   6. BKP fix condition SOA >> PAOM
                   7. BKP July 2005, included new SOA compounds see #3 above
*************************************************************************/
void bmain(float gas[], float aero[])
{
  /* global functions */
  extern void newt(float x[], int n, int * check,
                   void (*vecfunc)(int, float[], float[]));
  extern void TypeB(int n, float x[], float f[]);

  FILE *fout1;
  int neq = -1;                    /* no. of equations to be solved by Newt */
  int check;                       /* flag used by newt */
  /* need to allocate space for check to avoid overwriting x[1] */
  int idx;                         /* counter for aidx */
  int i;                           /* dummy counter */
  double guesssoa;                 /* if SOA >> PAOM, may need to reset intial
                                      guess for newt */
  for (i = 1; i <= NBSP; i ++)
    {
      cb[i] = gas[i] + aero[i];
    }

  if (Newtflag == 1)  /* use Newton to solve simultaneous equations*/
    {

      guesssoa = 0.;

      /* from this point on in newt, aero contains only non-zero elements */
      idx = 1;
      for (i = 1; i <= NBSP; i ++)
        {
          /* assign initial guess for PM phase conc */
          if (cb[i] > PRAC_ZERO)
            {
              aidxb[idx] = i;
              if (PAOM == 0.)
                {
                  if (VPB[i] > VPCrit) aero[idx] = 0.3 * cb[i];
                  else aero[idx] = 0.9 * cb[i];
                }
              else
                {
                  /* this condition fails if SOA >> PAOM BPK 10/13/00*/
                  if (VPB[i] / PAOM > VPCrit)
                    {
                      aero[idx] = 0.3 * cb[i] ;
                      guesssoa += aero[idx];
                    }
                  else
                    {
                      aero[idx] = 0.9 * cb[i];
                      guesssoa += aero[idx];
                    }
                }
              idx ++;
            }
        }
      neq = idx - 1;

      /* bkp reset initial guess similar to PAOM = 0 if guesssoa >> PAOM */
      /* 10/13/00 */
      if (10 * PAOM < guesssoa)
        {
          for (idx = 1; idx <= neq; idx++)
            {
              if (VPB[aidxb[idx]] > VPCrit) aero[idx] = 0.3 * cb[aidxb[idx]];
              else aero[idx] = 0.9 * cb[aidxb[idx]];
            }
        }
      /*   end reset initial guess if guesssoa >> PAOM */
      if (neq > 0)
        newt(aero, neq, &check, TypeB);

      if (check != 0)
        {
          bnrerrflag = 1;
        }
      idx = 1;
      for (i = 1; i <= NBSP; i ++)
        {
          if (i == aidxb[idx])
            {
              /* 8/1/00 add error trap for negative conc */
              if (aero[idx] < 0.0)
                {
                  aero[idx] = 0.0;
                }
              else if (aero[idx] > cb[i])
                {
                  aero[idx] = cb[i];
                }

              gas[i] = cb[i] - aero[idx];
              idx ++;
            }
          else
            {
              /* with the use of prac_zero, gas[i] not necessarily 0.0
                 when equation not solved gas[i] = 0.0; */
              gas[i] = cb[i];
            }
        }  /* end for */

      /* reassign aero to contain zero elements for output in main */
      for (i = 1; i <= NBSP; i++)
        aero[i] = cb[i] - gas[i];

    }
  else                /* do not use newt */
    {

      /* all concentration elements passed into TypeB */
      for (i = 1; i <= NBSP ; i ++)
        {
          aidxb[i] = i;
        }

      TypeB(NBSP, aero, NULL);
      /* solution is output in aero array */

      for (i = 1; i <= NBSP; i ++)
        gas[i] = cb[i] - aero[i];

    }

  return;
}




