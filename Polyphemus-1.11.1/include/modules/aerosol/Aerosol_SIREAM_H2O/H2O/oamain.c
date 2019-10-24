/* preprocessor DEFINES and INCLUDES */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "glodef.h"   /* define statements for C preprocessor */
#include "glo.h"

extern double PAOM;
extern double temperature, RH;
extern double negcharge;
extern double VP[NAMOL + 1] , VPAtorr[NAMOL + 1], HVAPA[NAMOL + 1], K[NAAERO + 1], \
  DRH[NAMOL + 1];
extern double VPB[NBSP + 1], HVAPB[NBSP + 1], K0A[NAMOL + 1], K0B[NBSP + 1];
extern int oligoflag;
extern int thermoflag;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(PAOM, temperature, RH, negcharge, VP, VPAtorr, \
                          HVAPA, K, DRH, VPB, HVAPB)
#endif

/**************************************************************************
Purpose:  Main for OA (Type A and B)
          Controls input/output

Arguments:
  Inputs : float *tempk    -- pointer to temperature in Kelvin
           float *rh       -- pointer to relative humidity in fraction
           float *worg     -- array of total g + p concentrations of 6 Type
                               A and 5 Type B aerosol compounds
                               in microgram / m**3
           float *gasorg   -- gas-phase concentration array
           float *partorg  -- particle-phase array (total mass concentration
                         of molecules and ions as molecules )
           (output gasorg and partorg should sum to worg)
           float *lwc      -- LWC in microgram / m**3
           float *protonconc - proton concentration in mole / g water
                  (microgram H+ / m3 air over microgram water / m3 air)
           int *tboaflag   -- flag to turn on type B oa module if == 1
             (type a is iterated, but no need to rerun type b)

  Outputs : float *DeltaLWC -- LWC associated with Type A OA in microgram / m3
            float *Organion -- organic anion concentration
                               (mole negativeve charge / m3)

Key calls: amain, bmain

Notes: 1. Type B input g and a, with option to calculate partition
          without using Newt, the iterative procedure by calculating
          K based on input PAOM and input a

       2. Parameters include file unifacparam.h used in unidriver

Revision History:  Prototype developed by Betty Pun, AER, Apr 00 Under EPRI

                   Further development by Betty Pun, AER, 2000-2001 Under
                   California ARB

                   Implemented in CMAQ, 2002 Under CARB

                   Extracted from CMAQ, June 2003

       July 2005 fixed aeros assignment when LWC > LWCZERO

Reference:  Pun, Griffin, Seigneur, Seinfeld, 2003, JGR, Vol 107, D17, 4333
            doi:10.1029/2001JD000542

****************************************************************************/
void oamain_(float *tempk, float *rh, float *worg, float *gasorg,
             float *partorg, float *mwaom_mix_loc, float *lwc,
             float *protonconc, float *Organion, float *DeltaLWC,
             int *tboaflag, int* ioligo, double *qsatref, double *tsatref,
             double *kpartref, double *drh, double *dhvap, int *ithermo)
{

  extern void amain(double * aero, double * aeros, double * LWCdelta2, int * ioligo);
  extern void bmain(float gasb[], float aerob[]);

  int i, j, ieq;                   /* counter */

  double aeros[NAMOL + 1]; /* solid species below DRH */
  double aero[NAAERO + 1]; /* aqueous species */
  double deltaLWC;        /* temporary variable for DeltaLWC */
  double *deltaLWC2;        /* temporary ptr variable for DeltaLWC */

  float aerob[(NBSP + 1)];
  /* type B particle concentrations ug/m3 air */
  /* for newt: initial guesses of Ai also contain final output */
  /* for direct solution: initial input Ai */
  float gasb[(NBSP + 1)];          /* type B gas-phase concentration */
  float tmp;         /* tmp var for partorg, gasorg type casted from double */
  FILE * fout;
  FILE * fout1;
  float *savpartorg;
  float *savlwc;

  VP[1] = qsatref[0];
  VP[2] = qsatref[1];
  VP[3] = qsatref[2];
  VP[4] = qsatref[3];
  VP[5] = qsatref[4];
  VP[6] = qsatref[5];
  VP[7] = qsatref[6];
  VP[8] = qsatref[7];
  VP[9] = qsatref[8];

  VPAtorr[1] = tsatref[0];
  VPAtorr[2] = tsatref[1];
  VPAtorr[3] = tsatref[2];
  VPAtorr[4] = tsatref[3];
  VPAtorr[5] = tsatref[4];
  VPAtorr[6] = tsatref[5];
  VPAtorr[7] = tsatref[6];
  VPAtorr[8] = tsatref[7];
  VPAtorr[9] = tsatref[8];

  K[1] = kpartref[0]; /* BiA2D */
  K[4] = kpartref[1]; /* BiA1D */
  K[6] = kpartref[2]; /* BiA0D */
  K[7] = kpartref[3]; /* AGLY */
  K[8] = kpartref[4]; /* AMGLY */
  K[9] = kpartref[5]; /* BiMT */
  K[10] = kpartref[6]; /* BiPER */
  K[11] = kpartref[7]; /* BiDER */
  K[12] = kpartref[8]; /* BiMGA */

  DRH[1] = drh[0];

  HVAPA[1] = dhvap[0] * 1.e-3;
  HVAPA[2] = dhvap[1] * 1.e-3;
  HVAPA[3] = dhvap[2] * 1.e-3;
  HVAPA[4] = dhvap[3] * 1.e-3;
  HVAPA[5] = dhvap[4] * 1.e-3;
  HVAPA[6] = dhvap[5] * 1.e-3;
  HVAPA[7] = dhvap[6] * 1.e-3;
  HVAPA[8] = dhvap[7] * 1.e-3;
  HVAPA[9] = dhvap[8] * 1.e-3;

  VPB[1] = tsatref[9];
  VPB[2] = tsatref[10];
  VPB[3] = tsatref[11];
  VPB[4] = tsatref[12];
  VPB[5] = tsatref[13];
  VPB[6] = tsatref[14];

  HVAPB[1] = dhvap[9] * 1.e-3;
  HVAPB[2] = dhvap[10] * 1.e-3;
  HVAPB[3] = dhvap[11] * 1.e-3;
  HVAPB[4] = dhvap[12] * 1.e-3;
  HVAPB[5] = dhvap[13] * 1.e-3;
  HVAPB[6] = dhvap[14] * 1.e-3;

  deltaLWC2 = calloc(sizeof(double), 1);

  savpartorg = partorg;

  temperature = (*tempk);
  RH = (*rh);
  LWC = (*lwc);
  acHP = (*protonconc);
  savlwc = lwc;
  oligoflag = (*ioligo);
  thermoflag = (*ithermo);

  /* initialize negcharge and deltaLWC */
  negcharge = 0.0;
  deltaLWC = 0.0;

  /* initialize error flag for numerical recipe */
  anrerrflag = 0;
  bnrerrflag = 0;

  if ((*tboaflag) == 0)
    {

      /* initialize aero */
      for (i = 1; i <= NAAERO; i ++) aero[i] = 0.0;

      /* initialize aeros */
      if (LWC < LWC_ZERO)
        for (i = 1; i <= NAMOL; i++) aeros[i] = partorg[i - 1];
      else
        for (i = 1; i <= NAMOL; i++) aeros[i] = 0.0;

      j = 0;           /* j counts the number of non-zero solute read */
      for (i = 1; i <= NAMOL; i ++)
        {
          totA[i] = worg[i - 1];
          if ((totA[i] > PRAC_ZERO) && (totA[i] < HUGE)) j++;
          else g[i] = totA[i];
          if (totA[i] >= HUGE)
            {
              fprintf(stderr, "totA[%d] too big, = %e\n", i, totA[i]);
              exit(1);
            }

          /* initial g for cases not requiring simultaneous solution, and not
             specified in amain.c (i.e., when 0 < totA < PRAC_ZERO */
        }

      if (j == 0)
        {
          /* printf("NO Type A organic solute in input\n"); */
          for (i = 1; i <= NAMOL; i ++) g[i] = totA[i];
        }
      else        /* add else to avoid calling Type A when no input 7/29/00 */
        {
          /* call Type A routines */
          amain(aero, aeros, deltaLWC2, ioligo);

          deltaLWC = *deltaLWC2;
        }

      partorg = savpartorg;
      if (anrerrflag == 0)
        {
          for (i = 1; i <= NAMOL; i++)
            {
              /* add negative conc error capture BKP 8/1/00 */
              if (g[i] < 0.0)
                g[i] = 0.0;
              else if (g[i] > totA[i])
                {
                  g[i] = totA[i];
                }
              /* partorg[i-1] = ((float) (totA[i]-g[i])); */
              tmp = ((float)(totA[i] - g[i]));
              partorg[i - 1] = tmp;
              tmp = ((float)(g[i]));
              /* gasorg[i-1] = ((float)(g[i])); */
              gasorg[i - 1] = tmp;

            }
          *DeltaLWC = ((float)(deltaLWC));
          *Organion = ((float)(negcharge));
        }
      else
        {
          /*  no change in gasorg, partorg for type A and deltalwc and organion = 0 */
        }

      /* Type A output */
      /* typeaoutput (tboaflag, aero, aeros, deltaLWC); */
    }

  if ((*tboaflag) == 1)
    {

      PAOM = worg[NAMOL + NBSP];
      MWaom_mix = (*mwaom_mix_loc);

      /* have to do full partition, can't do partapprox */

      for (i = (NAMOL + 1); i <= (NAMOL + NBSP); i ++)
        {
          gasb[i - NAMOL] = gasorg[i - 1];
          aerob[i - NAMOL] = partorg[i - 1];
        }

      bmain(gasb, aerob);
      /* inside bmain + subprograms, gas + aero are used, not gasb + aerob */

      /* Type B output: g[i], a[i] */
      /* typeboutput (PAOM, gasb, aerob); */

      if (bnrerrflag == 0)
        {
          for (i = NAMOL; i < (NAMOL + NBSP); i ++)
            {
              partorg[i] = ((float)(aerob[i - NAMOL + 1]));
              gasorg[i] = ((float)(gasb[i - NAMOL + 1]));
            }
        }     /* else don't change partorg and gasorg for type b */
    }       /* if tboaflag == 1 */
  free(deltaLWC2);
  return;
}
