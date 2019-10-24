#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "nr_double.h"
#include "glodef.h"

/* global variables */
extern double RH, LWC, temperature;
extern double totA[NAMOL + 1], acHP, Critsol, LWCTOL;
extern int naero, zsrflag, anrerrflag, saturationflag;
extern double g[NAMOL + 1];
extern int NK[NAMOL + 1] , aidx[NAAERO + 1];
extern double VP[NAMOL + 1] , MW[NAAERO + 1], K[NAAERO + 1], DRH[NAMOL + 1];
extern double Keff[NAAERO + 1], Koeffref[NAAERO + 1], pHoref;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(RH, LWC, temperature, totA, acHP, Critsol,	\
                          LWCTOL, naero, zsrflag, anrerrflag,		\
                          saturationflag, g, NK, aidx, VP, MW,		\
                          K, DRH, Keff, Koeffref, pHoref)
#endif

/***********************************************************************
Purpose: Starting point of Type A routine.  Performs 3 functions: (1)
         gas/particle partition of organic compounds based on available
   water or saturation, (2) calculate water associated with organic
   solutes, (3) calculate total organic anion concentrations,

Arguments: 1. *aero: initial guesses of Ai in same units as inputs ug/m3
       also contain final output from newt
           2. *aeros: solid concentrations (also used to temporarily store
                   input PM concentrations in LWC = 0 absorption case with
       no newton, and PM from absorption before water is added)

Return:    1. deltaLWC is the output for water associated with the organics
        in microgram/m3 air

Data Needed:
   Input read in main and passed unsing global variables
          totA[i] (ug/m3 air), acHP (mole/kg water),
                LWC (ug/m3 air), RH (between 0 and 1), temperature (K)

   Parameters: NK = no. of eq. relationship for each solute
               K = partition parameters H and K
                   Ai, Gi, W in ug per m3 air units
                         gamma in mole fraction units, change reference
       state to Henry's law
       {H+} in moles per kg solvent,
         MW = Molecular weights in the order defined below
              borrowing MW[0] to store MW(water)
         DRH = deliquescence humidities for molecules
               (still need data)
         VP = vapor pressure in units of mass (ug) per m3 air
              pure liquid-gas partitioning:
                    (1) deliquescent species
                          (2) input LWC = 0  from inorganic species 5-29-99

Key functional calls: newt_double; TypeA; ZSR; amainabs; saturation

Notes:
    Molecules (6: 1-3 are anthropogenic; 4-6 are biogenic products)
    1. Dimethyl hydroxy hexadiendioic acid (2 dissociations)
    2. C8 acid with aldehyde and carbonyl groups (single dissiciation)
    3. C7 aldehyde with carbonyl and hydroxy groups (non-dissociative)
    4. Pinic acid (2 dissociations)
    5. Norpinonic acid (single dissociation)
    6. C10 aldehyde with carbonyl (non-dissociative)

    Solutes (12)
    1. Dimethyl hydroxy hexadiendioic acid (H2A1)
    2. Dimethyl hydroxy hexadienoate acid (HA1-)
    3. Dimethyl hydroxy hexadiendioate (A1=)
    4. C8 acid with aldehyde and carbonyl groups (HA2)
    5. C8 carboxylate with aldehyde and carbonyl groups (A2-)
    6. C7 aldehyde with carbonyl and hydroxy groups (A3)
    7. Pinic acid
    8. Hydrogen pinate (HA4-)
    9. Pinate (A4=)
    10. Norpinonic acid (HA5)
    11. Norpinonate (A5-)
    12. C10 aldehyde with carbonyl (A6)

    Cases:
    1. LWC > 0
       A) RH > DRH[i] - aqueous phase (with ions)
       B) RH < DRH[i]
          (i)  totA > VP - gas phase = VP
                         - solid phase (no water associated with PM phase i)
          (ii) totA < VP - gas phase only
    2. LWC = 0 (option 1, saturation)
       A) RH > DRH[i]
          (i)  totA > VP - gas phase = VP
                   - aq phase molecules only = totA[i] - VP[i]
          (ii) totA < VP - gas phase only
       B) RH < DRH[i]
          (i)  totA > VP - gas phase = VP
                         - solid phase (no water associated with PM phase i)
          (ii) totA < VP - gas phase only
    2. LWC = 0 (option 2, absorption)
       A) RH > DRH - gas/liquid partition, liquid phase
                     associated with organic water (molecules, no ions)
       B) RH < DRH - gas/liquid partition, no water in liquid phase

Revisions: 1. Developed by Betty Pun, AER, Jan 99 under EPRI for prototype
              Type A module with 2 compounds: malic acid and glyoxalic acid
        using newt, the globally convergent multi-dimensional
        Newton's method to solve the non-linear simultaneous equations.
              NR convension index 1 .. n; most arrays run from 1 to NSP+1

     2. Added code May 99 to deal with LWC input = 0
        calculate gas-PM partition based on VPsat
        calculate water associated with organics (H2M RH < DRH)

     3. Under CARB funding, modified code October 99 to
        perform the partition of 6 compounds.  Removed hard-wired
              code regarding equations solved

     4. Modified to comply with Models-3 coding standard, Betty Pun,
        AER, November 99

     5. Included flag to solve for water content based on ZSR with
        binary solution characteristics stored in file.

     6. Combine Type A and B, Betty Pun Apr 00
     6a.   bkp 6/00 change criteria to deal with low LWC cases:
     old criteria was:
     if ((NK[i] >= 2) && (Keff[jeq]/acHP > Critsol)) {
     7. Add option to do absorption when LWC = 0 11/00
     8. For 3-D, add option to solve absorption based on fixed PM

     9. Changed surrogate compounds, notes updated July 2005.
****************************************************************************/
/* double amain (double *aero, double *aeros) */
void amain(double *aero, double *aeros, double *deltaLWC2, int* ioligo)
{
  FILE * fout1;
  /* external subroutines */
  extern void newt_double(double x[], int n, int * check,
                          void (*vecfunc)(int, double[], double[]));
  extern void TypeA(int n, double x[], double f[]);
  extern double ZSR(double x[]);
  extern void amainabs(double aeros[]);
  void saturation(double tot, double vp, double * pmconc, double * gasconc);
  extern double pow(double x, double y);

  int check;
  /* need to allocate space for check to avoid overwriting x[1] */
  int i, j;              /* dummy counter reused in several places */
  int ieq;               /* ieq counts the selected equations */
  int jeq;               /* jeq loops through all equations */
  double deltaLWC;       /* water associated with organics ug/m3 */

  double *jnk;

  /*
   * oligomerization correction on H & K partition parameters
   */

  double cHP = 1. * acHP;
  double cHPref = pow(10.0, -pHoref);
  /* H+ concentration in mol/L, needed to correct oligomerization constant */

  double correcHP = pow(cHP / cHPref, 1.91);
  /* Koeffref correction due to pH */

  /* if oligomerization, modifies Keff */
  if (*ioligo == 1)
    {
      for (i = 1; i <= NAAERO; i ++)
        Keff[i] = K[i] * (1.0 + Koeffref[i] * correcHP);
      /* only A0D is affected */
    }
  else
    for (i = 1; i <= NAAERO; i ++) Keff[i] = K[i];



  /* initialize check 7/29/00 */
  check = 0;

  /* initalize aidx regardless of LWC 9/19/00 */
  for (i = 1; i <= NAAERO; i ++) aidx[i] = 0;

  /* figure out what species and initial conditions */
  if (LWC > LWC_ZERO)
    {
      for (i = 1; i <= NAMOL; i++) aeros[i] = 0.0;
      naero = 0;
      ieq = 0;
      jeq = 0;
      for (i = 1; i <= NAMOL; i ++)
        {
          if ((totA[i] > PRAC_ZERO) && (RH > DRH[i]))
            {
              for (j = 1; j <= NK[i]; j ++)
                {
                  aidx[ieq + 1] = jeq + j;
                  ieq ++;
                }
            }
          jeq += NK[i];
        }                           /* end for i */

      if (ieq > 0) naero = ieq;

      /* initial guess for aq. phase conc of organics */

      ieq = 1;
      jeq = 1;  /* this time, jeq points at the first K related to species i */

      for (i = 1; i <= NAMOL; i ++)
        {
          if ((totA[i] > PRAC_ZERO) && (RH > DRH[i]))
            {
              for (j = 0; j < NK[i]; j++)
                {
                  if (j == 0)             /* Henry's law equation */
                    {
                      if ((NK[i] >= 2) && (Keff[jeq]*LWC > Critsol))
                        {

                          /* twice dissocaitve and high H: assume very soluble */
                          /* initial guess based on most solute in water */
                          /* let computer manipulate gas */
                          if (NK[i] == 3)
                            aero[ieq] = 0.999999 * totA[i] / (1.0 + MW[aidx[ieq]] /
                                                              MW[aidx[ieq + 1]] * Keff[jeq + 1] / acHP + MW[aidx[ieq]] /
                                                              MW[aidx[ieq + 2]] * Keff[jeq + 1] * Keff[jeq + 2] / acHP / acHP);
                          if (NK[i] == 2)
                            aero[ieq] = 0.999999 * totA[i] / (1.0 + MW[aidx[ieq]] /
                                                              MW[aidx[ieq + 1]] * Keff[jeq + 1] / acHP);
                        }
                      else                   /* not very soluble species */
                        {
                          aero[ieq] = totA[i] * Keff[jeq] * LWC / (1. + Keff[jeq] * LWC);
                        }
                    }                       /* end Henry's equation if j == 0 */
                  else
                    {
                      /* acid dissociation of a dissociative compound */
                      aero[ieq] = Keff[jeq + j] * aero[ieq - 1] / acHP;
                    }                       /* end if dissociation */
                  ieq ++;
                }                         /* end for j */
            }                           /* end if aqueous */
          jeq += NK[i];
        }                             /* end for i */

      if (naero != 0)
        {
          newt_double(aero, naero, &check, TypeA);
        }
      if (check != 0)
        {
          /* check not zero! Global minimum not found. */
          anrerrflag = 1;
        }
      else          /* newt solution found */
        {
          /* need to call TypeA again to calculate gamma and gas phase */
          jnk = calloc((naero + 1), sizeof(double));
          TypeA(naero, aero, jnk);
          free(jnk);
        }                                             /* end if check not zero */
      /* end LWC > 0 aqueous phase (totA > 0 and RH>DRH) */

      /* LWC > 0 but compounds < DRH are not associated with water: saturation */
      for (i = 1; i <= NAMOL; i ++)
        {
          if (RH < DRH[i]) saturation(totA[i], VP[i], &(aeros[i]), &(g[i]));
        }
    }
  else                                             /* LWC = 0.0 */
    {
      if (saturationflag == 1)
        {
          ieq = 0;
          jeq = 1;
          for (i = 1; i <= NAMOL; i ++)
            {
              if (totA[i] > PRAC_ZERO)
                {
                  if (RH > DRH[i])
                    {
                      saturation(totA[i], VP[i], &(aero[ieq + 1]), &(g[i]));
                      if (aero[ieq + 1] > 0.0)      /* new aq phase */
                        {
                          ieq ++;
                          aidx[ieq] = jeq;

                          for (j = 1; j < NK[i]; j ++)
                            {
                              ieq ++;
                              aero[ieq] = 0.0;
                              aidx[ieq] = jeq + j;
                            }                                   /* end for j */
                        }
                    }      /* end if RH > DRH */
                  else saturation(totA[i], VP[i], &(aeros[i]), &(g[i]));
                  /* check if dry particle formed */
                }                                         /* if totA > 0 */
              jeq += NK[i];
            }                                           /* end for i */
        }                                             /* end if saturationflag */
      else
        {
          amainabs(aeros);
          /* reset aidx for wet solution */
          for (ieq = 1; ieq <= NAAERO; ieq ++) aidx[ieq] = 0;
          ieq = 0;
          jeq = 1;
          for (i = 1; i <= NAMOL; i ++)
            {
              if ((RH > DRH[i]) && (aeros[i] > 0.0))    /* new aq phase */
                {
                  ieq ++;
                  aero[ieq] = aeros[i];
                  aeros[i] = 0.0;
                  aidx[ieq] = jeq;

                  for (j = 1; j < NK[i]; j ++)
                    {
                      ieq ++;
                      aero[ieq] = 0.0;
                      aidx[ieq] = jeq + j;
                    }                                     /* end for j */
                }                                       /* end if */
              jeq += NK[i];
            }                                         /* end for i */
        }
      naero = ieq;

    }                                                     /* LWC = 0.0 */

  /* calculate water content associated with Type A */
  /* call zsr only if no numerical error in solution 9/20/00 */

  if ((naero != 0) && (anrerrflag == 0))
    {
      /* take care of RH = 1 if zsrflag = 0 */
      if (zsrflag == 0)
        {
          if (RH >= 1.0) RH = 0.99;
        }

      deltaLWC = ZSR(aero);

      /* LWC in ug/m3, delta LWC in ug/m3 */

    }
  else deltaLWC = 0.0;

  *deltaLWC2 = deltaLWC;

}



void saturation(double tot, double vp, double *pmconc, double *gasconc)
{

  /**************************************************************
   purpose: given total amount of a compound,
            calculate partition based on saturation
            if tot > vp, gas = vp, pm = (tot - vp)
            if tot < vp, gas = tot pm = 0
   arguments: tot      total amount (microgram/m3 air) of compound
              vp       vapor pressure (microgram/m3 air) of compound
              *pmconc  pointer to output of particulate-phase concentration
              *gasconc pointer to output of gas-phase concentration
   history: 1. coded 11/28/00 BKP to replace code in 2 locations:
               - LWC > 0 and RH < DRH
         - LWC = 0 any RH

  ****************************************************************/


  if (tot > vp)
    {
      *pmconc = tot - vp;
      *gasconc = vp;
    }
  else
    {
      *gasconc = tot;
      *pmconc = 0.0;
    }
}
