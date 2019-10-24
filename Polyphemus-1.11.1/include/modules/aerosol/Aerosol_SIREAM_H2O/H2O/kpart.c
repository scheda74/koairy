#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


/* global variables */
extern double temperature;

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(temperature)
#endif

/*********************************************************************

Purpose: Kpart routine used in Type B particles calculation
         surrogate species approach with 5 condensing
         compound for Type B
         Reference: Pankow (1994) Atmos. Env. 28:189

         data needed:
         vapor pressure in Torr (subcooled if necessary)
   heat of vaporization in kJ/mol
         activity coefficients of condensing species
         calculated by Unifac, stored in global variable

         output KB[1..n] (particle phase/gas phase) in m3/ug

Preconditions: called by TypeB or Typeaabs

Revisions:  1. Developed by Betty Pun, AER, Jan 99 under EPRI for
               prototype SOA module with one compound octandecanoic acid
            2. Modified to comply with Models-3 coding standard, Betty Pun,
               AER, Nov 99, under ARB funding, incorporated into SOA module
               for 5 condensable compounds
      3. Modified July 2005 to add temperature correction
         use clausius-clapeyron equation.  Alternative method
         may be to estimate BP of surrogate and use
         Schwarzenbach et al. equation 4-33
***************************************************************************/
void Kpart(int n, double ac[], double mwom, double KB[], double VPB[], double HVAP[], double K0[], int oligoflag)
{
  /*
    parameters and constants used in Kpart
    may be moved to a data file
  */

  int i;
  float T;              /* local temperature variable */
  double R = 8.206e-5;  /* gas constant in units of m3-atm/mol/K */
  double RKJMOLK = 8.314e-3 ; /* gas constant in units of kJ/mol/K*/
  double Tref = 298.;   /* temperature at which group contribution method
                           applied to individual products */
  double VPCORR;

  T = temperature;
  for (i = 1; i <= n; i ++)
    {
      VPCORR = VPB[i] * exp(HVAP[i] * (1 / Tref - 1 / temperature) / RKJMOLK);
      KB[i] = 760 * R * 1e-6 * T / VPCORR / mwom / ac[i];
      if (oligoflag == 1)
        KB[i] *= 1 + K0[i];
    }
}




