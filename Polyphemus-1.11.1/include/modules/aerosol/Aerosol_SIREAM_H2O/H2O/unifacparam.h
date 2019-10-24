/*************************************************************************

Include file : unifacparma.h

Purpose: Unifac parameters (replace file read)

Include dependencies:  Included in Unidriver.c

Notes: Type A compounds + H2O at 25 C
       (Data of Lyman, Reehl, Rosenblatt, 1990)

       NMOL and NFUNC need to match DIMMOL and DIMFUN in unidriver.c

       Parameters inputted as in the original fortran input file, therefore
       a transpose is needed in the unidriver C program for matrices A and NU
       in order to pass them properly into the Fortran unifac routine.

       Orders of the functional groups:
         Groups                 Subgroups
         C                      CH3
         C                      CH2
         C                      CH1
         C                      C
         OH                     OH
         H2O                    H2O
         Carbonyls              CH3CO
         Aldehyde               HCO
         Carboxylic Acid        COOH
         Hydroxiperoxide        -O-O-H

revision History:  1. Developed by Betty Pun, AER, December, 1999
                under CARB funding
       2. Changed from original 5 Type A SOA + butandioic acid
                      to 6 anthropogenic and biogenic SOA for MADRID 1.5
          July 2005

**************************************************************************/



#ifndef UNIPARM_H
#define UNIPARM_H

/* no. of molecules */
int NMOL = 10;

/* no. of functional groups */
int NFUNC = 10;

/* Z = 10 is a fixed parameter in Unifac */
double Z = 10.0;

/* original file input has temperature,
   but temperature is in main input file now */

/* group volume parameters */
/* dimension of RG is the same as NFUNC */
double RG[DIMFUN] = {0.9011, 0.6744, 0.4469, 0.2195, 1.00, 0.92, 1.6724,
                     0.998, 1.3013, 1.3594
};

/* group surface area parameters */
/* dimension of QG is the same as NFUNC */
double QG[DIMFUN] = {0.8480, 0.5400, 0.2280, 0.0000, 1.20, 1.40, 1.4880,
                     0.948, 1.2240, 1.125
};

/* no. of groups in each molecule*/
int NU[DIMMOL][DIMFUN] =
  {
    {2, 2, 2, 1, 0, 0, 0, 0, 2, 0}, /* BiA2D */
    {2, 1, 2, 1, 0, 0, 1, 0, 1, 0}, /* BiA1D */
    {2, 2, 2, 1, 0, 0, 1, 1, 0, 0}, /* BiA0D */
    {0, 0, 0, 0, 0, 0, 0, 2, 0, 0}, /* Glyoxal */
    {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, /* Methyl Glyoxal */
    {1, 2, 1, 1, 4, 0, 0, 0, 0, 0}, /* BiMT */
    {1, 0, 1, 1, 2, 0, 0, 0, 0, 2}, /* BiPER */
    {1, 2, 1, 1, 4, 0, 0, 0, 0, 0}, /* BiDER */
    {1, 0, 1, 1, 2, 0, 0, 0, 1, 0}, /* BiMGA */
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}
  };

/* no. of groups in each molecule*/
double A[DIMFUN][DIMFUN] =
  {
    {0.0,      0.0    ,  0.0    ,  0.0    ,  986.5, 1318.00,  476.400,  677.000,  663.500, 977.56},
    {0.0,      0.0    ,  0.0    ,  0.0    ,  986.5, 1318.00,  476.400,  677.000,  663.500, 977.56},
    {0.0,      0.0    ,  0.0    ,  0.0    ,  986.5, 1318.00,  476.400,  677.000,  663.500, 977.56},
    {0.0,      0.0    ,  0.0    ,  0.0    ,  986.5, 1318.00,  476.400,  677.000,  663.500, 977.56},
    {156.4,  156.400,  156.400  , 156.400 ,  0.000,   353.5,   84.000,   -203.6,  199.00,  -330.28},
    {300.000,  300.000,  300.000,  300.000,  -229.10,   0.0,  -195.40,  -116.00,  -14.090, -314.18},
    {26.7600,  26.7600,  26.7600,  26.7600,  164.5, 472.500,  0.0    ,  -37.360,  669.400,  -350.58},
    {505.700,  505.700,  505.700,  505.700, 529.00, 480.800,  128.000,  0.0    ,  497.50,   -387.63},
    {315.300,  315.300,  315.300,  315.300, -151.00, -66.170,  -297.80,  -165.50,  0.0,     -501.23},
    { -23.233, -23.233 ,  -23.233,  -23.233, 342.92 , 795.55 , 380.94  , 408.88, 1479.0 , 0.0}
  };

#endif

/********************END unifacparam.h**********************************/


