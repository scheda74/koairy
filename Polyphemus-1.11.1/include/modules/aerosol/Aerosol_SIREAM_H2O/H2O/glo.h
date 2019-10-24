/* BKP July 2005 global variables for the OA module */

/* Control Variables */

int zsrflag = 1;   /* flag to use binary solution and zsr if = 1 */
/* zsrflag = 0 calls unifac and solves implicit
   equation for a.c. water = RH using newt1 */

int Newtflag = 0;      /* flag to use NEWT in Type B module
                          and type A module with absorption when = 1,
                          if = 0; don't use NEWT */

int saturationflag = 0;  /* saturationflag = 1 means to
                            use saturation to determine particulate-phase
                            concentration when inorganic particle is dry.
                            If = 0, use absorption */
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(zsrflag, Newtflag, saturationflag)
#endif

/* General Variables */

double temperature;
double PAOM;            /* as input to Type B, PAOM defined as
                           non-volatile organic (other OC) */
double VPCrit = 1.0e-8; /* Criterion for setting intital particle
                           conc, compared to VP in torr or VP
                           rated by PAOM */
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(zsrflag, Newtflag, saturationflag)
#endif

/* Type A variables */

double totA[NAMOL + 1], acHP, RH, LWC;

double g[NAMOL + 1];   /* results of gas phase concentrations */
double negcharge;      /* total mole per m3 air of -ve charges */
double KA[NAMOL + 1];  /* dry Type A absorption partition constants */

int naero;           /* no. of species in aerosol phase: molecule + ions */
int aidx[NAAERO + 1]; /* list of molecular or ionic species present and equations to use) */
int aidxmol[NAMOL + 1]; /* list of molecular species present (equations to use in absorption) */
int anrerrflag;      /* numerical recipe error flag for Type A */
int oligoflag;       /* Flag for oligomerization. */
int thermoflag;      /* Flag for organic thermodynamic model. */

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(totA, acHP, RH, LWC, g, negcharge, KA, naero, aidx, \
                          aidxmol, anrerrflag)
#endif

/* Type A parameters */

double Critsol = 100.0;
double LWCTOL = 0.01;

int NK[NAMOL + 1] = {0, 3, 2, 1, 1, 1, 1, 1, 1, 2};  /* no. of eq. relationship for
                                                        each solute, in the order
                                                        listed in amain.c */

                                                     /*
                                                      * list of hydrophilic gaseous species: A2D, A1D, A0D, GLY (glyoxal),
                                                      *                                      MGLY (methylglyoxal), MT, PER, DER, MGA
                                                      * list of ionic/solute aerosol species: A2D, A2D-, A2D=, A1D, A1D-, A0D, AGLY,
                                                      *                                       AMGLY, MT, PER, DER, MGA, MGA-
                                                      */

double K[NAAERO + 1] = {0.0, 6.25e-3, 3.95e-4, 7.70e-6, 2.73e-3, 6.52e-4,
                        4.82e-5, 6.56e-4, 5.78e-12, 0.8052, 0.1109, 2.8,
                        1.1281e-2, 1.0e-4
};
/* partition parameters H and K
   H is in units of
   [(aq. microgram / m3 air)/(microgram water / m3 air)]
   / (g. microgram / m3 )
   estimated based on Suzuki et al., 1992
   K is in units of mol/kg water (same as {H+}) with
   concentrations of molecules in ions in the same mass-based units
   K's of malic acid and glyoxalic acid used respectively for
   compounds that dissociate twice and once */

double Keff[NAAERO + 1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0
};
/* effective partition parameters H and K depending of oligomerization,
   taken from Kroll & Seinfeld 2005 Keff = K * ( 1 + Ko,eff ),
   only A0D is affected, same units as above */
/* COMPUTED AT RUN TIME */

double pHoref = 6.0;
double Koeffref[NAAERO + 1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0
};
/* oligomerization constant reaction at reference pH = 6
   only A0D is affected Koeffref = [oligomer]/[monomer] = 10% */

double MW[NAAERO + 1] = {18.0, 186.0, 185.0, 184.0, 170.0, 169.0, 168.0, 58.0,
                         72.0, 136.0, 168.0, 136.0, 120.0, 119.0
};
/* Molecular weights in the order defined in amain.c */
/* MW[0] stores MW of water */
double DRH[NAMOL + 1] = {0.0, 0.79, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
/* deliquescence humidities for molecules */
/* see Peng et al., EST, 35: 4495, 2001 */

double VP[NAMOL + 1] = {0.0, 1.43, 1.98, 2.44e3, 6.86e8, 8.51e8, 10.7, 30.4,
                        3.80, 90.4
};
/* vapor pressure in units of mass (microgram) per m3 air */

double VPAtorr[NAMOL + 1] = {0.0, 1.43e-7, 2.17e-7, 2.7e-4, 219.8, 219.8,
                             1.45e-6, 2.61e-6, 4.1e-7, 1.4e-5
};
/* vapor pressure in torr */
double HVAPA[NAMOL + 1] = {0., 109.0, 50.0, 50.0, 25.0, 38.0, 38.4, 38.4, 38.4,
                           43.2
};
/* enthalpy of vaporization in kJ/mol */
double GAMMAinf[NAMOL] = {148.0, 224.0, 1201.0, 5.1, 4.6, 1.6, 6.5, 1.6, 1.6};
/* Activity coefficient at infinite dilution. */
double K0A[NAMOL + 1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.2};
/* Additional absorption due to oligomerization in the organic phase
   when no water Kp, eff = Kp * (1 + K0) */

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(Critsol, LWCTOL, NK, K, Keff, pHoref, Koeffref, MW, \
                          DRH, VP, VPAtorr, HVAPA)
#endif

/*
 * list of hydrophobic aerosol species :
 * AnBlP, AnBmP, BiBlP, BiBmP, BiNGA, trinitrate
 */

/* Type B variables */

double cb[NBSP + 1];    /* total Type B condensable concentrations */
double KB[NBSP + 1];    /* Type B partition constants */
int aidxb[NBSP + 1];    /* index for Type B equations to solve */
int bnrerrflag;         /* numerical recipe error flag for type b */
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(cb, KB, aidxb, bnrerrflag)
#endif

/* Type B parameters */

double MWB[NBSP + 1] = {0.0, 167.0, 152.0, 298.0, 236.0, 165.0, 272.0, 215.0};
double VPB[NBSP + 1] = {0.0, 6.8e-8, 8.4e-6, 6.0e-10, 3.0e-7, 1.39e-5, 1.45e-6,
                        2.5e-6
}; /* torr */

/* delta H vaporization added July 2005*/
double HVAPB[NBSP + 1] = {0.0, 50.0, 50.0, 175.0, 175.0, 43.2, 38.4, 50.0}; /* KJ/MOL */
double K0B[NBSP + 1] = {0.0, 0.0, 0.0, 0.0, 0.0, 64.2, 0.0, 0.0};
/* Additional absorption due to oligomerization in the organic phase
   Kp, eff = Kp * (1 + K0) */

/* the following are Type B absorbing medium parameters used in
   subroutine typeb */

double xaom[(NBSPAOM + 1)] = {0.0, 0.4, 0.05, 0.15, 0.12, 0.28};
/* mole fraction of non-volatile organics */

/* double fom = 0.1; */

double MWaom = 280.0;
double MWaom_mix;
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#pragma omp threadprivate(MWB, VPB, HVAPB, xaom, MWaom, MWaom_mix)
#endif
