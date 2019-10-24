C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
C     Author(s): Maryline Tombette
C
C     This file is part of the Simple Aqueous model (SIMPLE_AQUEOUS), a
C     component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

C     Aqueous-phase parameters and variables.

C     Number of gaseous species considered for the Simple Aqueous model.
      integer ns_aq
      parameter (ns_aq = 6)

C     Chosen index of species in the variable "gas".
      integer igso2, igh2o2, ignh3
      integer ighno3, igo3, igco2
      parameter (igso2  = 1)
      parameter (ignh3  = 2)
      parameter (ighno3 = 3)
      parameter (igh2o2 = 4)
      parameter (igo3   = 5)
      parameter (igco2  = 6)

C     Reference temperature for the comptation of Henry's constants.
      double precision temp_ref
      parameter(temp_ref = 298.d0)

C     Henry constants at 298K (in mol L-1 atm-1)
C     and heat of dissolution (in kcal mol-1)
C     for resp. SO2, NH3, HNO3, H2O2, O3 and CO2
      double precision chenry(ns_aq), dhhenry(ns_aq)
      data chenry /1.23d0, 62.d0, 2.1d5, 7.45d4, 1.13d-2, 3.4d-2/
      data dhhenry /-6.25d0, -8.17d0, -17.27d0, -14.5d0, -5.04d0,
     &             -4.85d0 /

C     Kinetic constants at 298K (in mol L-1) and delta H/R (in kcal mol-1)
C     For dissociation reactions (resp.):
C     H2O -> H+ + OH-
C     SO2.H2O -> HSO3- + H+
C     HSO3- -> H+ + SO3--
C     H2SO4(aq) -> H+ + HSO4-
C     HSO4- -> H+ + SO4--
C     NH3.H2O -> NH4+ + OH-
C     HNO3 -> NO3- + H+
C     CO2.H2O -> HCO3- + H+
C     HCO3- -> H+ + CO3--
      integer nreact_dissoc
      parameter(nreact_dissoc = 9)
      double precision ckdissoc(nreact_dissoc)
      double precision dhkdissoc(nreact_dissoc)
      data ckdissoc /1.d-14, 1.3d-2, 6.6d-8,  1.0d3, 1.02d-2,
     &     1.7d-5, 15.4d0, 4.3d-7, 4.7d-11/
      data dhkdissoc/-6710.d0, 1960.d0, 1500.d0, 0.d0, 2720.d0,
     &                  -450.d0, 8700.d0, -1000.d0, -1760.d0/

C     Kinetic constants at 298K and delta H
C     For sulfur oxydation reactions (resp.):
C     S(IV) + O3 -> S(IV) + O2 (3 constantes)
C     S(IV) + H2O2 -> S(VI) + H2O (2 constantes)
      integer nreact_oxydation
      parameter(nreact_oxydation = 5)
      double precision ckoxydation(nreact_oxydation)
      double precision dhkoxydation(nreact_oxydation)
      data ckoxydation / 2.4d4, 3.7d5, 1.5d9, 7.5d7, 13.0d0/
      data dhkoxydation / 0.0d0, 10.99d0, 10.48d0, 8.79d0, 0.0d0/

C     Activation diameter (dry) in \mu m.
      double precision dactiv
      parameter (dactiv = 0.2d0)

C     Separation diameter for two sections (dry) in \mu m.
      double precision dsep
      parameter (dsep = 1.5D0)

C     Parameters for bimodal activated distribution:
C     DM: mean diameter; SD: variance.
      double precision DM1_aq, DM2_aq, SD1_aq, SD2_aq
      parameter(DM1_aq = 0.4d0)
      parameter(DM2_aq = 2.5d0)
      parameter(SD1_aq = 1.8d0)
      parameter(SD2_aq = 2.15d0)

