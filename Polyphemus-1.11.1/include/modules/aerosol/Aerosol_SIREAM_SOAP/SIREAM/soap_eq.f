C-----------------------------------------------------------------------
C     Copyright (C) 2017, CEREA- ENPC - EDF R&D
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
C     component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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

      SUBROUTINE SOAP_EQ(nesp_aer, nbin_aer,
     &     neq,q,lwcorg, lwc, rh, ionic, proton, 
     &     temp, aero, gas, liquid, psoap_config, psurrogate)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes ...
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     FLAG: whether to solved hydrophilic (=0) or hydrophobic (=1) species.
C     PROTON: hydronium ion concentration ([\mu g.m^-3]).
C     LWC: total liquid water content ([\mu g.m^-3]).
C     RH: relative humidity 0< <1 ([]).
C     TEMP: temperature ([Kelvin]).
C     IOLIGO: flag for oligomerization (true if =1)
C
C     -- INPUT/OUTPUT VARIABLES
C
C     AERO: aerosol bulk concentration ([\mu g.m^-3]).
C     GAS: gas concentration ([\mu g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C     ORGANION: organic ions ([\mu mol.m^-3]).
C     LWCORG: organic liquid water content ([\mu g.m^-3]).
C     CHP: hydronium ion concentration in water ([mol.L^-1]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2017: Youngseob Kim, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'time.inc'

      INTEGER nesp_aer, neq, nbin_aer
      double precision q(neq)
      DOUBLE PRECISION aero(nesp_aer),gas(nesp_aer)
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
      DOUBLE PRECISION liquid(12)

      INTEGER i,j
      integer psoap_config, psurrogate

      double precision DSD(nbin_aer),csol(nbin_aer)

      csol = 0.D0
      DSD = 0.D0    

c    Calculate the concentration of hydronium ion in water
c    microg/m3(=micromol/m3) / microg/m3 (H+ molar mass: 1 g/mol) 
c     = micromol/microg * 1000 
c     = mol/kg = mol/L (Water density: 1 kg/L) 
      chp = proton / lwc * 1.0e3

      CALL soap_main(lwc, rh, temp, ionic, chp, lwcorg, 
     &     psoap_config, psurrogate, 
     &     DT2, DSD, csol, liquid,
     &     nesp_aer, neq, q, aero, gas)

C     In case there is no gas-phase species.
C     For instance, CB05 mechanism doesn't have GLY for PGLY.
C     If gaseoues species don't exist, gas(j) can't be a gas-phase
C     concentration of the species and it must be set to zero.
      DO i = 1,nesp_aec
         j = aec_species(i)
         IF (aerosol_species_interact(j).LT.0) THEN
            aero(j) = aero(j) + gas(j)
            gas(j) = 0.0
         ENDIF
      ENDDO

      END
