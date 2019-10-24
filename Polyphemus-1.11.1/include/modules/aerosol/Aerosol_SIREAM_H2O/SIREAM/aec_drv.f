C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
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

      SUBROUTINE AEC_DRV(nesp_aer, flag, aero, gas, proton, lwc,
     &     organion, watorg, rh, temp, ioligo,
     &     ithermo)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the equilibrium between gas and
C     particle phase for oganic species using the AEC partioning
C     model (Pun et al 2001).
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
C     WATORG: organic liquid water content ([\mu g.m^-3]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2010/03/05: In case that there is no gas-phase species related,
C     See below (Youngseob KIM)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2007: Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'vara.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'varp.inc'

      INTEGER nesp_aer
      INTEGER flag, ioligo, ithermo
      DOUBLE PRECISION aero(nesp_aer),gas(nesp_aer)
      DOUBLE PRECISION aero_pom(nesp_pom), gas_pom(nesp_pom)
      DOUBLE PRECISION proton, lwc, organion, watorg, rh, temp
      DOUBLE PRECISION worg_poa,mwaom_mix

      REAL worg(nesp_aec + 1)
      REAL gasorg(nesp_aec), partorg(nesp_aec)
      REAL frh,ftempk,fdeltalwc,forganion
      REAL fprotonconc,flwc,fmwaom_mix

      INTEGER i,j

      DOUBLE PRECISION qsatref_loc(nesp_aec)
      DOUBLE PRECISION tsatref_loc(nesp_aec)
      DOUBLE PRECISION kpartref_loc(nesp_aec)
      DOUBLE PRECISION drh_loc(nesp_aec)
      DOUBLE PRECISION dhvap_loc(nesp_aec)

      DO i = 1,nesp_aec
         j = aec_species(i)
         qsatref_loc(i)=QSATREF(j)
         tsatref_loc(i)=TSATREF(j)
         kpartref_loc(i)=KPARTREF(j)
         drh_loc(i)=DRH(j)
         dhvap_loc(i)=DHVAP(j)
      ENDDO

C     zero init
      DO i=1,nesp_aec
         gasorg(i) = 0.0
         partorg(i) = 0.0
         worg(i) = 0.0
      ENDDO
      worg(nesp_aec + 1) = 0.0

      frh = REAL(rh)
      ftempk = REAL(temp)
      fdeltalwc = 0.0
      forganion = 0.0
      fprotonconc = 0.0
      flwc = 0.0
      fmwaom_mix = 0.0
      worg_poa = 0.0

C     total liquid water content
      IF (lwc.GT.0.d0) THEN
         flwc = REAL(lwc)
!     microg/m3(=micromol/m3) / microg/m3
!     = micromol/microg * 1000 = mole /kg = mol/L
         fprotonconc = REAL(proton / lwc * 1.0e3)
      ENDIF

C     concentrations in microg.m-3
      DO i = 1,nesp_aec
         j = aec_species(i)
         partorg(i) = REAL(aero(j))
         gasorg(i) = REAL(gas(j))

         worg(i) = partorg(i) + gasorg(i)
      ENDDO

      mwaom_mix = 0.0
      DO i = 1, nesp_pom
         j = poa_species(i)
         worg_poa = worg_poa + aero(j)
         mwaom_mix = mwaom_mix + EMW(j) * aero(j)
      ENDDO

C     primary organic mass + additional soa
      DO i=1,nesp_pankow        ! Pankow species
         j = pankow_species(i)
         worg_poa = worg_poa + aero(j)
      ENDDO

      worg(nesp_aec + 1) = REAL(worg_poa)

C     compute mean poa molweight in g/mol
      IF (worg_poa.GT.0.d0) THEN
         DO i=1,nesp_pankow
            j = pankow_species(i)
            mwaom_mix = mwaom_mix + EMW(j) * aero(j)
         ENDDO

         mwaom_mix = mwaom_mix / worg_poa
      ELSE
         mwaom_mix = EMW(poa_species(1))
      ENDIF
      mwaom_mix =  mwaom_mix * 1.d-06 ! from microg/mol to g/mol

      fmwaom_mix = REAL(mwaom_mix)

ccccccccccccccccccccccccccccccccccc
c     ftempk in K
c     frh 0< adim <1
c     worg in microg.m-3
c     gasorg in microg.m-3
c     partorg in microg.m-3
c     flwc in microg.m-3
c     fprotoconc in  mole / g water
c     fdeltalwc in microg.m-3
c     forganion in micromol.m-3
ccccccccccccccccccccccccccccccccccc

      CALL oamain(ftempk, frh, worg, gasorg, partorg, fmwaom_mix,
     &     flwc, fprotonconc, forganion, fdeltalwc, flag, ioligo,
     &     qsatref_loc, tsatref_loc, kpartref_loc, drh_loc, dhvap_loc,
     &     ithermo)

C     mol neg charge in micromol.m-3
      organion = DBLE(forganion)

C     AEC own liquid water content in microg.m-3
      watorg = DBLE(fdeltalwc)

C     Aqueous phase total liquid water content in microg.m-3
      lwc = lwc + DBLE(fdeltalwc)

C     Give back concentrations
      DO i = 1,nesp_aec
         j = aec_species(i)
         aero(j) = DBLE(partorg(i))
         gas(j) = DBLE(gasorg(i))
      ENDDO

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
