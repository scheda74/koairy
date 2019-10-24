C-----------------------------------------------------------------------
C     Copyright (C) 2012-2013, ENPC
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

      SUBROUTINE POA_DRV(nesp_aer,aero, gas, kpart, temp, dhvap)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the equilibrium between gas and particle
C     phase using an absorption partioning model (Poa, 1994a, 1994b)
C     for primary semivolatile organic species and their oxidation
C     products (Couvidat et al, JGR, 2012)."
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     kpart: partition coefficient ([m^3.\mu g^-1]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     AERO: aerosol bulk concentration ([\mu g.m^-3]).
C     GAS: gas concentration ([\mu g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2012: Florian Couvidat, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'CONST_A.INC'
      INCLUDE 'param.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'vara.inc'
      INCLUDE 'varp.inc'

      INTEGER nesp_aer

      DOUBLE PRECISION aero(nesp_aer), gas(nesp_aer), kpart(nesp_aer)
      DOUBLE PRECISION dhvap(nesp_aer)

      DOUBLE PRECISION ctot(nesp_pom), caer(nesp_pom)
      DOUBLE PRECISION cgas(nesp_pom), kpart2(nesp_pom)
      DOUBLE PRECISION emw2(nesp_pom), temp

      INTEGER i,j
      DOUBLE PRECISION totom
      DOUBLE PRECISION a, b, c, deter, q, paom,soam

C     Fill concentration vectors
      paom = 0.0
      soam = 0.0
      DO i = 1, nesp_pom
         j = poa_species(i)
         caer(i) = aero(j)
         cgas(i) = gas(j)
         emw2(i) = EMW(j)
         paom = paom + caer(i)
      ENDDO

      IF(paom < 0.1) THEN
         paom = 0.1
      ENDIF

C     Set secondary organic mass quantity in microg/m3
      DO i = NHYDRO + 1, nesp_aec ! only dry AEC species
         j = aec_species(i)
         soam = soam + aero(j)
      ENDDO

      DO i = 1, nesp_pankow
         j = pankow_species(i)
         soam = soam + aero(j)
      ENDDO

C     aero in microg/m3
      DO i = 1, nesp_pom
         j = poa_species(i)
         kpart2(i) = kpart(j)
      ENDDO

      DO i = 1, nesp_pom
         ctot(i) = caer(i) + cgas(i)
      ENDDO

      DO j = 1, NITER_POA
         totom = paom + soam
         paom = 0
         DO i = 1, nesp_pom
            caer(i) = ctot(i) * kpart2(i) * totom /
     &           (1 + kpart2(i) * totom)
            paom = paom + caer(i)
         ENDDO
      ENDDO

      DO j = 1, nesp_pom
         cgas(j) = ctot(j) - caer(j)
      ENDDO

      DO i = 1, nesp_pom
         j = poa_species(i)
         aero(j) = caer(i)
         gas(j) = cgas(i)
      ENDDO

      END
