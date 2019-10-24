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

      SUBROUTINE PANKOW_DRV(nesp_aer,aero, gas, kpart)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the equilibrium between gas and particle
C     phase using an absorption partioning model (Pankow, 1994a, 1994b)
C     for organic species which are not managed by AEC model (Pun et al 2001).
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
C     2007: Edouard Debry, CEREA.
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

      DOUBLE PRECISION ctot(nesp_pankow), caer(nesp_pankow)
      DOUBLE PRECISION cgas(nesp_pankow), kpart2(nesp_pankow)
      DOUBLE PRECISION emw2(nesp_pankow)

      INTEGER i,j
      DOUBLE PRECISION totmol, totmol2
      DOUBLE PRECISION a, b, c, deter, q, paom

C     Fill concentration vectors
      DO i = 1,nesp_pankow
         j = pankow_species(i)
         caer(i) = aero(j)
         cgas(i) = gas(j)
         emw2(i) = EMW(j)
      ENDDO

C     Set primary organic molar quantity in mol/m3
      paom = 0.0
      DO i = 1, nesp_pom
         j = poa_species(i)
         paom = paom + aero(j) / EMW(j)
      ENDDO
      DO i = NHYDRO + 1, nesp_aec ! only dry AEC species
         j = aec_species(i)
         paom = paom + aero(j) / EMW(j)
      ENDDO

C     aero in microg/m3, EMW in microg/mol
      DO i = 1,nesp_pankow
         j = pankow_species(i)
         kpart2(i) = KPART(j)
      ENDDO

      DO i = 1,nesp_pankow
         ctot(i) = caer(i) + cgas(i)
      ENDDO

      DO j=1,NITER_PKW
         totmol = paom

         DO i=1,nesp_pankow
            totmol = totmol + caer(i) / emw2(i)
         ENDDO

         DO i=1,nesp_pankow
            totmol2 = totmol - caer(i) / emw2(i)

            a = 1.D0 / emw2(i)
            b = (1.D0 / kpart2(i) - ctot(i)) / emw2(i) + totmol2
            c = - ctot(i) * totmol2

            deter = b * b - 4.D0 * a * c
            IF (deter.LT.0.d0) STOP 'pankow_drv.f: deter < 0'

            q= - 0.5D0 * ( b + DSIGN(1.D0,b) * DSQRT(deter))
            caer(i) = DMAX1(q / a, c / q)
         ENDDO
      ENDDO

      DO j = 1,nesp_pankow
         cgas(j) = ctot(j) - caer(j)
      ENDDO

      DO i = 1,nesp_pankow
         j = pankow_species(i)
         aero(j) = caer(i)
         gas(j) = cgas(i)
      ENDDO

      END
