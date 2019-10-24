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

      SUBROUTINE ISOROPIA_DRV(nesp_aer,
     &     aero, gas, organion, watorg, proton, lwc, rh, temp)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the equilibrium between inorganic aerosols
C     and gas-phase (forward mode), taking in account organic liquid
C     water content and organic ions.
C     It calls ISORROPIA by Nenes et al.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     ORGANION: organic ions ([\mu mol.m^-3]).
C     WATORG: organic liquid water content ([\mu g.m^-3]).
C     RH: relative humidity 0< <1 ([]).
C     TEMP: temperature ([Kelvin]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     AERO: aerosol bulk concentration ([\mu g.m^-3]).
C     GAS: gas concentration ([\mu g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C     PROTON: hydronium ion concentration ([\mu g.m^-3]).
C     LWC: total liquid water content ([\mu g.m^-3]).
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

      INCLUDE 'param.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'vara.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'imw.inc'
c     INCLUDE 'smw.inc'

      INTEGER nesp_aer
      DOUBLE PRECISION aero(nesp_aer), gas(nesp_aer)
      DOUBLE PRECISION organion, watorg, proton
      DOUBLE PRECISION lwc, rh, temp

      DOUBLE PRECISION wi(5),w(5),gas2(3),cntrl(2), other(6)
      DOUBLE PRECISION liquid(12),solid(9)
      DOUBLE PRECISION organion2, watorg2, ionic, gammaH
      INTEGER i,idx             !,j


C     mol neg charge in mol.m-3 */
      organion2 = organion * 1.D-6

C     organic water content converted from
C     microg/m3 (aec output) to kg/m3 (isorropia input)
      watorg2 = watorg * 1.D-9

      cntrl(1) = 0.D0
      cntrl(2) = MTSBL

C     concentration in microg.m-3
      gas(ENa) = 0.D0

      DO i=1,5
         idx = isorropia_species(i)
         wi(i) = aero(idx) + gas(idx)
      ENDDO

C     conversion unit for isorropia needed in mol.m-3
      DO i=1,5
         idx = isorropia_species(i)
         wi(i) = wi(i) / EMW(idx) ! microg.m-3 / microg.mol-1 = mol.m-3
      ENDDO

C     call isorropia fortran routine
      CALL ISOROPIA(wi, rh, temp, cntrl, w, gas2,
     &     liquid, solid, other, organion2, watorg2)

c     DO i=IH,IOH
c     liquid(i) = liquid(i) * IMW(i) ! microg.m-3
c     ENDDO

c     DO i=SNaNO3,SLC
c     j = i - SNaNO3 + 1
c     solid(j) = solid(j) * SMW(i) ! microg.m-3
c     ENDDO

c     DO i=3,nesp_isorropia
c     j = isorropia_species(i)
c     gas2(i-2) = gas2(i-2) * EMW(j) ! microg.m-3
c     ENDDO

C     Isorropia own liquid water content
C     liquid(IH2O) - watorg ! microg.m-3

C     Aqueous phase total liquid water content and pH (proton) concentration
      lwc = liquid(IH2O) * IMW(IH2O) ! microg.m-3
      ionic = other(5)
!     gammaH is in microg.m-3 which is equivalent to micromol.m-3
      gammaH = 10**(-0.511 * (298.0 / temp)**1.5
     &     * sqrt(ionic) / (1 + sqrt(ionic)))
      proton = liquid(IH) * IMW(IH) * gammaH

      gas(ESO4) = 0.D0
      DO i=3,5
         idx = isorropia_species(i)
         gas(idx) = gas2(i-2) * EMW(idx)
      ENDDO

      aero(ESO4) = w(2) * EMW(ESO4)
      DO i=3,5
         idx = isorropia_species(i)
         aero(idx) = (w(i) - gas2(i-2)) * EMW(idx)
      ENDDO

      END
