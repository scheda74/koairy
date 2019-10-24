C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
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

      SUBROUTINE EQORG(qext,qgeqi)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the equilibrium between aerosols and
C     gas phase for organics. It provides the gas-phase equilibrium
C     concentrations for organics.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     QEXT : external aerosol concentration ([µg.m-3]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     QGEQI : gas-phase equilibrium concentration ([µg.m-3]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'emw.inc'

      DOUBLE PRECISION qext(NEXT),qgeqi(NEXT)

      INTEGER jesp
      DOUBLE PRECISION xorg(EARO1:EPOA)
      DOUBLE PRECISION worg(EARO1:EPOA)
      DOUBLE PRECISION wtorg


C     Molar fraction of organics

      wtorg=TINYA
      DO jesp=EARO1,EPOA
         worg(jesp)=qext(jesp)/EMW(jesp) ! mol.m-3
         wtorg=wtorg+worg(jesp) ! mol.m-3
      END DO

      DO jesp=EARO1,EPOA
         xorg(jesp)=worg(jesp)/wtorg ! adim
      END DO

C     Partial pressure of ideal mixing

      DO jesp=EARO1,EBiBmP
         qgeqi(jesp)=QSAT(jesp)*xorg(jesp) ! µg.m-3
      END DO

      END

