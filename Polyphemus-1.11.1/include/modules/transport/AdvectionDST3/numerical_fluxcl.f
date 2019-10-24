C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
C
C     This file is part of the air quality modeling system Polyphemus.
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


COD Unit from the initial code : numerical_flux



COD Compilation unit : numerical_fluxcl
COD Derivative of unit :  numerical_flux
COD Dummys:  nu c_0 c_1 c_2 flux
COD Active IN   dummys:  c_0 c_1 c_2 flux
COD Active OUT  dummys:  flux
COD Dependencies between IN and OUT:
COD flux <--   c_0 c_1 c_2 flux


      SUBROUTINE NUMERICAL_FLUXCL (NU, C_0, C_1, C_2, FLUX, C_0CCL,
     : C_1CCL, C_2CCL, FLUXCCL)

      IMPLICIT NONE
      DOUBLE PRECISION PHICCL
      DOUBLE PRECISION M1
      DOUBLE PRECISION FLUXCCL
      DOUBLE PRECISION FLUX
      DOUBLE PRECISION SD01SCCL
      DOUBLE PRECISION M2
      LOGICAL TEST11
      DOUBLE PRECISION RD02R
      DOUBLE PRECISION DIFF_C2CCL
      DOUBLE PRECISION PHI
      DOUBLE PRECISION C_1CCL
      DOUBLE PRECISION M2CCL
      LOGICAL TEST14
      LOGICAL TEST20
      DOUBLE PRECISION SD01S
      DOUBLE PRECISION DIFF_C1CCL
      LOGICAL TEST6
      DOUBLE PRECISION DIFF_C1
      DOUBLE PRECISION C_0CCL
      DOUBLE PRECISION C_0
      DOUBLE PRECISION DIFF_C2
      DOUBLE PRECISION NU
      DOUBLE PRECISION M1CCL
      DOUBLE PRECISION C_1
      DOUBLE PRECISION C_2
      DOUBLE PRECISION SGNDIFF_C2
      LOGICAL TEST25
      DOUBLE PRECISION C_2CCL
      DOUBLE PRECISION SAVE1
      DOUBLE PRECISION SAVE15
      DOUBLE PRECISION SAVE16
      DOUBLE PRECISION SAVE22
      DOUBLE PRECISION SAVE23
      DOUBLE PRECISION SAVE26
      DOUBLE PRECISION SAVE27
      DOUBLE PRECISION SAVE28
      DOUBLE PRECISION SAVE8
      DOUBLE PRECISION SAVE9

C
C Initializations of local variables
C

      DIFF_C2CCL = 0d0
      DIFF_C1CCL = 0d0
      M2CCL = 0d0
      M1CCL = 0d0
      SD01SCCL = 0d0
      PHICCL = 0d0
C
C Trajectory
C

      SAVE1 = FLUX
      FLUX = NU*C_1
      DIFF_C1 = C_1-C_0
      DIFF_C2 = C_2-C_1
      RD02R = 1.d0
      TEST6 = DIFF_C2.GE.0d0
      IF (TEST6) THEN
        SGNDIFF_C2 = RD02R
      ELSE
        SAVE8 = SGNDIFF_C2
        SGNDIFF_C2 = -RD02R
      END IF
      SAVE9 = DIFF_C1
      DIFF_C1 = DIFF_C1*SGNDIFF_C2
      TEST11 = DIFF_C1.GT.0.d0
      IF (TEST11) THEN
        SD01S = DIFF_C2
        TEST14 = SD01S.GE.0d0
        IF (TEST14) THEN
          SAVE15 = DIFF_C2
          DIFF_C2 = SD01S
        ELSE
          SAVE16 = DIFF_C2
          DIFF_C2 = -SD01S
        END IF
        M1 = ((1.d0-NU)/6.d0)*((2.d0-NU)*DIFF_C2+(1.d0+NU)*DIFF_C1)
        M2 = (1.d0-NU)*DIFF_C1
        TEST20 = DIFF_C2.LE.M1
        IF (TEST20) THEN
          PHI = DIFF_C2
        ELSE
          SAVE22 = PHI
          PHI = M1
        END IF
        SAVE23 = SD01S
        SD01S = PHI*NU
        TEST25 = SD01S.LE.M2
        IF (TEST25) THEN
          SAVE26 = PHI
          PHI = SD01S
        ELSE
          SAVE27 = PHI
          PHI = M2
        END IF
        SAVE28 = FLUX
        FLUX = FLUX+PHI*SGNDIFF_C2
      END IF
C
C Transposed linear forms
C

      IF (TEST11) THEN
        FLUX = SAVE28
        PHICCL = PHICCL+FLUXCCL*SGNDIFF_C2
        IF (TEST25) THEN
          PHI = SAVE26
          SD01SCCL = SD01SCCL+PHICCL
          PHICCL = 0d0
        ELSE
          PHI = SAVE27
          M2CCL = M2CCL+PHICCL
          PHICCL = 0d0
        END IF
        SD01S = SAVE23
        PHICCL = PHICCL+SD01SCCL*NU
        SD01SCCL = 0d0
        IF (TEST20) THEN
          DIFF_C2CCL = DIFF_C2CCL+PHICCL
          PHICCL = 0d0
        ELSE
          PHI = SAVE22
          M1CCL = M1CCL+PHICCL
          PHICCL = 0d0
        END IF
        DIFF_C1CCL = DIFF_C1CCL+M2CCL*(1.d0-NU)
        M2CCL = 0d0
        DIFF_C2CCL = DIFF_C2CCL+M1CCL*((2.d0-NU)*((1.d0-NU)/6.d0))
        DIFF_C1CCL = DIFF_C1CCL+M1CCL*((1.d0+NU)*((1.d0-NU)/6.d0))
        M1CCL = 0d0
        IF (TEST14) THEN
          DIFF_C2 = SAVE15
          SD01SCCL = SD01SCCL+DIFF_C2CCL
          DIFF_C2CCL = 0d0
        ELSE
          DIFF_C2 = SAVE16
          SD01SCCL = SD01SCCL-DIFF_C2CCL
          DIFF_C2CCL = 0d0
        END IF
        DIFF_C2CCL = DIFF_C2CCL+SD01SCCL
        SD01SCCL = 0d0
      END IF
      DIFF_C1 = SAVE9
      DIFF_C1CCL = DIFF_C1CCL*SGNDIFF_C2
      IF (TEST6) THEN
      ELSE
        SGNDIFF_C2 = SAVE8
      END IF
      C_2CCL = C_2CCL+DIFF_C2CCL
      C_1CCL = C_1CCL-DIFF_C2CCL
      DIFF_C2CCL = 0d0
      C_1CCL = C_1CCL+DIFF_C1CCL
      C_0CCL = C_0CCL-DIFF_C1CCL
      DIFF_C1CCL = 0d0
      FLUX = SAVE1
      C_1CCL = C_1CCL+FLUXCCL*NU
      FLUXCCL = 0d0
      END
