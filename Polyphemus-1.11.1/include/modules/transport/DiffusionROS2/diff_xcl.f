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


COD Compilation unit : diff_xcl
COD Derivative of unit :  diff_x
COD Dummys:  nx ny nz ts tf dkx dcx dmx rho zc
COD Active IN   dummys:  zc
COD Active OUT  dummys:  zc
COD Dependencies between IN and OUT:
COD zc <--   zc


      SUBROUTINE DIFF_XCL (NX, NY, NZ, TS, TF, DKX, DCX, DMX, RHO, ZC,
     : ZCCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NX
      INTEGER JJ
      INTEGER NY
      INTEGER JK
      INTEGER NZ
      LOGICAL TEST2
      DOUBLE PRECISION TS
      DOUBLE PRECISION TF
      LOGICAL LITERDIFFX
      INTEGER NNN1
      DOUBLE PRECISION DLMATXU(NX)
      DOUBLE PRECISION RHO(NX,NY,NZ)
      DOUBLE PRECISION ZCXCCL(NX)
      DOUBLE PRECISION DLMATXL(NX)
      DOUBLE PRECISION ZKX(NX+1)
      DOUBLE PRECISION DMX(NX)
      DOUBLE PRECISION DKX(NX+1,NY,NZ)
      DOUBLE PRECISION ZCCCL(NX,NY,NZ)
      DOUBLE PRECISION RHOX(NX)
      DOUBLE PRECISION DLMATX(NX)
      DOUBLE PRECISION ZCX(NX)
      DOUBLE PRECISION DCX(NX)
      DOUBLE PRECISION ZC(NX,NY,NZ)
      DOUBLE PRECISION SAVE10(NY,NZ,NX)
      LOGICAL SAVE11(NY,NZ)
      DOUBLE PRECISION SAVE12(NX,NY,NZ)
      INTEGER SAVE13(NY,NZ)
      INTEGER SAVE14(NY,NZ)
      INTEGER SAVE15(NY,NZ)
      INTEGER SAVE16(NZ)
      DOUBLE PRECISION SAVE4(NX,NY,NZ)
      DOUBLE PRECISION SAVE5(NX,NY,NZ)
      DOUBLE PRECISION SAVE6(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE7(NY,NZ,NX)
      DOUBLE PRECISION SAVE8(NY,NZ,NX)
      DOUBLE PRECISION SAVE9(NY,NZ,NX)

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NX
         DLMATX(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLMATXL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLMATXU(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         RHOX(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         ZCX(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX + 1
         ZKX(NNN1) = 0d0
      END DO

      JI = 0
      JJ = 0

C
C Initializations of local variables
C

      DO NNN1 = 1, NX
         ZCXCCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      TEST2 = NX.NE.1
      IF (TEST2) THEN
        LITERDIFFX = .TRUE.
        DO JK = 1, NZ
           SAVE16(JK) = JJ
           DO JJ = 1, NY
              SAVE15(JJ,JK) = JI
              DO JI = 1, NX
                 SAVE4(JI,JJ,JK) = ZCX(JI)
                 ZCX(JI) = ZC(JI,JJ,JK)
                 SAVE5(JI,JJ,JK) = RHOX(JI)
                 RHOX(JI) = RHO(JI,JJ,JK)
              END DO
              SAVE14(JJ,JK) = JI
              DO JI = 1, NX+1
                 SAVE6(JI,JJ,JK) = ZKX(JI)
                 ZKX(JI) = DKX(JI,JJ,JK)
              END DO
              DO NNN1 = 1, NX
                 SAVE7(JJ,JK,NNN1) = DLMATXL(NNN1)
              END DO
              DO NNN1 = 1, NX
                 SAVE8(JJ,JK,NNN1) = DLMATX(NNN1)
              END DO
              DO NNN1 = 1, NX
                 SAVE9(JJ,JK,NNN1) = DLMATXU(NNN1)
              END DO
              DO NNN1 = 1, NX
                 SAVE10(JJ,JK,NNN1) = ZCX(NNN1)
              END DO
              CALL ROSDIFF_X(NX, DCX, DMX, DLMATXL, DLMATX, DLMATXU,
     :         ZKX, RHOX, ZCX, TS, TF, LITERDIFFX)
              SAVE11(JJ,JK) = LITERDIFFX
              LITERDIFFX = .FALSE.
              SAVE13(JJ,JK) = JI
              DO JI = 1, NX
                 SAVE12(JI,JJ,JK) = ZC(JI,JJ,JK)
                 ZC(JI,JJ,JK) = ZCX(JI)
              END DO
           END DO
        END DO
      END IF
C
C Transposed linear forms
C

      IF (TEST2) THEN
        DO JK = NZ, 1, -1
           DO JJ = NY, 1, -1
              DO JI = NX, 1, -1
                 ZC(JI,JJ,JK) = SAVE12(JI,JJ,JK)
                 ZCXCCL(JI) = ZCXCCL(JI)+ZCCCL(JI,JJ,JK)
                 ZCCCL(JI,JJ,JK) = 0d0
              END DO
              JI = SAVE13(JJ,JK)
              LITERDIFFX = SAVE11(JJ,JK)
              DO NNN1 = NX, 1, -1
                 ZCX(NNN1) = SAVE10(JJ,JK,NNN1)
              END DO
              DO NNN1 = NX, 1, -1
                 DLMATXU(NNN1) = SAVE9(JJ,JK,NNN1)
              END DO
              DO NNN1 = NX, 1, -1
                 DLMATX(NNN1) = SAVE8(JJ,JK,NNN1)
              END DO
              DO NNN1 = NX, 1, -1
                 DLMATXL(NNN1) = SAVE7(JJ,JK,NNN1)
              END DO
              CALL ROSDIFF_XCL(NX, DCX, DMX, DLMATXL, DLMATX, DLMATXU,
     :         ZKX, RHOX, ZCX, TS, TF, LITERDIFFX, ZCXCCL)
              DO JI = NX+1, 1, -1
                 ZKX(JI) = SAVE6(JI,JJ,JK)
              END DO
              JI = SAVE14(JJ,JK)
              DO JI = NX, 1, -1
                 RHOX(JI) = SAVE5(JI,JJ,JK)
                 ZCX(JI) = SAVE4(JI,JJ,JK)
                 ZCCCL(JI,JJ,JK) = ZCCCL(JI,JJ,JK)+ZCXCCL(JI)
                 ZCXCCL(JI) = 0d0
              END DO
              JI = SAVE15(JJ,JK)
           END DO
           JJ = SAVE16(JK)
        END DO
      END IF
      END



COD Compilation unit : fexdiff_xcl
COD Derivative of unit :  fexdiff_x
COD Dummys:  nx dcx dmx dlconc dlkx rhox dlr
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlr
COD Dependencies between IN and OUT:
COD dlr <--   dlconc


      SUBROUTINE FEXDIFF_XCL (NX, DCX, DMX, DLCONC, DLKX, RHOX, DLR,
     : DLCONCCCL, DLRCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NX
      INTEGER NNN1
      DOUBLE PRECISION FLUXCCL(NX+1)
      DOUBLE PRECISION FLUX(NX+1)
      DOUBLE PRECISION DLKX(NX+1)
      LOGICAL TEST2(NX+1)
      DOUBLE PRECISION DMX(NX)
      DOUBLE PRECISION DLRCCL(NX)
      LOGICAL TEST5(NX+1)
      DOUBLE PRECISION DLR(NX)
      DOUBLE PRECISION RHOX(NX)
      DOUBLE PRECISION DCX(NX)
      DOUBLE PRECISION DLCONC(NX)
      DOUBLE PRECISION DLCONCCCL(NX)
      DOUBLE PRECISION SAVE6(NX+1)
      DOUBLE PRECISION SAVE7(NX+1)
      DOUBLE PRECISION SAVE8(NX)
      INTEGER SAVE9

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NX+1
         FLUX(NNN1) = 0d0
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NX+1
         FLUXCCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      DO JI = 1, NX+1
         TEST2(JI) = JI.EQ.1
         IF (TEST2(JI)) THEN
           FLUX(JI) = 0.d0
         ELSE
           TEST5(JI) = JI.EQ.NX+1
           IF (TEST5(JI)) THEN
             SAVE6(JI) = FLUX(JI)
             FLUX(JI) = 0.d0
           ELSE
             SAVE7(JI) = FLUX(JI)
             FLUX(JI) =
     :        (DLKX(JI)/DCX(JI))*(DLCONC(JI)/RHOX(JI)-DLCONC(JI-1)/RHOX
     :       (JI-1))
           END IF
         END IF
      END DO
      SAVE9 = JI
      DO JI = 1, NX
         SAVE8(JI) = DLR(JI)
         DLR(JI) = (FLUX(JI+1)-FLUX(JI))/DMX(JI)
      END DO
C
C Transposed linear forms
C

      DO JI = NX, 1, -1
         DLR(JI) = SAVE8(JI)
         FLUXCCL(JI+1) = FLUXCCL(JI+1)+DLRCCL(JI)*(1d0/DMX(JI))
         FLUXCCL(JI) = FLUXCCL(JI)-DLRCCL(JI)*(1d0/DMX(JI))
         DLRCCL(JI) = 0d0
      END DO
      JI = SAVE9
      DO JI = NX+1, 1, -1
         IF (TEST2(JI)) THEN
           FLUXCCL(JI) = 0d0
         ELSE
           IF (TEST5(JI)) THEN
             FLUX(JI) = SAVE6(JI)
             FLUXCCL(JI) = 0d0
           ELSE
             FLUX(JI) = SAVE7(JI)
             DLCONCCCL(JI) =
     :        DLCONCCCL(JI)+FLUXCCL(JI)*((1d0/RHOX(JI))*(DLKX(JI)/DCX
     :       (JI)))
             DLCONCCCL(JI-1) =
     :        DLCONCCCL(JI-1)-FLUXCCL(JI)*((1d0/RHOX(JI-1))*(DLKX(JI)/
     :       DCX(JI)))
             FLUXCCL(JI) = 0d0
           END IF
         END IF
      END DO
      END



COD Unit from the initial code : rosdiff_x



COD Compilation unit : rosdiff_xcl
COD Derivative of unit :  rosdiff_x
COD Dummys:  nx dcx dmx dlmatxl dlmatx dlmatxu dlkx rhox dlconc ts tf literdiffx
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlconc
COD Dependencies between IN and OUT:
COD dlconc <--   dlconc


      SUBROUTINE ROSDIFF_XCL (NX, DCX, DMX, DLMATXL, DLMATX, DLMATXU,
     : DLKX, RHOX, DLCONC, TS, TF, LITERDIFFX, DLCONCCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NX
      DOUBLE PRECISION THRESHOLD
      DOUBLE PRECISION DLSTEP
      DOUBLE PRECISION TS
      DOUBLE PRECISION TF
      DOUBLE PRECISION SD01S
      LOGICAL TEST7
      DOUBLE PRECISION IGAMMA
      LOGICAL LITERDIFFX
      INTEGER NNN1
      LOGICAL TEST26(NX)
      DOUBLE PRECISION DLB2CCL(NX)
      DOUBLE PRECISION DLMATXU(NX)
      DOUBLE PRECISION DLDXU(NX)
      DOUBLE PRECISION DLDX(NX)
      DOUBLE PRECISION DLKX(NX+1)
      DOUBLE PRECISION DLK2CCL(NX)
      DOUBLE PRECISION DLCONCBISCCL(NX)
      DOUBLE PRECISION DLCONC_NEW(NX)
      DOUBLE PRECISION DLMATXL(NX)
      DOUBLE PRECISION DLDXL(NX)
      DOUBLE PRECISION DLB1CCL(NX)
      DOUBLE PRECISION DLCONCBIS(NX)
      DOUBLE PRECISION DMX(NX)
      DOUBLE PRECISION DLK1CCL(NX)
      DOUBLE PRECISION DLB1(NX)
      DOUBLE PRECISION RHOX(NX)
      DOUBLE PRECISION DLMATX(NX)
      DOUBLE PRECISION DLCONC_NEWCCL(NX)
      DOUBLE PRECISION DLB2(NX)
      LOGICAL TEST18(NX)
      DOUBLE PRECISION DLK1(NX)
      DOUBLE PRECISION DCX(NX)
      DOUBLE PRECISION DLK2(NX)
      DOUBLE PRECISION DLCONC(NX)
      DOUBLE PRECISION DLCONCCCL(NX)
      DOUBLE PRECISION SAVE11(NX)
      DOUBLE PRECISION SAVE12(NX)
      DOUBLE PRECISION SAVE13(NX)
      DOUBLE PRECISION SAVE19(NX)
      DOUBLE PRECISION SAVE20(NX)
      DOUBLE PRECISION SAVE22(NX)
      DOUBLE PRECISION SAVE27(NX)
      DOUBLE PRECISION SAVE28(NX)
      INTEGER SAVE29
      INTEGER SAVE30
      INTEGER SAVE31
      DOUBLE PRECISION SAVE9(NX)

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NX
         DLDX(NNN1) = 0d0
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NX
         DLK1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLB2CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLB1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLCONC_NEWCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLCONCBISCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NX
         DLK2CCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      THRESHOLD = 0.d0
      DLSTEP = TF-TS
      SD01S = DSQRT(2.d0)
      IGAMMA = 1.d0+1.d0/SD01S
      CALL FEXDIFF_X(NX, DCX, DMX, DLCONC, DLKX, RHOX, DLB1)
      TEST7 = LITERDIFFX
      IF (TEST7) THEN
        DO NNN1 = 1, NX
           SAVE9(NNN1) = DLDX(NNN1)
        END DO
        CALL JACDDIFFDC_X(NX, DCX, DMX, DLKX, RHOX, DLDXL, DLDX, DLDXU)
        DO JI = 1, NX
           SAVE11(JI) = DLMATX(JI)
           DLMATX(JI) = 1.d0-IGAMMA*DLSTEP*DLDX(JI)
           SAVE12(JI) = DLMATXU(JI)
           DLMATXU(JI) = -IGAMMA*DLSTEP*DLDXU(JI)
           SAVE13(JI) = DLMATXL(JI)
           DLMATXL(JI) = -IGAMMA*DLSTEP*DLDXL(JI)
        END DO
      END IF
      CALL SOLVTRIDIAG_X(NX, DLMATXL, DLMATX, DLMATXU, DLK1, DLB1)
      SAVE31 = JI
      DO JI = 1, NX
         DLCONCBIS(JI) = DLCONC(JI)+DLSTEP*DLK1(JI)
         TEST18(JI) = DLCONCBIS(JI).LE.THRESHOLD
         IF (TEST18(JI)) THEN
           SAVE19(JI) = DLCONCBIS(JI)
           DLCONCBIS(JI) = THRESHOLD
           SAVE20(JI) = DLK1(JI)
           DLK1(JI) = (DLCONCBIS(JI)-DLCONC(JI))/DLSTEP
         END IF
      END DO
      CALL FEXDIFF_X(NX, DCX, DMX, DLCONCBIS, DLKX, RHOX, DLB2)
      SAVE30 = JI
      DO JI = 1, NX
         SAVE22(JI) = DLB2(JI)
         DLB2(JI) = DLB2(JI)-2.d0*DLK1(JI)
      END DO
      CALL SOLVTRIDIAG_X(NX, DLMATXL, DLMATX, DLMATXU, DLK2, DLB2)
      SAVE29 = JI
      DO JI = 1, NX
         DLCONC_NEW(JI) =
     :    DLCONC(JI)+1.5d0*DLSTEP*DLK1(JI)+0.5d0*DLSTEP*DLK2(JI)
         TEST26(JI) = DLCONC_NEW(JI).LE.THRESHOLD
         IF (TEST26(JI)) THEN
           SAVE27(JI) = DLCONC(JI)
           DLCONC(JI) = THRESHOLD
         ELSE
           SAVE28(JI) = DLCONC(JI)
           DLCONC(JI) = DLCONC_NEW(JI)
         END IF
      END DO
C
C Transposed linear forms
C

      DO JI = NX, 1, -1
         IF (TEST26(JI)) THEN
           DLCONC(JI) = SAVE27(JI)
           DLCONCCCL(JI) = 0d0
         ELSE
           DLCONC(JI) = SAVE28(JI)
           DLCONC_NEWCCL(JI) = DLCONC_NEWCCL(JI)+DLCONCCCL(JI)
           DLCONCCCL(JI) = 0d0
         END IF
         DLCONCCCL(JI) = DLCONCCCL(JI)+DLCONC_NEWCCL(JI)
         DLK1CCL(JI) = DLK1CCL(JI)+DLCONC_NEWCCL(JI)*(DLSTEP*1.5d0)
         DLK2CCL(JI) = DLK2CCL(JI)+DLCONC_NEWCCL(JI)*(DLSTEP*0.5d0)
         DLCONC_NEWCCL(JI) = 0d0
      END DO
      JI = SAVE29
      CALL SOLVTRIDIAG_XCL(NX, DLMATXL, DLMATX, DLMATXU, DLK2, DLB2,
     : DLK2CCL, DLB2CCL)
      DO JI = NX, 1, -1
         DLB2(JI) = SAVE22(JI)
         DLK1CCL(JI) = DLK1CCL(JI)-DLB2CCL(JI)*2.d0
      END DO
      JI = SAVE30
      CALL FEXDIFF_XCL(NX, DCX, DMX, DLCONCBIS, DLKX, RHOX, DLB2,
     : DLCONCBISCCL, DLB2CCL)
      DO JI = NX, 1, -1
         IF (TEST18(JI)) THEN
           DLK1(JI) = SAVE20(JI)
           DLCONCBISCCL(JI) = DLCONCBISCCL(JI)+DLK1CCL(JI)*(1d0/DLSTEP)
           DLCONCCCL(JI) = DLCONCCCL(JI)-DLK1CCL(JI)*(1d0/DLSTEP)
           DLK1CCL(JI) = 0d0
           DLCONCBIS(JI) = SAVE19(JI)
           DLCONCBISCCL(JI) = 0d0
         END IF
         DLCONCCCL(JI) = DLCONCCCL(JI)+DLCONCBISCCL(JI)
         DLK1CCL(JI) = DLK1CCL(JI)+DLCONCBISCCL(JI)*DLSTEP
         DLCONCBISCCL(JI) = 0d0
      END DO
      JI = SAVE31
      CALL SOLVTRIDIAG_XCL(NX, DLMATXL, DLMATX, DLMATXU, DLK1, DLB1,
     : DLK1CCL, DLB1CCL)
      IF (TEST7) THEN
        DO JI = NX, 1, -1
           DLMATXL(JI) = SAVE13(JI)
           DLMATXU(JI) = SAVE12(JI)
           DLMATX(JI) = SAVE11(JI)
        END DO
        DO NNN1 = NX, 1, -1
           DLDX(NNN1) = SAVE9(NNN1)
        END DO
      END IF
      CALL FEXDIFF_XCL(NX, DCX, DMX, DLCONC, DLKX, RHOX, DLB1,
     : DLCONCCCL, DLB1CCL)
      END



COD Unit from the initial code : jacddiffdc_x



COD Unit from the initial code : solvtridiag_x



COD Unit from the initial code : fexdiff_x



COD Compilation unit : solvtridiag_xcl
COD Derivative of unit :  solvtridiag_x
COD Dummys:  nx dlmxl dlmx dlmxu dlx dlb
COD Active IN   dummys:  dlb
COD Active OUT  dummys:  dlx
COD Dependencies between IN and OUT:
COD dlx <--   dlb


      SUBROUTINE SOLVTRIDIAG_XCL (NX, DLMXL, DLMX, DLMXU, DLX, DLB,
     : DLXCCL, DLBCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NX
      DOUBLE PRECISION BET
      DOUBLE PRECISION DLX(NX)
      DOUBLE PRECISION DLXCCL(NX)
      DOUBLE PRECISION DLMX(NX)
      LOGICAL TEST11(NX)
      DOUBLE PRECISION DLMXL(NX)
      DOUBLE PRECISION DLB(NX)
      DOUBLE PRECISION GAM(NX)
      LOGICAL TEST4(NX)
      DOUBLE PRECISION DLBCCL(NX)
      LOGICAL TEST8(NX)
      DOUBLE PRECISION DLMXU(NX)
      DOUBLE PRECISION SAVE12(NX)
      INTEGER SAVE13
      DOUBLE PRECISION SAVE2
      DOUBLE PRECISION SAVE6(NX)
      DOUBLE PRECISION SAVE9(NX)

C
C Trajectory
C

      BET = DLMX(1)
      SAVE2 = DLX(1)
      DLX(1) = DLB(1)/BET
      DO JI = 1, NX
         TEST4(JI) = JI.NE.1
         IF (TEST4(JI)) THEN
           GAM(JI) = DLMXU(JI-1)/BET
           SAVE6(JI) = BET
           BET = DLMX(JI)-DLMXL(JI)*GAM(JI)
           TEST8(JI) = BET.EQ.0.d0
           IF (TEST8(JI)) THEN
             WRITE (*, *) 'solvdiff failed'
           END IF
           SAVE9(JI) = DLX(JI)
           DLX(JI) = (DLB(JI)-DLMXL(JI)*DLX(JI-1))/BET
         END IF
      END DO
      SAVE13 = JI
      DO JI = NX, 1, -1
         TEST11(JI) = JI.NE.NX
         IF (TEST11(JI)) THEN
           SAVE12(JI) = DLX(JI)
           DLX(JI) = DLX(JI)-GAM(JI+1)*DLX(JI+1)
         END IF
      END DO
C
C Transposed linear forms
C

      DO JI = 1, NX
         IF (TEST11(JI)) THEN
           DLX(JI) = SAVE12(JI)
           DLXCCL(JI+1) = DLXCCL(JI+1)-DLXCCL(JI)*GAM(JI+1)
         END IF
      END DO
      JI = SAVE13
      DO JI = NX, 1, -1
         IF (TEST4(JI)) THEN
           DLX(JI) = SAVE9(JI)
           DLBCCL(JI) = DLBCCL(JI)+DLXCCL(JI)*(1d0/BET)
           DLXCCL(JI-1) = DLXCCL(JI-1)-DLXCCL(JI)*((1d0/BET)*DLMXL(JI))
           DLXCCL(JI) = 0d0
           BET = SAVE6(JI)
         END IF
      END DO
      DLX(1) = SAVE2
      DLBCCL(1) = DLBCCL(1)+DLXCCL(1)*(1d0/BET)
      DLXCCL(1) = 0d0
      END



