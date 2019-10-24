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


COD Compilation unit : diff_ycl
COD Derivative of unit :  diff_y
COD Dummys:  nx ny nz ts tf dky dcy dmy rho zc
COD Active IN   dummys:  zc
COD Active OUT  dummys:  zc
COD Dependencies between IN and OUT:
COD zc <--   zc


      SUBROUTINE DIFF_YCL (NX, NY, NZ, TS, TF, DKY, DCY, DMY, RHO, ZC,
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
      LOGICAL LITERDIFFY
      INTEGER NNN1
      DOUBLE PRECISION RHO(NX,NY,NZ)
      DOUBLE PRECISION DLMATYU(NY)
      DOUBLE PRECISION ZKY(NY+1)
      DOUBLE PRECISION DMY(NY)
      DOUBLE PRECISION ZCCCL(NX,NY,NZ)
      DOUBLE PRECISION DKY(NX,NY+1,NZ)
      DOUBLE PRECISION RHOY(NY)
      DOUBLE PRECISION DLMATY(NY)
      DOUBLE PRECISION DLMATYL(NY)
      DOUBLE PRECISION ZCY(NY)
      DOUBLE PRECISION ZCYCCL(NY)
      DOUBLE PRECISION ZC(NX,NY,NZ)
      DOUBLE PRECISION DCY(NY)
      DOUBLE PRECISION SAVE10(NX,NZ,NY)
      LOGICAL SAVE11(NX,NZ)
      DOUBLE PRECISION SAVE12(NY,NX,NZ)
      INTEGER SAVE13(NX,NZ)
      INTEGER SAVE14(NX,NZ)
      INTEGER SAVE15(NX,NZ)
      INTEGER SAVE16(NZ)
      DOUBLE PRECISION SAVE4(NY,NX,NZ)
      DOUBLE PRECISION SAVE5(NY,NX,NZ)
      DOUBLE PRECISION SAVE6(NY+1,NX,NZ)
      DOUBLE PRECISION SAVE7(NX,NZ,NY)
      DOUBLE PRECISION SAVE8(NX,NZ,NY)
      DOUBLE PRECISION SAVE9(NX,NZ,NY)

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NY
         DLMATY(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLMATYL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLMATYU(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         RHOY(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         ZCY(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY + 1
         ZKY(NNN1) = 0d0
      END DO

      JI = 0
      JJ = 0

C
C Initializations of local variables
C

      DO NNN1 = 1, NY
         ZCYCCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      TEST2 = NY.NE.1
      IF (TEST2) THEN
        LITERDIFFY = .TRUE.
        DO JK = 1, NZ
           SAVE16(JK) = JI
           DO JI = 1, NX
              SAVE15(JI,JK) = JJ
              DO JJ = 1, NY
                 SAVE4(JJ,JI,JK) = ZCY(JJ)
                 ZCY(JJ) = ZC(JI,JJ,JK)
                 SAVE5(JJ,JI,JK) = RHOY(JJ)
                 RHOY(JJ) = RHO(JI,JJ,JK)
              END DO
              SAVE14(JI,JK) = JJ
              DO JJ = 1, NY+1
                 SAVE6(JJ,JI,JK) = ZKY(JJ)
                 ZKY(JJ) = DKY(JI,JJ,JK)
              END DO
              DO NNN1 = 1, NY
                 SAVE7(JI,JK,NNN1) = DLMATYL(NNN1)
              END DO
              DO NNN1 = 1, NY
                 SAVE8(JI,JK,NNN1) = DLMATY(NNN1)
              END DO
              DO NNN1 = 1, NY
                 SAVE9(JI,JK,NNN1) = DLMATYU(NNN1)
              END DO
              DO NNN1 = 1, NY
                 SAVE10(JI,JK,NNN1) = ZCY(NNN1)
              END DO
              CALL ROSDIFF_Y(NY, DCY, DMY, DLMATYL, DLMATY, DLMATYU,
     :         ZKY, RHOY, ZCY, TS, TF, LITERDIFFY)
              SAVE11(JI,JK) = LITERDIFFY
              LITERDIFFY = .FALSE.
              SAVE13(JI,JK) = JJ
              DO JJ = 1, NY
                 SAVE12(JJ,JI,JK) = ZC(JI,JJ,JK)
                 ZC(JI,JJ,JK) = ZCY(JJ)
              END DO
           END DO
        END DO
      END IF
C
C Transposed linear forms
C

      IF (TEST2) THEN
        DO JK = NZ, 1, -1
           DO JI = NX, 1, -1
              DO JJ = NY, 1, -1
                 ZC(JI,JJ,JK) = SAVE12(JJ,JI,JK)
                 ZCYCCL(JJ) = ZCYCCL(JJ)+ZCCCL(JI,JJ,JK)
                 ZCCCL(JI,JJ,JK) = 0d0
              END DO
              JJ = SAVE13(JI,JK)
              LITERDIFFY = SAVE11(JI,JK)
              DO NNN1 = NY, 1, -1
                 ZCY(NNN1) = SAVE10(JI,JK,NNN1)
              END DO
              DO NNN1 = NY, 1, -1
                 DLMATYU(NNN1) = SAVE9(JI,JK,NNN1)
              END DO
              DO NNN1 = NY, 1, -1
                 DLMATY(NNN1) = SAVE8(JI,JK,NNN1)
              END DO
              DO NNN1 = NY, 1, -1
                 DLMATYL(NNN1) = SAVE7(JI,JK,NNN1)
              END DO
              CALL ROSDIFF_YCL(NY, DCY, DMY, DLMATYL, DLMATY, DLMATYU,
     :         ZKY, RHOY, ZCY, TS, TF, LITERDIFFY, ZCYCCL)
              DO JJ = NY+1, 1, -1
                 ZKY(JJ) = SAVE6(JJ,JI,JK)
              END DO
              JJ = SAVE14(JI,JK)
              DO JJ = NY, 1, -1
                 RHOY(JJ) = SAVE5(JJ,JI,JK)
                 ZCY(JJ) = SAVE4(JJ,JI,JK)
                 ZCCCL(JI,JJ,JK) = ZCCCL(JI,JJ,JK)+ZCYCCL(JJ)
                 ZCYCCL(JJ) = 0d0
              END DO
              JJ = SAVE15(JI,JK)
           END DO
           JI = SAVE16(JK)
        END DO
      END IF
      END



COD Compilation unit : fexdiff_ycl
COD Derivative of unit :  fexdiff_y
COD Dummys:  ny dcy dmy dlconc dlky rhoy dlr
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlr
COD Dependencies between IN and OUT:
COD dlr <--   dlconc


      SUBROUTINE FEXDIFF_YCL (NY, DCY, DMY, DLCONC, DLKY, RHOY, DLR,
     : DLCONCCCL, DLRCCL)

      IMPLICIT NONE
      INTEGER JJ
      INTEGER NY
      INTEGER NNN1
      DOUBLE PRECISION FLUXCCL(NY+1)
      DOUBLE PRECISION FLUX(NY+1)
      DOUBLE PRECISION DLKY(NY+1)
      LOGICAL TEST2(NY+1)
      DOUBLE PRECISION DLRCCL(NY)
      DOUBLE PRECISION DMY(NY)
      LOGICAL TEST5(NY+1)
      DOUBLE PRECISION DLR(NY)
      DOUBLE PRECISION RHOY(NY)
      DOUBLE PRECISION DLCONC(NY)
      DOUBLE PRECISION DLCONCCCL(NY)
      DOUBLE PRECISION DCY(NY)
      DOUBLE PRECISION SAVE6(NY+1)
      DOUBLE PRECISION SAVE7(NY+1)
      DOUBLE PRECISION SAVE8(NY)
      INTEGER SAVE9

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NY+1
         FLUX(NNN1) = 0d0
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NY+1
         FLUXCCL(NNN1) = 0d0
      END DO

C
C Trajectory
C

      DO JJ = 1, NY+1
         TEST2(JJ) = JJ.EQ.1
         IF (TEST2(JJ)) THEN
           FLUX(JJ) = 0.d0
         ELSE
           TEST5(JJ) = JJ.EQ.NY+1
           IF (TEST5(JJ)) THEN
             SAVE6(JJ) = FLUX(JJ)
             FLUX(JJ) = 0.d0
           ELSE
             SAVE7(JJ) = FLUX(JJ)
             FLUX(JJ) =
     :        (DLKY(JJ)/DCY(JJ))*(DLCONC(JJ)/RHOY(JJ)-DLCONC(JJ-1)/RHOY
     :       (JJ-1))
           END IF
         END IF
      END DO
      SAVE9 = JJ
      DO JJ = 1, NY
         SAVE8(JJ) = DLR(JJ)
         DLR(JJ) = (FLUX(JJ+1)-FLUX(JJ))/DMY(JJ)
      END DO
C
C Transposed linear forms
C

      DO JJ = NY, 1, -1
         DLR(JJ) = SAVE8(JJ)
         FLUXCCL(JJ+1) = FLUXCCL(JJ+1)+DLRCCL(JJ)*(1d0/DMY(JJ))
         FLUXCCL(JJ) = FLUXCCL(JJ)-DLRCCL(JJ)*(1d0/DMY(JJ))
         DLRCCL(JJ) = 0d0
      END DO
      JJ = SAVE9
      DO JJ = NY+1, 1, -1
         IF (TEST2(JJ)) THEN
           FLUXCCL(JJ) = 0d0
         ELSE
           IF (TEST5(JJ)) THEN
             FLUX(JJ) = SAVE6(JJ)
             FLUXCCL(JJ) = 0d0
           ELSE
             FLUX(JJ) = SAVE7(JJ)
             DLCONCCCL(JJ) =
     :        DLCONCCCL(JJ)+FLUXCCL(JJ)*((1d0/RHOY(JJ))*(DLKY(JJ)/DCY
     :       (JJ)))
             DLCONCCCL(JJ-1) =
     :        DLCONCCCL(JJ-1)-FLUXCCL(JJ)*((1d0/RHOY(JJ-1))*(DLKY(JJ)/
     :       DCY(JJ)))
             FLUXCCL(JJ) = 0d0
           END IF
         END IF
      END DO
      END



COD Unit from the initial code : rosdiff_y



COD Compilation unit : rosdiff_ycl
COD Derivative of unit :  rosdiff_y
COD Dummys:  ny dcy dmy dlmatyl dlmaty dlmatyu dlky rhoy dlconc ts tf literdiffy
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlconc
COD Dependencies between IN and OUT:
COD dlconc <--   dlconc


      SUBROUTINE ROSDIFF_YCL (NY, DCY, DMY, DLMATYL, DLMATY, DLMATYU,
     : DLKY, RHOY, DLCONC, TS, TF, LITERDIFFY, DLCONCCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NY
      DOUBLE PRECISION THRESHOLD
      DOUBLE PRECISION DLSTEP
      DOUBLE PRECISION TS
      DOUBLE PRECISION TF
      DOUBLE PRECISION SD01S
      LOGICAL TEST7
      DOUBLE PRECISION IGAMMA
      LOGICAL LITERDIFFY
      INTEGER NNN1
      LOGICAL TEST26(NY)
      DOUBLE PRECISION DLB2CCL(NY)
      DOUBLE PRECISION DLDY(NY)
      DOUBLE PRECISION DLK2CCL(NY)
      DOUBLE PRECISION DLKY(NY+1)
      DOUBLE PRECISION DLCONCBISCCL(NY)
      DOUBLE PRECISION DLCONC_NEW(NY)
      DOUBLE PRECISION DLB1CCL(NY)
      DOUBLE PRECISION DLCONCBIS(NY)
      DOUBLE PRECISION DLMATYU(NY)
      DOUBLE PRECISION DLDYU(NY)
      DOUBLE PRECISION DMY(NY)
      DOUBLE PRECISION DLK1CCL(NY)
      DOUBLE PRECISION DLB1(NY)
      DOUBLE PRECISION DLCONC_NEWCCL(NY)
      DOUBLE PRECISION DLB2(NY)
      LOGICAL TEST18(NY)
      DOUBLE PRECISION RHOY(NY)
      DOUBLE PRECISION DLMATY(NY)
      DOUBLE PRECISION DLK1(NY)
      DOUBLE PRECISION DLMATYL(NY)
      DOUBLE PRECISION DLDYL(NY)
      DOUBLE PRECISION DLK2(NY)
      DOUBLE PRECISION DLCONC(NY)
      DOUBLE PRECISION DLCONCCCL(NY)
      DOUBLE PRECISION DCY(NY)
      DOUBLE PRECISION SAVE11(NY)
      DOUBLE PRECISION SAVE12(NY)
      DOUBLE PRECISION SAVE13(NY)
      DOUBLE PRECISION SAVE19(NY)
      DOUBLE PRECISION SAVE20(NY)
      DOUBLE PRECISION SAVE22(NY)
      DOUBLE PRECISION SAVE27(NY)
      DOUBLE PRECISION SAVE28(NY)
      INTEGER SAVE29
      INTEGER SAVE30
      INTEGER SAVE31
      DOUBLE PRECISION SAVE9(NY)

C
C Initializations of uninitialized variables of double precision
C

      DO NNN1 = 1, NY
         DLDY(NNN1) = 0d0
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NY
         DLK1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLB2CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLB1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLCONC_NEWCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLCONCBISCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NY
         DLK2CCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      THRESHOLD = 0.d0
      DLSTEP = TF-TS
      SD01S = DSQRT(2.d0)
      IGAMMA = 1.d0+1.d0/SD01S
      CALL FEXDIFF_Y(NY, DCY, DMY, DLCONC, DLKY, RHOY, DLB1)
      TEST7 = LITERDIFFY
      IF (TEST7) THEN
        DO NNN1 = 1, NY
           SAVE9(NNN1) = DLDY(NNN1)
        END DO
        CALL JACDDIFFDC_Y(NY, DCY, DMY, DLKY, RHOY, DLDYL, DLDY, DLDYU)
        DO JI = 1, NY
           SAVE11(JI) = DLMATY(JI)
           DLMATY(JI) = 1.d0-IGAMMA*DLSTEP*DLDY(JI)
           SAVE12(JI) = DLMATYU(JI)
           DLMATYU(JI) = -IGAMMA*DLSTEP*DLDYU(JI)
           SAVE13(JI) = DLMATYL(JI)
           DLMATYL(JI) = -IGAMMA*DLSTEP*DLDYL(JI)
        END DO
      END IF
      CALL SOLVTRIDIAG_Y(NY, DLMATYL, DLMATY, DLMATYU, DLK1, DLB1)
      SAVE31 = JI
      DO JI = 1, NY
         DLCONCBIS(JI) = DLCONC(JI)+DLSTEP*DLK1(JI)
         TEST18(JI) = DLCONCBIS(JI).LE.THRESHOLD
         IF (TEST18(JI)) THEN
           SAVE19(JI) = DLCONCBIS(JI)
           DLCONCBIS(JI) = THRESHOLD
           SAVE20(JI) = DLK1(JI)
           DLK1(JI) = (DLCONCBIS(JI)-DLCONC(JI))/DLSTEP
         END IF
      END DO
      CALL FEXDIFF_Y(NY, DCY, DMY, DLCONCBIS, DLKY, RHOY, DLB2)
      SAVE30 = JI
      DO JI = 1, NY
         SAVE22(JI) = DLB2(JI)
         DLB2(JI) = DLB2(JI)-2.d0*DLK1(JI)
      END DO
      CALL SOLVTRIDIAG_Y(NY, DLMATYL, DLMATY, DLMATYU, DLK2, DLB2)
      SAVE29 = JI
      DO JI = 1, NY
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

      DO JI = NY, 1, -1
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
      CALL SOLVTRIDIAG_YCL(NY, DLMATYL, DLMATY, DLMATYU, DLK2, DLB2,
     : DLK2CCL, DLB2CCL)
      DO JI = NY, 1, -1
         DLB2(JI) = SAVE22(JI)
         DLK1CCL(JI) = DLK1CCL(JI)-DLB2CCL(JI)*2.d0
      END DO
      JI = SAVE30
      CALL FEXDIFF_YCL(NY, DCY, DMY, DLCONCBIS, DLKY, RHOY, DLB2,
     : DLCONCBISCCL, DLB2CCL)
      DO JI = NY, 1, -1
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
      CALL SOLVTRIDIAG_YCL(NY, DLMATYL, DLMATY, DLMATYU, DLK1, DLB1,
     : DLK1CCL, DLB1CCL)
      IF (TEST7) THEN
        DO JI = NY, 1, -1
           DLMATYL(JI) = SAVE13(JI)
           DLMATYU(JI) = SAVE12(JI)
           DLMATY(JI) = SAVE11(JI)
        END DO
        DO NNN1 = NY, 1, -1
           DLDY(NNN1) = SAVE9(NNN1)
        END DO
      END IF
      CALL FEXDIFF_YCL(NY, DCY, DMY, DLCONC, DLKY, RHOY, DLB1,
     : DLCONCCCL, DLB1CCL)
      END



COD Unit from the initial code : jacddiffdc_y



COD Unit from the initial code : solvtridiag_y



COD Unit from the initial code : fexdiff_y



COD Compilation unit : solvtridiag_ycl
COD Derivative of unit :  solvtridiag_y
COD Dummys:  ny dlmyl dlmy dlmyu dlx dlb
COD Active IN   dummys:  dlb
COD Active OUT  dummys:  dlx
COD Dependencies between IN and OUT:
COD dlx <--   dlb


      SUBROUTINE SOLVTRIDIAG_YCL (NY, DLMYL, DLMY, DLMYU, DLX, DLB,
     : DLXCCL, DLBCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NY
      DOUBLE PRECISION BET
      DOUBLE PRECISION DLX(NY)
      DOUBLE PRECISION DLXCCL(NY)
      LOGICAL TEST11(NY)
      DOUBLE PRECISION DLMY(NY)
      DOUBLE PRECISION DLB(NY)
      DOUBLE PRECISION GAM(NY)
      DOUBLE PRECISION DLMYU(NY)
      LOGICAL TEST4(NY)
      DOUBLE PRECISION DLBCCL(NY)
      DOUBLE PRECISION DLMYL(NY)
      LOGICAL TEST8(NY)
      DOUBLE PRECISION SAVE12(NY)
      INTEGER SAVE13
      DOUBLE PRECISION SAVE2
      DOUBLE PRECISION SAVE6(NY)
      DOUBLE PRECISION SAVE9(NY)

C
C Trajectory
C

      BET = DLMY(1)
      SAVE2 = DLX(1)
      DLX(1) = DLB(1)/BET
      DO JI = 1, NY
         TEST4(JI) = JI.NE.1
         IF (TEST4(JI)) THEN
           GAM(JI) = DLMYU(JI-1)/BET
           SAVE6(JI) = BET
           BET = DLMY(JI)-DLMYL(JI)*GAM(JI)
           TEST8(JI) = BET.EQ.0.d0
           IF (TEST8(JI)) THEN
             WRITE (*, *) 'solvdiff failed'
           END IF
           SAVE9(JI) = DLX(JI)
           DLX(JI) = (DLB(JI)-DLMYL(JI)*DLX(JI-1))/BET
         END IF
      END DO
      SAVE13 = JI
      DO JI = NY, 1, -1
         TEST11(JI) = JI.NE.NY
         IF (TEST11(JI)) THEN
           SAVE12(JI) = DLX(JI)
           DLX(JI) = DLX(JI)-GAM(JI+1)*DLX(JI+1)
         END IF
      END DO
C
C Transposed linear forms
C

      DO JI = 1, NY
         IF (TEST11(JI)) THEN
           DLX(JI) = SAVE12(JI)
           DLXCCL(JI+1) = DLXCCL(JI+1)-DLXCCL(JI)*GAM(JI+1)
         END IF
      END DO
      JI = SAVE13
      DO JI = NY, 1, -1
         IF (TEST4(JI)) THEN
           DLX(JI) = SAVE9(JI)
           DLBCCL(JI) = DLBCCL(JI)+DLXCCL(JI)*(1d0/BET)
           DLXCCL(JI-1) = DLXCCL(JI-1)-DLXCCL(JI)*((1d0/BET)*DLMYL(JI))
           DLXCCL(JI) = 0d0
           BET = SAVE6(JI)
         END IF
      END DO
      DLX(1) = SAVE2
      DLBCCL(1) = DLBCCL(1)+DLXCCL(1)*(1d0/BET)
      DLXCCL(1) = 0d0
      END



