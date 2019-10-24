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


COD Compilation unit : diff_zcl
COD Derivative of unit :  diff_z
COD Dummys:  nx ny nz ts tf rkz cdep cemis rkzf cdepf cemisf dcz dmz rho c
COD Active IN   dummys:  c
COD Active OUT  dummys:  c
COD Dependencies between IN and OUT:
COD c <--   c


      SUBROUTINE DIFF_ZCL (NX, NY, NZ, TS, TF, RKZ, CDEP, CEMIS, RKZF,
     : CDEPF, CEMISF, DCZ, DMZ, RHO, C, CCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NX
      INTEGER JJ
      INTEGER NY
      INTEGER JK
      INTEGER NZ
      DOUBLE PRECISION DCEMIS
      DOUBLE PRECISION DCDEP
      DOUBLE PRECISION TS
      DOUBLE PRECISION TF
      DOUBLE PRECISION DCEMISF
      DOUBLE PRECISION DCDEPF
      INTEGER NNN1
      DOUBLE PRECISION DKZF(NZ+1)
      DOUBLE PRECISION DCZ(NZ)
      DOUBLE PRECISION CCCL(NX,NY,NZ)
      DOUBLE PRECISION CDEP(NX,NY)
      DOUBLE PRECISION RHO(NX,NY,NZ)
      DOUBLE PRECISION CDEPF(NX,NY)
      DOUBLE PRECISION RKZ(NX,NY,NZ+1)
      DOUBLE PRECISION ZCZCCL(NZ)
      DOUBLE PRECISION DMZ(NZ)
      DOUBLE PRECISION C(NX,NY,NZ)
      DOUBLE PRECISION RKZF(NX,NY,NZ+1)
      DOUBLE PRECISION DKZ(NZ+1)
      DOUBLE PRECISION RHOZ(NZ)
      DOUBLE PRECISION CEMIS(NX,NY)
      DOUBLE PRECISION ZCZ(NZ)
      DOUBLE PRECISION CEMISF(NX,NY)
      DOUBLE PRECISION SAVE1(NZ,NX,NY)
      DOUBLE PRECISION SAVE10(NZ,NX,NY)
      INTEGER SAVE11(NX,NY)
      INTEGER SAVE12(NX,NY)
      INTEGER SAVE13(NX,NY)
      INTEGER SAVE14(NY)
      DOUBLE PRECISION SAVE2(NZ,NX,NY)
      DOUBLE PRECISION SAVE3(NZ+1,NX,NY)
      DOUBLE PRECISION SAVE4(NZ+1,NX,NY)
      DOUBLE PRECISION SAVE5(NX,NY)
      DOUBLE PRECISION SAVE6(NX,NY)
      DOUBLE PRECISION SAVE7(NX,NY)
      DOUBLE PRECISION SAVE8(NX,NY)
      DOUBLE PRECISION SAVE9(NX,NY,NZ)


C
C Initializations of uninitialized variables
C

      DCDEP = 0d0
      DCDEPF = 0d0
      DCEMIS = 0d0
      DCEMISF = 0d0
      DO NNN1 = 1, NZ + 1
         DKZ(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ + 1
         DKZF(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         RHOZ(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         ZCZ(NNN1) = 0d0
      END DO

      JI = 0
      JK = 0

C
C Initializations of local variables
C

      DO NNN1 = 1, NZ
         ZCZCCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      DO JJ = 1, NY
         SAVE14(JJ) = JI
         DO JI = 1, NX
            SAVE13(JI,JJ) = JK
            DO JK = 1, NZ
               SAVE1(JK,JI,JJ) = ZCZ(JK)
               ZCZ(JK) = C(JI,JJ,JK)
               SAVE2(JK,JI,JJ) = RHOZ(JK)
               RHOZ(JK) = RHO(JI,JJ,JK)
            END DO
            SAVE12(JI,JJ) = JK
            DO JK = 1, NZ+1
               SAVE3(JK,JI,JJ) = DKZ(JK)
               DKZ(JK) = RKZ(JI,JJ,JK)
               SAVE4(JK,JI,JJ) = DKZF(JK)
               DKZF(JK) = RKZF(JI,JJ,JK)
            END DO
            SAVE5(JI,JJ) = DCDEP
            DCDEP = CDEP(JI,JJ)
            SAVE6(JI,JJ) = DCDEPF
            DCDEPF = CDEPF(JI,JJ)
            SAVE7(JI,JJ) = DCEMIS
            DCEMIS = CEMIS(JI,JJ)
            SAVE8(JI,JJ) = DCEMISF
            DCEMISF = CEMISF(JI,JJ)
            DO NNN1 = 1, NZ
               SAVE9(JI,JJ,NNN1) = ZCZ(NNN1)
            END DO
            CALL ROSDIFF_Z(NZ, ZCZ, TS, TF, DKZ, DKZF, DMZ, DCZ,
     :       DCDEP, DCEMIS, DCDEPF, DCEMISF, RHOZ)
            SAVE11(JI,JJ) = JK
            DO JK = 1, NZ
               SAVE10(JK,JI,JJ) = C(JI,JJ,JK)
               C(JI,JJ,JK) = ZCZ(JK)
            END DO
         END DO
      END DO
C
C Transposed linear forms
C

      DO JJ = NY, 1, -1
         DO JI = NX, 1, -1
            DO JK = NZ, 1, -1
               C(JI,JJ,JK) = SAVE10(JK,JI,JJ)
               ZCZCCL(JK) = ZCZCCL(JK)+CCCL(JI,JJ,JK)
               CCCL(JI,JJ,JK) = 0d0
            END DO
            JK = SAVE11(JI,JJ)
            DO NNN1 = NZ, 1, -1
               ZCZ(NNN1) = SAVE9(JI,JJ,NNN1)
            END DO
            CALL ROSDIFF_ZCL(NZ, ZCZ, TS, TF, DKZ, DKZF, DMZ, DCZ,
     :       DCDEP, DCEMIS, DCDEPF, DCEMISF, RHOZ, ZCZCCL)
            DCEMISF = SAVE8(JI,JJ)
            DCEMIS = SAVE7(JI,JJ)
            DCDEPF = SAVE6(JI,JJ)
            DCDEP = SAVE5(JI,JJ)
            DO JK = NZ+1, 1, -1
               DKZF(JK) = SAVE4(JK,JI,JJ)
               DKZ(JK) = SAVE3(JK,JI,JJ)
            END DO
            JK = SAVE12(JI,JJ)
            DO JK = NZ, 1, -1
               RHOZ(JK) = SAVE2(JK,JI,JJ)
               ZCZ(JK) = SAVE1(JK,JI,JJ)
               CCCL(JI,JJ,JK) = CCCL(JI,JJ,JK)+ZCZCCL(JK)
               ZCZCCL(JK) = 0d0
            END DO
            JK = SAVE13(JI,JJ)
         END DO
         JI = SAVE14(JJ)
      END DO
      END



COD Unit from the initial code : rosdiff_z



COD Compilation unit : fexdiff_zcl
COD Derivative of unit :  fexdiff_z
COD Dummys:  nz dlconc dlr dlkz zdm zdc dcdep dcemis rhoz
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlr
COD Dependencies between IN and OUT:
COD dlr <--   dlconc


      SUBROUTINE FEXDIFF_ZCL (NZ, DLCONC, DLR, DLKZ, ZDM, ZDC, DCDEP,
     : DCEMIS, RHOZ, DLCONCCCL, DLRCCL)

      IMPLICIT NONE
      INTEGER JK
      INTEGER NZ
      DOUBLE PRECISION DCEMIS
      DOUBLE PRECISION DCDEP
      INTEGER NNN1
      DOUBLE PRECISION FLUXCCL(NZ+1)
      DOUBLE PRECISION FLUX(NZ+1)
      LOGICAL TEST2(NZ+1)
      DOUBLE PRECISION DLKZ(NZ+1)
      DOUBLE PRECISION ZDM(NZ)
      DOUBLE PRECISION DLRCCL(NZ)
      LOGICAL TEST5(NZ+1)
      DOUBLE PRECISION DLR(NZ)
      DOUBLE PRECISION ZDC(NZ)
      DOUBLE PRECISION RHOZ(NZ)
      DOUBLE PRECISION DLCONC(NZ)
      DOUBLE PRECISION DLCONCCCL(NZ)
      DOUBLE PRECISION SAVE6(NZ+1)
      DOUBLE PRECISION SAVE7(NZ+1)
      DOUBLE PRECISION SAVE8(NZ)
      INTEGER SAVE9

C
C Initializations of uninitialized variables
C

      DO NNN1 = 1, NZ+1
         FLUX(NNN1) = 0d0
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NZ+1
         FLUXCCL(NNN1) = 0d0
      END DO

C
C Trajectory
C

      DO JK = 1, NZ+1
         TEST2(JK) = JK.EQ.1
         IF (TEST2(JK)) THEN
           FLUX(JK) = DCDEP*DLCONC(JK)-DCEMIS
         ELSE
           TEST5(JK) = JK.EQ.NZ+1
           IF (TEST5(JK)) THEN
             SAVE6(JK) = FLUX(JK)
             FLUX(JK) = 0.d0
           ELSE
             SAVE7(JK) = FLUX(JK)
             FLUX(JK) =
     :        (DLKZ(JK)/ZDC(JK))*(DLCONC(JK)/RHOZ(JK)-DLCONC(JK-1)/RHOZ
     :       (JK-1))
           END IF
         END IF
      END DO
      SAVE9 = JK
      DO JK = 1, NZ
         SAVE8(JK) = DLR(JK)
         DLR(JK) = (FLUX(JK+1)-FLUX(JK))/ZDM(JK)
      END DO
C
C Transposed linear forms
C

      DO JK = NZ, 1, -1
         DLR(JK) = SAVE8(JK)
         FLUXCCL(JK+1) = FLUXCCL(JK+1)+DLRCCL(JK)*(1d0/ZDM(JK))
         FLUXCCL(JK) = FLUXCCL(JK)-DLRCCL(JK)*(1d0/ZDM(JK))
         DLRCCL(JK) = 0d0
      END DO
      JK = SAVE9
      DO JK = NZ+1, 1, -1
         IF (TEST2(JK)) THEN
           DLCONCCCL(JK) = DLCONCCCL(JK)+FLUXCCL(JK)*DCDEP
           FLUXCCL(JK) = 0d0
         ELSE
           IF (TEST5(JK)) THEN
             FLUX(JK) = SAVE6(JK)
             FLUXCCL(JK) = 0d0
           ELSE
             FLUX(JK) = SAVE7(JK)
             DLCONCCCL(JK) =
     :        DLCONCCCL(JK)+FLUXCCL(JK)*((1d0/RHOZ(JK))*(DLKZ(JK)/ZDC
     :       (JK)))
             DLCONCCCL(JK-1) =
     :        DLCONCCCL(JK-1)-FLUXCCL(JK)*((1d0/RHOZ(JK-1))*(DLKZ(JK)/
     :       ZDC(JK)))
             FLUXCCL(JK) = 0d0
           END IF
         END IF
      END DO
      END



COD Unit from the initial code : fexdiff_z



COD Compilation unit : rosdiff_zcl
COD Derivative of unit :  rosdiff_z
COD Dummys:  nz dlconc ts tf dkz dkzf zdm zdc dcdep dcemis dcdepf dcemisf rhoz
COD Active IN   dummys:  dlconc
COD Active OUT  dummys:  dlconc
COD Dependencies between IN and OUT:
COD dlconc <--   dlconc


      SUBROUTINE ROSDIFF_ZCL (NZ, DLCONC, TS, TF, DKZ, DKZF, ZDM, ZDC,
     : DCDEP, DCEMIS, DCDEPF, DCEMISF, RHOZ, DLCONCCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NZ
      DOUBLE PRECISION DCEMIS
      DOUBLE PRECISION THRESHOLD
      DOUBLE PRECISION DCDEP
      DOUBLE PRECISION DLSTEP
      DOUBLE PRECISION TS
      DOUBLE PRECISION TF
      DOUBLE PRECISION DCEMISF
      DOUBLE PRECISION SD01S
      DOUBLE PRECISION DCDEPF
      DOUBLE PRECISION IGAMMA
      INTEGER NNN1
      DOUBLE PRECISION DLB2CCL(NZ)
      DOUBLE PRECISION DKZF(NZ+1)
      DOUBLE PRECISION DLDZ(NZ)
      DOUBLE PRECISION DLDZL(NZ)
      DOUBLE PRECISION DLMATZL(NZ)
      DOUBLE PRECISION DLK2CCL(NZ)
      DOUBLE PRECISION DLCONCBISCCL(NZ)
      DOUBLE PRECISION DLCONC_NEW(NZ)
      DOUBLE PRECISION DLB1CCL(NZ)
      DOUBLE PRECISION ZDM(NZ)
      DOUBLE PRECISION DLCONCBIS(NZ)
      LOGICAL TEST15(NZ)
      DOUBLE PRECISION DLK1CCL(NZ)
      DOUBLE PRECISION ZDC(NZ)
      DOUBLE PRECISION DLB1(NZ)
      DOUBLE PRECISION DLCONC_NEWCCL(NZ)
      DOUBLE PRECISION DKZ(NZ+1)
      DOUBLE PRECISION DLB2(NZ)
      DOUBLE PRECISION DLK1(NZ)
      LOGICAL TEST23(NZ)
      DOUBLE PRECISION RHOZ(NZ)
      DOUBLE PRECISION DLMATZ(NZ)
      DOUBLE PRECISION DLK2(NZ)
      DOUBLE PRECISION DLCONC(NZ)
      DOUBLE PRECISION DLCONCCCL(NZ)
      DOUBLE PRECISION DLMATZU(NZ)
      DOUBLE PRECISION DLDZU(NZ)
      DOUBLE PRECISION SAVE16(NZ)
      DOUBLE PRECISION SAVE17(NZ)
      DOUBLE PRECISION SAVE19(NZ)
      DOUBLE PRECISION SAVE24(NZ)
      DOUBLE PRECISION SAVE25(NZ)
      INTEGER SAVE26
      INTEGER SAVE27
      INTEGER SAVE28

C
C Initializations of local variables
C

      DO NNN1 = 1, NZ
         DLK1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         DLB2CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         DLB1CCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         DLCONC_NEWCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         DLCONCBISCCL(NNN1) = 0d0
      END DO
      DO NNN1 = 1, NZ
         DLK2CCL(NNN1) = 0d0
      END DO
C
C Trajectory
C

      THRESHOLD = 0.d0
      DLSTEP = TF-TS
      SD01S = DSQRT(2.d0)
      IGAMMA = 1.d0+1.d0/SD01S
      CALL FEXDIFF_Z(NZ, DLCONC, DLB1, DKZ, ZDM, ZDC, DCDEP, DCEMIS,
     : RHOZ)
      CALL JACDDIFFDC_Z(NZ, DLDZ, DLDZU, DLDZL, DKZ, ZDM, ZDC, DCDEP,
     : RHOZ)
      DO JI = 1, NZ
         DLMATZ(JI) = 1.d0-IGAMMA*DLSTEP*DLDZ(JI)
         DLMATZU(JI) = -IGAMMA*DLSTEP*DLDZU(JI)
         DLMATZL(JI) = -IGAMMA*DLSTEP*DLDZL(JI)
      END DO
      CALL SOLVTRIDIAG_Z(NZ, DLMATZL, DLMATZ, DLMATZU, DLK1, DLB1)
      SAVE28 = JI
      DO JI = 1, NZ
         DLCONCBIS(JI) = DLCONC(JI)+DLSTEP*DLK1(JI)
         TEST15(JI) = DLCONCBIS(JI).LE.THRESHOLD
         IF (TEST15(JI)) THEN
           SAVE16(JI) = DLCONCBIS(JI)
           DLCONCBIS(JI) = THRESHOLD
           SAVE17(JI) = DLK1(JI)
           DLK1(JI) = (DLCONCBIS(JI)-DLCONC(JI))/DLSTEP
         END IF
      END DO
      CALL FEXDIFF_Z(NZ, DLCONCBIS, DLB2, DKZF, ZDM, ZDC, DCDEPF,
     : DCEMISF, RHOZ)
      SAVE27 = JI
      DO JI = 1, NZ
         SAVE19(JI) = DLB2(JI)
         DLB2(JI) = DLB2(JI)-2.d0*DLK1(JI)
      END DO
      CALL SOLVTRIDIAG_Z(NZ, DLMATZL, DLMATZ, DLMATZU, DLK2, DLB2)
      SAVE26 = JI
      DO JI = 1, NZ
         DLCONC_NEW(JI) =
     :    DLCONC(JI)+1.5d0*DLSTEP*DLK1(JI)+0.5d0*DLSTEP*DLK2(JI)
         TEST23(JI) = DLCONC_NEW(JI).LE.THRESHOLD
         IF (TEST23(JI)) THEN
           SAVE24(JI) = DLCONC(JI)
           DLCONC(JI) = THRESHOLD
         ELSE
           SAVE25(JI) = DLCONC(JI)
           DLCONC(JI) = DLCONC_NEW(JI)
         END IF
      END DO
C
C Transposed linear forms
C

      DO JI = NZ, 1, -1
         IF (TEST23(JI)) THEN
           DLCONC(JI) = SAVE24(JI)
           DLCONCCCL(JI) = 0d0
         ELSE
           DLCONC(JI) = SAVE25(JI)
           DLCONC_NEWCCL(JI) = DLCONC_NEWCCL(JI)+DLCONCCCL(JI)
           DLCONCCCL(JI) = 0d0
         END IF
         DLCONCCCL(JI) = DLCONCCCL(JI)+DLCONC_NEWCCL(JI)
         DLK1CCL(JI) = DLK1CCL(JI)+DLCONC_NEWCCL(JI)*(DLSTEP*1.5d0)
         DLK2CCL(JI) = DLK2CCL(JI)+DLCONC_NEWCCL(JI)*(DLSTEP*0.5d0)
         DLCONC_NEWCCL(JI) = 0d0
      END DO
      JI = SAVE26
      CALL SOLVTRIDIAG_ZCL(NZ, DLMATZL, DLMATZ, DLMATZU, DLK2, DLB2,
     : DLK2CCL, DLB2CCL)
      DO JI = NZ, 1, -1
         DLB2(JI) = SAVE19(JI)
         DLK1CCL(JI) = DLK1CCL(JI)-DLB2CCL(JI)*2.d0
      END DO
      JI = SAVE27
      CALL FEXDIFF_ZCL(NZ, DLCONCBIS, DLB2, DKZF, ZDM, ZDC, DCDEPF,
     : DCEMISF, RHOZ, DLCONCBISCCL, DLB2CCL)
      DO JI = NZ, 1, -1
         IF (TEST15(JI)) THEN
           DLK1(JI) = SAVE17(JI)
           DLCONCBISCCL(JI) = DLCONCBISCCL(JI)+DLK1CCL(JI)*(1d0/DLSTEP)
           DLCONCCCL(JI) = DLCONCCCL(JI)-DLK1CCL(JI)*(1d0/DLSTEP)
           DLK1CCL(JI) = 0d0
           DLCONCBIS(JI) = SAVE16(JI)
           DLCONCBISCCL(JI) = 0d0
         END IF
         DLCONCCCL(JI) = DLCONCCCL(JI)+DLCONCBISCCL(JI)
         DLK1CCL(JI) = DLK1CCL(JI)+DLCONCBISCCL(JI)*DLSTEP
         DLCONCBISCCL(JI) = 0d0
      END DO
      JI = SAVE28
      CALL SOLVTRIDIAG_ZCL(NZ, DLMATZL, DLMATZ, DLMATZU, DLK1, DLB1,
     : DLK1CCL, DLB1CCL)
      CALL FEXDIFF_ZCL(NZ, DLCONC, DLB1, DKZ, ZDM, ZDC, DCDEP, DCEMIS,
     : RHOZ, DLCONCCCL, DLB1CCL)
      END



COD Unit from the initial code : jacddiffdc_z



COD Unit from the initial code : solvtridiag_z



COD Compilation unit : solvtridiag_zcl
COD Derivative of unit :  solvtridiag_z
COD Dummys:  nz dlmzl dlmz dlmzu dlx dlb
COD Active IN   dummys:  dlb
COD Active OUT  dummys:  dlx
COD Dependencies between IN and OUT:
COD dlx <--   dlb


      SUBROUTINE SOLVTRIDIAG_ZCL (NZ, DLMZL, DLMZ, DLMZU, DLX, DLB,
     : DLXCCL, DLBCCL)

      IMPLICIT NONE
      INTEGER JI
      INTEGER NZ
      DOUBLE PRECISION BET
      DOUBLE PRECISION DLX(NZ)
      DOUBLE PRECISION DLMZL(NZ)
      DOUBLE PRECISION DLXCCL(NZ)
      LOGICAL TEST11(NZ)
      DOUBLE PRECISION DLMZ(NZ)
      DOUBLE PRECISION DLB(NZ)
      DOUBLE PRECISION GAM(NZ)
      LOGICAL TEST4(NZ)
      DOUBLE PRECISION DLBCCL(NZ)
      LOGICAL TEST8(NZ)
      DOUBLE PRECISION DLMZU(NZ)
      DOUBLE PRECISION SAVE12(NZ)
      INTEGER SAVE13
      DOUBLE PRECISION SAVE2
      DOUBLE PRECISION SAVE6(NZ)
      DOUBLE PRECISION SAVE9(NZ)

C
C Trajectory
C

      BET = DLMZ(1)
      SAVE2 = DLX(1)
      DLX(1) = DLB(1)/BET
      DO JI = 1, NZ
         TEST4(JI) = JI.NE.1
         IF (TEST4(JI)) THEN
           GAM(JI) = DLMZU(JI-1)/BET
           SAVE6(JI) = BET
           BET = DLMZ(JI)-DLMZL(JI)*GAM(JI)
           TEST8(JI) = BET.EQ.0.d0
           IF (TEST8(JI)) THEN
             WRITE (*, *) 'solvdiff failed'
           END IF
           SAVE9(JI) = DLX(JI)
           DLX(JI) = (DLB(JI)-DLMZL(JI)*DLX(JI-1))/BET
         END IF
      END DO
      SAVE13 = JI
      DO JI = NZ, 1, -1
         TEST11(JI) = JI.NE.NZ
         IF (TEST11(JI)) THEN
           SAVE12(JI) = DLX(JI)
           DLX(JI) = DLX(JI)-GAM(JI+1)*DLX(JI+1)
         END IF
      END DO
C
C Transposed linear forms
C

      DO JI = 1, NZ
         IF (TEST11(JI)) THEN
           DLX(JI) = SAVE12(JI)
           DLXCCL(JI+1) = DLXCCL(JI+1)-DLXCCL(JI)*GAM(JI+1)
         END IF
      END DO
      JI = SAVE13
      DO JI = NZ, 1, -1
         IF (TEST4(JI)) THEN
           DLX(JI) = SAVE9(JI)
           DLBCCL(JI) = DLBCCL(JI)+DLXCCL(JI)*(1d0/BET)
           DLXCCL(JI-1) = DLXCCL(JI-1)-DLXCCL(JI)*((1d0/BET)*DLMZL(JI))
           DLXCCL(JI) = 0d0
           BET = SAVE6(JI)
         END IF
      END DO
      DLX(1) = SAVE2
      DLBCCL(1) = DLBCCL(1)+DLXCCL(1)*(1d0/BET)
      DLXCCL(1) = 0d0
      END



