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


COD Compilation unit : advectioncl


COD Derivative of unit :  advection
COD Dummys:  nx ny nz dmx dmy dmz u v w clesp zclx zcly zclz dt c
COD Active IN   dummys:  c
COD Active OUT  dummys:  c zclx zcly zclz
COD Dependencies between IN and OUT:
COD c <--   c
COD zclx <--   c
COD zcly <--   c
COD zclz <--   c


      SUBROUTINE ADVECTIONCL (NX, NY, NZ, DMX, DMY, DMZ, U, V, W,
     : CLESP, ZCLX, ZCLY, ZCLZ, DT, C, ZCLXCCL, ZCLYCCL, ZCLZCCL,
     : CCCL)

      IMPLICIT NONE
      DOUBLE PRECISION FX_1CCL
      INTEGER NNN3
      INTEGER JI
      INTEGER NX
      INTEGER JJ
      INTEGER NY
      DOUBLE PRECISION DT
      DOUBLE PRECISION NU_P
      DOUBLE PRECISION FY_1CCL
      INTEGER J
      DOUBLE PRECISION FZ_0
      INTEGER JK
      INTEGER NZ
      INTEGER K
      DOUBLE PRECISION FZ_1
      DOUBLE PRECISION C_1CCL
      DOUBLE PRECISION FZ_1CCL
      LOGICAL TEST2
      DOUBLE PRECISION FX_0CCL
      DOUBLE PRECISION FY_0
      INTEGER CLESP
      DOUBLE PRECISION FY_1
      DOUBLE PRECISION WIND
      DOUBLE PRECISION C_3CCL
      DOUBLE PRECISION FY_0CCL
      DOUBLE PRECISION SD01S
      DOUBLE PRECISION FX_0
      DOUBLE PRECISION FX_1
      DOUBLE PRECISION C_0CCL
      DOUBLE PRECISION FZ_0CCL
      DOUBLE PRECISION C_0
      DOUBLE PRECISION C_1
      DOUBLE PRECISION C_2
      INTEGER NNN1
      DOUBLE PRECISION C_2CCL
      DOUBLE PRECISION C_3
      DOUBLE PRECISION NU_M
      INTEGER NNN2
      LOGICAL TEST51(NX+1,NY,NZ)
      LOGICAL TEST89(NX,NY+1,NZ)
      DOUBLE PRECISION FLUXYM(NX,NY+1,NZ)
      DOUBLE PRECISION ZCLXCCL(2,NY,NZ)
      DOUBLE PRECISION U(NX+1,NY,NZ)
      LOGICAL TEST106(NX,NY+1,NZ)
      LOGICAL TEST148(NX,NY,NZ+1)
      LOGICAL TEST169(NX,NY,NZ)
      DOUBLE PRECISION CCCL(NX,NY,NZ)
      DOUBLE PRECISION FLUXZPCCL(NX,NY,NZ+1)
      DOUBLE PRECISION V(NX,NY+1,NZ)
      LOGICAL TEST132(NX,NY,NZ+1)
      LOGICAL TEST94(NX,NY+1,NZ)
      DOUBLE PRECISION W(NX,NY,NZ+1)
      LOGICAL TEST108(NX,NY+1,NZ)
      LOGICAL TEST129(NX,NY,NZ+1)
      LOGICAL TEST11(NX,NZ)
      LOGICAL TEST49(NX+1,NY,NZ)
      DOUBLE PRECISION ZCLZCCL(NX,NY)
      DOUBLE PRECISION FLUXZMCCL(NX,NY,NZ+1)
      DOUBLE PRECISION FLUXYP(NX,NY+1,NZ)
      LOGICAL TEST135(NX,NY,NZ+1)
      LOGICAL TEST156(NX,NY,NZ+1)
      DOUBLE PRECISION ZCLX(2,NY,NZ)
      LOGICAL TEST76(NX+1,NY,NZ)
      DOUBLE PRECISION FLUXYPCCL(NX,NY+1,NZ)
      LOGICAL TEST140(NX,NY,NZ+1)
      DOUBLE PRECISION ZCLY(NX,2,NZ)
      LOGICAL TEST14(NX,NZ)
      LOGICAL TEST98(NX,NY+1,NZ)
      DOUBLE PRECISION DMX(NX)
      LOGICAL TEST116(NX,NY+1,NZ)
      DOUBLE PRECISION ZCLZ(NX,NY)
      LOGICAL TEST4(NY,NZ)
      LOGICAL TEST40(NX+1,NY,NZ)
      LOGICAL TEST57(NX+1,NY,NZ)
      DOUBLE PRECISION FLUXZM(NX,NY,NZ+1)
      DOUBLE PRECISION FLUXYMCCL(NX,NY+1,NZ)
      DOUBLE PRECISION DMY(NY)
      LOGICAL TEST100(NX,NY+1,NZ)
      LOGICAL TEST138(NX,NY,NZ+1)
      LOGICAL TEST37(NX+1,NY,NZ)
      DOUBLE PRECISION FLUXXM(NX+1,NY,NZ)
      DOUBLE PRECISION DMZ(NZ)
      DOUBLE PRECISION ZCLYCCL(NX,2,NZ)
      LOGICAL TEST59(NX+1,NY,NZ)
      DOUBLE PRECISION C(NX,NY,NZ)
      LOGICAL TEST7(NY,NZ)
      LOGICAL TEST18(NX,NY)
      LOGICAL TEST43(NX+1,NY,NZ)
      LOGICAL TEST64(NX+1,NY,NZ)
      DOUBLE PRECISION FLUXZP(NX,NY,NZ+1)
      DOUBLE PRECISION FLUXXPCCL(NX+1,NY,NZ)
      LOGICAL TEST86(NX,NY+1,NZ)
      DOUBLE PRECISION FLUXXP(NX+1,NY,NZ)
      LOGICAL TEST146(NX,NY,NZ+1)
      LOGICAL TEST45(NX+1,NY,NZ)
      LOGICAL TEST70(NX+1,NY,NZ)
      DOUBLE PRECISION FLUXXMCCL(NX+1,NY,NZ)
      LOGICAL TEST92(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE102(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE103(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE104(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE109(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE110(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE111(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE112(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE113(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE114(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE117(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE118(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE119(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE12(NX,NZ)
      DOUBLE PRECISION SAVE120(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE121(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE122(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE123(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE124(NX,NY+1,NZ)
      INTEGER SAVE125(NY+1,NZ)
      INTEGER SAVE126(NZ)
      DOUBLE PRECISION SAVE127(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE130(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE133(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE142(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE143(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE144(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE149(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE15(NX,NZ)
      DOUBLE PRECISION SAVE150(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE151(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE152(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE153(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE154(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE157(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE158(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE159(NX,NY,NZ+1)
      INTEGER SAVE16(NZ)
      DOUBLE PRECISION SAVE160(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE161(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE162(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE163(NX,NY,NZ+1)
      DOUBLE PRECISION SAVE164(NX,NY,NZ+1)
      INTEGER SAVE165(NY,NZ+1)
      INTEGER SAVE166(NZ+1)
      DOUBLE PRECISION SAVE167(NX,NY,NZ)
      DOUBLE PRECISION SAVE170(NX,NY,NZ)
      INTEGER SAVE171(NY,NZ)
      INTEGER SAVE172(NZ)
      INTEGER SAVE173
      INTEGER SAVE174
      INTEGER SAVE175
      DOUBLE PRECISION SAVE19(NX,NY)
      INTEGER SAVE20(NY)
      INTEGER SAVE21
      INTEGER SAVE22
      DOUBLE PRECISION SAVE24(NY,NZ)
      DOUBLE PRECISION SAVE25(NY,NZ)
      INTEGER SAVE26(NZ)
      DOUBLE PRECISION SAVE27(NX,NZ)
      DOUBLE PRECISION SAVE28(NX,NZ)
      INTEGER SAVE29(NZ)
      DOUBLE PRECISION SAVE30(NX,NY)
      INTEGER SAVE31(NY)
      INTEGER SAVE32
      INTEGER SAVE33
      INTEGER SAVE34
      DOUBLE PRECISION SAVE35(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE38(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE41(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE47(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE5(NY,NZ)
      DOUBLE PRECISION SAVE53(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE54(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE55(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE60(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE61(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE62(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE65(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE66(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE67(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE68(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE71(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE72(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE73(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE74(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE77(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE78(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE79(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE8(NY,NZ)
      DOUBLE PRECISION SAVE80(NX+1,NY,NZ)
      DOUBLE PRECISION SAVE81(NX+1,NY,NZ)
      INTEGER SAVE82(NY,NZ)
      INTEGER SAVE83(NZ)
      DOUBLE PRECISION SAVE84(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE87(NX,NY+1,NZ)
      INTEGER SAVE9(NZ)
      DOUBLE PRECISION SAVE90(NX,NY+1,NZ)
      DOUBLE PRECISION SAVE96(NX,NY+1,NZ)

      DOUBLE PRECISION C_RES(NX,NY,NZ)

C
C Initializations of uninitialized variables
C

      C_0 = 0d0
      C_1 = 0d0
      C_2 = 0d0
      C_3 = 0d0
      FX_0 = 0d0
      FX_1 = 0d0
      FY_0 = 0d0
      FY_1 = 0d0
      FZ_0 = 0d0
      FZ_1 = 0d0
      NU_M = 0d0
      NU_P = 0d0
      SD01S = 0d0
      WIND = 0d0

      J = 0
      JI = 0
      JJ = 0

      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX+1
               FLUXXM(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX+1
               FLUXXP(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY+1
            DO NNN3 = 1, NX
               FLUXYM(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY+1
            DO NNN3 = 1, NX
               FLUXYP(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ+1
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX
               FLUXZM(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ+1
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX
               FLUXZP(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO

C
C Initializations of local variables
C

      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY+1
            DO NNN3 = 1, NX
               FLUXYPCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      FX_1CCL = 0d0
      FX_0CCL = 0d0
      C_3CCL = 0d0
      C_2CCL = 0d0
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY+1
            DO NNN3 = 1, NX
               FLUXYMCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      FY_1CCL = 0d0
      C_1CCL = 0d0
      FY_0CCL = 0d0
      C_0CCL = 0d0
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX+1
               FLUXXPCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      FZ_1CCL = 0d0
      DO NNN1 = 1, NZ+1
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX
               FLUXZPCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      FZ_0CCL = 0d0
      DO NNN1 = 1, NZ
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX+1
               FLUXXMCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
      DO NNN1 = 1, NZ+1
         DO NNN2 = 1, NY
            DO NNN3 = 1, NX
               FLUXZMCCL(NNN3,NNN2,NNN1) = 0d0
            END DO
         END DO
      END DO
C
C Trajectory
C

      TEST2 = CLESP.EQ.1
      IF (TEST2) THEN
        DO K = 1, NZ
           SAVE9(K) = J
           DO J = 1, NY
              TEST4(J,K) = U(1,J,K).LE.0.d0
              IF (TEST4(J,K)) THEN
                SAVE5(J,K) = ZCLX(1,J,K)
                ZCLX(1,J,K) = C(1,J,K)
              END IF
              TEST7(J,K) = U(NX+1,J,K).GE.0.d0
              IF (TEST7(J,K)) THEN
                SAVE8(J,K) = ZCLX(2,J,K)
                ZCLX(2,J,K) = C(NX,J,K)
              END IF
           END DO
        END DO
        SAVE22 = K
        DO K = 1, NZ
           SAVE16(K) = J
           DO J = 1, NX
              TEST11(J,K) = V(J,1,K).LE.0.d0
              IF (TEST11(J,K)) THEN
                SAVE12(J,K) = ZCLY(J,1,K)
                ZCLY(J,1,K) = C(J,1,K)
              END IF
              TEST14(J,K) = V(J,NY+1,K).GE.0.d0
              IF (TEST14(J,K)) THEN
                SAVE15(J,K) = ZCLY(J,2,K)
                ZCLY(J,2,K) = C(J,NY,K)
              END IF
           END DO
        END DO
        SAVE21 = K
        DO K = 1, NY
           SAVE20(K) = J
           DO J = 1, NX
              TEST18(J,K) = W(J,K,NZ+1).GE.0.d0
              IF (TEST18(J,K)) THEN
                SAVE19(J,K) = ZCLZ(J,K)
                ZCLZ(J,K) = C(J,K,NZ)
              END IF
           END DO
        END DO
      ELSE
        SAVE34 = K
        DO K = 1, NZ
           SAVE26(K) = J
           DO J = 1, NY
              SAVE24(J,K) = ZCLX(1,J,K)
              ZCLX(1,J,K) = C(1,J,K)
              SAVE25(J,K) = ZCLX(2,J,K)
              ZCLX(2,J,K) = C(NX,J,K)
           END DO
        END DO
        SAVE33 = K
        DO K = 1, NZ
           SAVE29(K) = J
           DO J = 1, NX
              SAVE27(J,K) = ZCLY(J,1,K)
              ZCLY(J,1,K) = C(J,1,K)
              SAVE28(J,K) = ZCLY(J,2,K)
              ZCLY(J,2,K) = C(J,NY,K)
           END DO
        END DO
        SAVE32 = K
        DO K = 1, NY
           SAVE31(K) = J
           DO J = 1, NX
              SAVE30(J,K) = ZCLZ(J,K)
              ZCLZ(J,K) = C(J,K,NZ)
           END DO
        END DO
      END IF
      DO JK = 1, NZ
         SAVE83(JK) = JJ
         DO JJ = 1, NY
            SAVE82(JJ,JK) = JI
            DO JI = 1, NX+1
               SAVE35(JI,JJ,JK) = WIND
               WIND = U(JI,JJ,JK)
               TEST37(JI,JJ,JK) = JI.NE.(NX+1)
               IF (TEST37(JI,JJ,JK)) THEN
                 SAVE38(JI,JJ,JK) = NU_M
                 NU_M = (WIND*DT)/DMX(JI)
               END IF
               TEST40(JI,JJ,JK) = JI.NE.1
               IF (TEST40(JI,JJ,JK)) THEN
                 SAVE41(JI,JJ,JK) = NU_P
                 NU_P = (WIND*DT)/DMX(JI-1)
               END IF
               TEST43(JI,JJ,JK) = JI.EQ.1
               IF (TEST43(JI,JJ,JK)) THEN
                 TEST45(JI,JJ,JK) = WIND.GE.0.d0
                 IF (TEST45(JI,JJ,JK)) THEN
                   FLUXXM(JI,JJ,JK) = NU_M*ZCLX(1,JJ,JK)
                 ELSE
                   SAVE47(JI,JJ,JK) = FLUXXM(JI,JJ,JK)
                   FLUXXM(JI,JJ,JK) = NU_M*C(JI,JJ,JK)
                 END IF
               ELSE
                 TEST49(JI,JJ,JK) = JI.EQ.(NX+1)
                 IF (TEST49(JI,JJ,JK)) THEN
                   TEST51(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST51(JI,JJ,JK)) THEN
                     FLUXXP(JI,JJ,JK) = NU_P*C(JI-1,JJ,JK)
                   ELSE
                     SAVE53(JI,JJ,JK) = FLUXXP(JI,JJ,JK)
                     FLUXXP(JI,JJ,JK) = NU_P*ZCLX(2,JJ,JK)
                   END IF
                 ELSE
                   SAVE54(JI,JJ,JK) = C_1
                   C_1 = C(JI-1,JJ,JK)
                   SAVE55(JI,JJ,JK) = C_2
                   C_2 = C(JI,JJ,JK)
                   TEST57(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST57(JI,JJ,JK)) THEN
                     TEST59(JI,JJ,JK) = JI.EQ.2
                     IF (TEST59(JI,JJ,JK)) THEN
                       SAVE60(JI,JJ,JK) = C_0
                       C_0 = ZCLX(1,JJ,JK)
                     ELSE
                       SAVE61(JI,JJ,JK) = C_0
                       C_0 = C(JI-2,JJ,JK)
                     END IF
                     SAVE62(JI,JJ,JK) = FX_0
                     CALL NUMERICAL_FLUX(NU_M, C_0, C_1, C_2, FX_0)
                     TEST64(JI,JJ,JK) = NU_P.NE.NU_M
                     IF (TEST64(JI,JJ,JK)) THEN
                       SAVE65(JI,JJ,JK) = FX_1
                       CALL NUMERICAL_FLUX(NU_P, C_0, C_1, C_2, FX_1)
                     ELSE
                       SAVE66(JI,JJ,JK) = FX_1
                       FX_1 = FX_0
                     END IF
                     SAVE67(JI,JJ,JK) = FLUXXM(JI,JJ,JK)
                     FLUXXM(JI,JJ,JK) = FX_0
                     SAVE68(JI,JJ,JK) = FLUXXP(JI,JJ,JK)
                     FLUXXP(JI,JJ,JK) = FX_1
                   ELSE
                     TEST70(JI,JJ,JK) = JI.EQ.NX
                     IF (TEST70(JI,JJ,JK)) THEN
                       SAVE71(JI,JJ,JK) = C_3
                       C_3 = ZCLX(2,JJ,JK)
                     ELSE
                       SAVE72(JI,JJ,JK) = C_3
                       C_3 = C(JI+1,JJ,JK)
                     END IF
                     SAVE73(JI,JJ,JK) = SD01S
                     SD01S = -NU_M
                     SAVE74(JI,JJ,JK) = FX_0
                     CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FX_0)
                     TEST76(JI,JJ,JK) = NU_P.NE.NU_M
                     IF (TEST76(JI,JJ,JK)) THEN
                       SAVE77(JI,JJ,JK) = SD01S
                       SD01S = -NU_P
                       SAVE78(JI,JJ,JK) = FX_1
                       CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FX_1)
                     ELSE
                       SAVE79(JI,JJ,JK) = FX_1
                       FX_1 = FX_0
                     END IF
                     SAVE80(JI,JJ,JK) = FLUXXM(JI,JJ,JK)
                     FLUXXM(JI,JJ,JK) = -FX_0
                     SAVE81(JI,JJ,JK) = FLUXXP(JI,JJ,JK)
                     FLUXXP(JI,JJ,JK) = -FX_1
                   END IF
                 END IF
               END IF
            END DO
         END DO
      END DO
      SAVE175 = JK
      DO JK = 1, NZ
         SAVE126(JK) = JJ
         DO JJ = 1, NY+1
            SAVE125(JJ,JK) = JI
            DO JI = 1, NX
               SAVE84(JI,JJ,JK) = WIND
               WIND = V(JI,JJ,JK)
               TEST86(JI,JJ,JK) = JJ.NE.(NY+1)
               IF (TEST86(JI,JJ,JK)) THEN
                 SAVE87(JI,JJ,JK) = NU_M
                 NU_M = (WIND*DT)/DMY(JJ)
               END IF
               TEST89(JI,JJ,JK) = JJ.NE.1
               IF (TEST89(JI,JJ,JK)) THEN
                 SAVE90(JI,JJ,JK) = NU_P
                 NU_P = (WIND*DT)/DMY(JJ-1)
               END IF
               TEST92(JI,JJ,JK) = JJ.EQ.1
               IF (TEST92(JI,JJ,JK)) THEN
                 TEST94(JI,JJ,JK) = WIND.GE.0.d0
                 IF (TEST94(JI,JJ,JK)) THEN
                   FLUXYM(JI,JJ,JK) = NU_M*ZCLY(JI,1,JK)
                 ELSE
                   SAVE96(JI,JJ,JK) = FLUXYM(JI,JJ,JK)
                   FLUXYM(JI,JJ,JK) = NU_M*C(JI,JJ,JK)
                 END IF
               ELSE
                 TEST98(JI,JJ,JK) = JJ.EQ.(NY+1)
                 IF (TEST98(JI,JJ,JK)) THEN
                   TEST100(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST100(JI,JJ,JK)) THEN
                     FLUXYP(JI,JJ,JK) = NU_P*C(JI,JJ-1,JK)
                   ELSE
                     SAVE102(JI,JJ,JK) = FLUXYP(JI,JJ,JK)
                     FLUXYP(JI,JJ,JK) = NU_P*ZCLY(JI,2,JK)
                   END IF
                 ELSE
                   SAVE103(JI,JJ,JK) = C_1
                   C_1 = C(JI,JJ-1,JK)
                   SAVE104(JI,JJ,JK) = C_2
                   C_2 = C(JI,JJ,JK)
                   TEST106(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST106(JI,JJ,JK)) THEN
                     TEST108(JI,JJ,JK) = JJ.EQ.2
                     IF (TEST108(JI,JJ,JK)) THEN
                       SAVE109(JI,JJ,JK) = C_0
                       C_0 = ZCLY(JI,1,JK)
                     ELSE
                       SAVE110(JI,JJ,JK) = C_0
                       C_0 = C(JI,JJ-2,JK)
                     END IF
                     SAVE111(JI,JJ,JK) = FY_0
                     CALL NUMERICAL_FLUX(NU_M, C_0, C_1, C_2, FY_0)
                     SAVE112(JI,JJ,JK) = FY_1
                     CALL NUMERICAL_FLUX(NU_P, C_0, C_1, C_2, FY_1)
                     SAVE113(JI,JJ,JK) = FLUXYM(JI,JJ,JK)
                     FLUXYM(JI,JJ,JK) = FY_0
                     SAVE114(JI,JJ,JK) = FLUXYP(JI,JJ,JK)
                     FLUXYP(JI,JJ,JK) = FY_1
                   ELSE
                     TEST116(JI,JJ,JK) = JJ.EQ.NY
                     IF (TEST116(JI,JJ,JK)) THEN
                       SAVE117(JI,JJ,JK) = C_3
                       C_3 = ZCLY(JI,2,JK)
                     ELSE
                       SAVE118(JI,JJ,JK) = C_3
                       C_3 = C(JI,JJ+1,JK)
                     END IF
                     SAVE119(JI,JJ,JK) = SD01S
                     SD01S = -NU_M
                     SAVE120(JI,JJ,JK) = FY_0
                     CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FY_0)
                     SAVE121(JI,JJ,JK) = SD01S
                     SD01S = -NU_P
                     SAVE122(JI,JJ,JK) = FY_1
                     CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FY_1)
                     SAVE123(JI,JJ,JK) = FLUXYM(JI,JJ,JK)
                     FLUXYM(JI,JJ,JK) = -FY_0
                     SAVE124(JI,JJ,JK) = FLUXYP(JI,JJ,JK)
                     FLUXYP(JI,JJ,JK) = -FY_1
                   END IF
                 END IF
               END IF
            END DO
         END DO
      END DO
      SAVE174 = JK
      DO JK = 1, NZ+1
         SAVE166(JK) = JJ
         DO JJ = 1, NY
            SAVE165(JJ,JK) = JI
            DO JI = 1, NX
               SAVE127(JI,JJ,JK) = WIND
               WIND = W(JI,JJ,JK)
               TEST129(JI,JJ,JK) = JK.NE.(NZ+1)
               IF (TEST129(JI,JJ,JK)) THEN
                 SAVE130(JI,JJ,JK) = NU_M
                 NU_M = (WIND*DT)/DMZ(JK)
               END IF
               TEST132(JI,JJ,JK) = JK.NE.1
               IF (TEST132(JI,JJ,JK)) THEN
                 SAVE133(JI,JJ,JK) = NU_P
                 NU_P = (WIND*DT)/DMZ(JK-1)
               END IF
               TEST135(JI,JJ,JK) = JK.EQ.1
               IF (TEST135(JI,JJ,JK)) THEN
                 FLUXZM(JI,JJ,JK) = 0.d0
               ELSE
                 TEST138(JI,JJ,JK) = JK.EQ.(NZ+1)
                 IF (TEST138(JI,JJ,JK)) THEN
                   TEST140(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST140(JI,JJ,JK)) THEN
                     FLUXZP(JI,JJ,JK) = NU_P*C(JI,JJ,JK-1)
                   ELSE
                     SAVE142(JI,JJ,JK) = FLUXZP(JI,JJ,JK)
                     FLUXZP(JI,JJ,JK) = NU_P*ZCLZ(JI,JJ)
                   END IF
                 ELSE
                   SAVE143(JI,JJ,JK) = C_1
                   C_1 = C(JI,JJ,JK-1)
                   SAVE144(JI,JJ,JK) = C_2
                   C_2 = C(JI,JJ,JK)
                   TEST146(JI,JJ,JK) = WIND.GE.0.d0
                   IF (TEST146(JI,JJ,JK)) THEN
                     TEST148(JI,JJ,JK) = JK.EQ.2
                     IF (TEST148(JI,JJ,JK)) THEN
                       SAVE149(JI,JJ,JK) = C_0
                       C_0 = C_1
                     ELSE
                       SAVE150(JI,JJ,JK) = C_0
                       C_0 = C(JI,JJ,JK-2)
                     END IF
                     SAVE151(JI,JJ,JK) = FZ_0
                     CALL NUMERICAL_FLUX(NU_M, C_0, C_1, C_2, FZ_0)
                     SAVE152(JI,JJ,JK) = FZ_1
                     CALL NUMERICAL_FLUX(NU_P, C_0, C_1, C_2, FZ_1)
                     SAVE153(JI,JJ,JK) = FLUXZM(JI,JJ,JK)
                     FLUXZM(JI,JJ,JK) = FZ_0
                     SAVE154(JI,JJ,JK) = FLUXZP(JI,JJ,JK)
                     FLUXZP(JI,JJ,JK) = FZ_1
                   ELSE
                     TEST156(JI,JJ,JK) = JK.EQ.NZ
                     IF (TEST156(JI,JJ,JK)) THEN
                       SAVE157(JI,JJ,JK) = C_3
                       C_3 = ZCLZ(JI,JJ)
                     ELSE
                       SAVE158(JI,JJ,JK) = C_3
                       C_3 = C(JI,JJ,JK+1)
                     END IF
                     SAVE159(JI,JJ,JK) = SD01S
                     SD01S = -NU_M
                     SAVE160(JI,JJ,JK) = FZ_0
                     CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FZ_0)
                     SAVE161(JI,JJ,JK) = SD01S
                     SD01S = -NU_P
                     SAVE162(JI,JJ,JK) = FZ_1
                     CALL NUMERICAL_FLUX(SD01S, C_3, C_2, C_1, FZ_1)
                     SAVE163(JI,JJ,JK) = FLUXZM(JI,JJ,JK)
                     FLUXZM(JI,JJ,JK) = -FZ_0
                     SAVE164(JI,JJ,JK) = FLUXZP(JI,JJ,JK)
                     FLUXZP(JI,JJ,JK) = -FZ_1
                   END IF
                 END IF
               END IF
            END DO
         END DO
      END DO
      SAVE173 = JK
      DO JK = 1, NZ
         SAVE172(JK) = JJ
         DO JJ = 1, NY
            SAVE171(JJ,JK) = JI
            DO JI = 1, NX
               SAVE167(JI,JJ,JK) = C(JI,JJ,JK)
               C(JI,JJ,JK) =
     :          C(JI,JJ,JK)-FLUXXP(JI+1,JJ,JK)+FLUXXM(JI,JJ,JK)-FLUXYP(
     :         JI,JJ+1,JK)+FLUXYM(JI,JJ,JK)-FLUXZP(JI,JJ,JK+1)+FLUXZM
     :         (JI,JJ,JK)
               TEST169(JI,JJ,JK) = C(JI,JJ,JK).LT.0.d0
               IF (TEST169(JI,JJ,JK)) THEN
                 SAVE170(JI,JJ,JK) = C(JI,JJ,JK)
                 C(JI,JJ,JK) = 0.d0
               END IF
            END DO
         END DO
      END DO

C
C Save forward results.
C
      DO JK = 1, NZ
         DO  JJ = 1, NY
            DO JI = 1, NX
               C_RES(JI, JJ, JK) = C(JI, JJ, JK)
            END DO
         END DO
      END DO

C
C Transposed linear forms
C

      DO JK = NZ, 1, -1
         DO JJ = NY, 1, -1
            DO JI = NX, 1, -1
               IF (TEST169(JI,JJ,JK)) THEN
                 C(JI,JJ,JK) = SAVE170(JI,JJ,JK)
                 CCCL(JI,JJ,JK) = 0d0
               END IF
               C(JI,JJ,JK) = SAVE167(JI,JJ,JK)
               FLUXXPCCL(JI+1,JJ,JK) =
     :          FLUXXPCCL(JI+1,JJ,JK)-CCCL(JI,JJ,JK)
               FLUXXMCCL(JI,JJ,JK) = FLUXXMCCL(JI,JJ,JK)+CCCL(JI,JJ,JK)
               FLUXYPCCL(JI,JJ+1,JK) =
     :          FLUXYPCCL(JI,JJ+1,JK)-CCCL(JI,JJ,JK)
               FLUXYMCCL(JI,JJ,JK) = FLUXYMCCL(JI,JJ,JK)+CCCL(JI,JJ,JK)
               FLUXZPCCL(JI,JJ,JK+1) =
     :          FLUXZPCCL(JI,JJ,JK+1)-CCCL(JI,JJ,JK)
               FLUXZMCCL(JI,JJ,JK) = FLUXZMCCL(JI,JJ,JK)+CCCL(JI,JJ,JK)
            END DO
            JI = SAVE171(JJ,JK)
         END DO
         JJ = SAVE172(JK)
      END DO
      JK = SAVE173
      DO JK = NZ+1, 1, -1
         DO JJ = NY, 1, -1
            DO JI = NX, 1, -1
               IF (TEST135(JI,JJ,JK)) THEN
                 FLUXZMCCL(JI,JJ,JK) = 0d0
               ELSE
                 IF (TEST138(JI,JJ,JK)) THEN
                   IF (TEST140(JI,JJ,JK)) THEN
                     CCCL(JI,JJ,JK-1) =
     :                CCCL(JI,JJ,JK-1)+FLUXZPCCL(JI,JJ,JK)*NU_P
                     FLUXZPCCL(JI,JJ,JK) = 0d0
                   ELSE
                     FLUXZP(JI,JJ,JK) = SAVE142(JI,JJ,JK)
                     ZCLZCCL(JI,JJ) =
     :                ZCLZCCL(JI,JJ)+FLUXZPCCL(JI,JJ,JK)*NU_P
                     FLUXZPCCL(JI,JJ,JK) = 0d0
                   END IF
                 ELSE
                   IF (TEST146(JI,JJ,JK)) THEN
                     FLUXZP(JI,JJ,JK) = SAVE154(JI,JJ,JK)
                     FZ_1CCL = FZ_1CCL+FLUXZPCCL(JI,JJ,JK)
                     FLUXZPCCL(JI,JJ,JK) = 0d0
                     FLUXZM(JI,JJ,JK) = SAVE153(JI,JJ,JK)
                     FZ_0CCL = FZ_0CCL+FLUXZMCCL(JI,JJ,JK)
                     FLUXZMCCL(JI,JJ,JK) = 0d0
                     FZ_1 = SAVE152(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(NU_P, C_0, C_1, C_2, FZ_1,
     :                C_0CCL, C_1CCL, C_2CCL, FZ_1CCL)
                     FZ_0 = SAVE151(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(NU_M, C_0, C_1, C_2, FZ_0,
     :                C_0CCL, C_1CCL, C_2CCL, FZ_0CCL)
                     IF (TEST148(JI,JJ,JK)) THEN
                       C_0 = SAVE149(JI,JJ,JK)
                       C_1CCL = C_1CCL+C_0CCL
                       C_0CCL = 0d0
                     ELSE
                       C_0 = SAVE150(JI,JJ,JK)
                       CCCL(JI,JJ,JK-2) = CCCL(JI,JJ,JK-2)+C_0CCL
                       C_0CCL = 0d0
                     END IF
                   ELSE
                     FLUXZP(JI,JJ,JK) = SAVE164(JI,JJ,JK)
                     FZ_1CCL = FZ_1CCL-FLUXZPCCL(JI,JJ,JK)
                     FLUXZPCCL(JI,JJ,JK) = 0d0
                     FLUXZM(JI,JJ,JK) = SAVE163(JI,JJ,JK)
                     FZ_0CCL = FZ_0CCL-FLUXZMCCL(JI,JJ,JK)
                     FLUXZMCCL(JI,JJ,JK) = 0d0
                     FZ_1 = SAVE162(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1, FZ_1,
     :                C_3CCL, C_2CCL, C_1CCL, FZ_1CCL)
                     SD01S = SAVE161(JI,JJ,JK)
                     FZ_0 = SAVE160(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1, FZ_0,
     :                C_3CCL, C_2CCL, C_1CCL, FZ_0CCL)
                     SD01S = SAVE159(JI,JJ,JK)
                     IF (TEST156(JI,JJ,JK)) THEN
                       C_3 = SAVE157(JI,JJ,JK)
                       ZCLZCCL(JI,JJ) = ZCLZCCL(JI,JJ)+C_3CCL
                       C_3CCL = 0d0
                     ELSE
                       C_3 = SAVE158(JI,JJ,JK)
                       CCCL(JI,JJ,JK+1) = CCCL(JI,JJ,JK+1)+C_3CCL
                       C_3CCL = 0d0
                     END IF
                   END IF
                   C_2 = SAVE144(JI,JJ,JK)
                   CCCL(JI,JJ,JK) = CCCL(JI,JJ,JK)+C_2CCL
                   C_2CCL = 0d0
                   C_1 = SAVE143(JI,JJ,JK)
                   CCCL(JI,JJ,JK-1) = CCCL(JI,JJ,JK-1)+C_1CCL
                   C_1CCL = 0d0
                 END IF
               END IF
               IF (TEST132(JI,JJ,JK)) THEN
                 NU_P = SAVE133(JI,JJ,JK)
               END IF
               IF (TEST129(JI,JJ,JK)) THEN
                 NU_M = SAVE130(JI,JJ,JK)
               END IF
               WIND = SAVE127(JI,JJ,JK)
            END DO
            JI = SAVE165(JJ,JK)
         END DO
         JJ = SAVE166(JK)
      END DO
      JK = SAVE174
      DO JK = NZ, 1, -1
         DO JJ = NY+1, 1, -1
            DO JI = NX, 1, -1
               IF (TEST92(JI,JJ,JK)) THEN
                 IF (TEST94(JI,JJ,JK)) THEN
                   ZCLYCCL(JI,1,JK) =
     :              ZCLYCCL(JI,1,JK)+FLUXYMCCL(JI,JJ,JK)*NU_M
                   FLUXYMCCL(JI,JJ,JK) = 0d0
                 ELSE
                   FLUXYM(JI,JJ,JK) = SAVE96(JI,JJ,JK)
                   CCCL(JI,JJ,JK) =
     :              CCCL(JI,JJ,JK)+FLUXYMCCL(JI,JJ,JK)*NU_M
                   FLUXYMCCL(JI,JJ,JK) = 0d0
                 END IF
               ELSE
                 IF (TEST98(JI,JJ,JK)) THEN
                   IF (TEST100(JI,JJ,JK)) THEN
                     CCCL(JI,JJ-1,JK) =
     :                CCCL(JI,JJ-1,JK)+FLUXYPCCL(JI,JJ,JK)*NU_P
                     FLUXYPCCL(JI,JJ,JK) = 0d0
                   ELSE
                     FLUXYP(JI,JJ,JK) = SAVE102(JI,JJ,JK)
                     ZCLYCCL(JI,2,JK) =
     :                ZCLYCCL(JI,2,JK)+FLUXYPCCL(JI,JJ,JK)*NU_P
                     FLUXYPCCL(JI,JJ,JK) = 0d0
                   END IF
                 ELSE
                   IF (TEST106(JI,JJ,JK)) THEN
                     FLUXYP(JI,JJ,JK) = SAVE114(JI,JJ,JK)
                     FY_1CCL = FY_1CCL+FLUXYPCCL(JI,JJ,JK)
                     FLUXYPCCL(JI,JJ,JK) = 0d0
                     FLUXYM(JI,JJ,JK) = SAVE113(JI,JJ,JK)
                     FY_0CCL = FY_0CCL+FLUXYMCCL(JI,JJ,JK)
                     FLUXYMCCL(JI,JJ,JK) = 0d0
                     FY_1 = SAVE112(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(NU_P, C_0, C_1, C_2, FY_1,
     :                C_0CCL, C_1CCL, C_2CCL, FY_1CCL)
                     FY_0 = SAVE111(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(NU_M, C_0, C_1, C_2, FY_0,
     :                C_0CCL, C_1CCL, C_2CCL, FY_0CCL)
                     IF (TEST108(JI,JJ,JK)) THEN
                       C_0 = SAVE109(JI,JJ,JK)
                       ZCLYCCL(JI,1,JK) = ZCLYCCL(JI,1,JK)+C_0CCL
                       C_0CCL = 0d0
                     ELSE
                       C_0 = SAVE110(JI,JJ,JK)
                       CCCL(JI,JJ-2,JK) = CCCL(JI,JJ-2,JK)+C_0CCL
                       C_0CCL = 0d0
                     END IF
                   ELSE
                     FLUXYP(JI,JJ,JK) = SAVE124(JI,JJ,JK)
                     FY_1CCL = FY_1CCL-FLUXYPCCL(JI,JJ,JK)
                     FLUXYPCCL(JI,JJ,JK) = 0d0
                     FLUXYM(JI,JJ,JK) = SAVE123(JI,JJ,JK)
                     FY_0CCL = FY_0CCL-FLUXYMCCL(JI,JJ,JK)
                     FLUXYMCCL(JI,JJ,JK) = 0d0
                     FY_1 = SAVE122(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1, FY_1,
     :                C_3CCL, C_2CCL, C_1CCL, FY_1CCL)
                     SD01S = SAVE121(JI,JJ,JK)
                     FY_0 = SAVE120(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1, FY_0,
     :                C_3CCL, C_2CCL, C_1CCL, FY_0CCL)
                     SD01S = SAVE119(JI,JJ,JK)
                     IF (TEST116(JI,JJ,JK)) THEN
                       C_3 = SAVE117(JI,JJ,JK)
                       ZCLYCCL(JI,2,JK) = ZCLYCCL(JI,2,JK)+C_3CCL
                       C_3CCL = 0d0
                     ELSE
                       C_3 = SAVE118(JI,JJ,JK)
                       CCCL(JI,JJ+1,JK) = CCCL(JI,JJ+1,JK)+C_3CCL
                       C_3CCL = 0d0
                     END IF
                   END IF
                   C_2 = SAVE104(JI,JJ,JK)
                   CCCL(JI,JJ,JK) = CCCL(JI,JJ,JK)+C_2CCL
                   C_2CCL = 0d0
                   C_1 = SAVE103(JI,JJ,JK)
                   CCCL(JI,JJ-1,JK) = CCCL(JI,JJ-1,JK)+C_1CCL
                   C_1CCL = 0d0
                 END IF
               END IF
               IF (TEST89(JI,JJ,JK)) THEN
                 NU_P = SAVE90(JI,JJ,JK)
               END IF
               IF (TEST86(JI,JJ,JK)) THEN
                 NU_M = SAVE87(JI,JJ,JK)
               END IF
               WIND = SAVE84(JI,JJ,JK)
            END DO
            JI = SAVE125(JJ,JK)
         END DO
         JJ = SAVE126(JK)
      END DO
      JK = SAVE175
      DO JK = NZ, 1, -1
         DO JJ = NY, 1, -1
            DO JI = NX+1, 1, -1
               IF (TEST43(JI,JJ,JK)) THEN
                 IF (TEST45(JI,JJ,JK)) THEN
                   ZCLXCCL(1,JJ,JK) =
     :              ZCLXCCL(1,JJ,JK)+FLUXXMCCL(JI,JJ,JK)*NU_M
                   FLUXXMCCL(JI,JJ,JK) = 0d0
                 ELSE
                   FLUXXM(JI,JJ,JK) = SAVE47(JI,JJ,JK)
                   CCCL(JI,JJ,JK) =
     :              CCCL(JI,JJ,JK)+FLUXXMCCL(JI,JJ,JK)*NU_M
                   FLUXXMCCL(JI,JJ,JK) = 0d0
                 END IF
               ELSE
                 IF (TEST49(JI,JJ,JK)) THEN
                   IF (TEST51(JI,JJ,JK)) THEN
                     CCCL(JI-1,JJ,JK) =
     :                CCCL(JI-1,JJ,JK)+FLUXXPCCL(JI,JJ,JK)*NU_P
                     FLUXXPCCL(JI,JJ,JK) = 0d0
                   ELSE
                     FLUXXP(JI,JJ,JK) = SAVE53(JI,JJ,JK)
                     ZCLXCCL(2,JJ,JK) =
     :                ZCLXCCL(2,JJ,JK)+FLUXXPCCL(JI,JJ,JK)*NU_P
                     FLUXXPCCL(JI,JJ,JK) = 0d0
                   END IF
                 ELSE
                   IF (TEST57(JI,JJ,JK)) THEN
                     FLUXXP(JI,JJ,JK) = SAVE68(JI,JJ,JK)
                     FX_1CCL = FX_1CCL+FLUXXPCCL(JI,JJ,JK)
                     FLUXXPCCL(JI,JJ,JK) = 0d0
                     FLUXXM(JI,JJ,JK) = SAVE67(JI,JJ,JK)
                     FX_0CCL = FX_0CCL+FLUXXMCCL(JI,JJ,JK)
                     FLUXXMCCL(JI,JJ,JK) = 0d0
                     IF (TEST64(JI,JJ,JK)) THEN
                       FX_1 = SAVE65(JI,JJ,JK)
                       CALL NUMERICAL_FLUXCL(NU_P, C_0, C_1, C_2,
     :                  FX_1, C_0CCL, C_1CCL, C_2CCL, FX_1CCL)
                     ELSE
                       FX_1 = SAVE66(JI,JJ,JK)
                       FX_0CCL = FX_0CCL+FX_1CCL
                       FX_1CCL = 0d0
                     END IF
                     FX_0 = SAVE62(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(NU_M, C_0, C_1, C_2, FX_0,
     :                C_0CCL, C_1CCL, C_2CCL, FX_0CCL)
                     IF (TEST59(JI,JJ,JK)) THEN
                       C_0 = SAVE60(JI,JJ,JK)
                       ZCLXCCL(1,JJ,JK) = ZCLXCCL(1,JJ,JK)+C_0CCL
                       C_0CCL = 0d0
                     ELSE
                       C_0 = SAVE61(JI,JJ,JK)
                       CCCL(JI-2,JJ,JK) = CCCL(JI-2,JJ,JK)+C_0CCL
                       C_0CCL = 0d0
                     END IF
                   ELSE
                     FLUXXP(JI,JJ,JK) = SAVE81(JI,JJ,JK)
                     FX_1CCL = FX_1CCL-FLUXXPCCL(JI,JJ,JK)
                     FLUXXPCCL(JI,JJ,JK) = 0d0
                     FLUXXM(JI,JJ,JK) = SAVE80(JI,JJ,JK)
                     FX_0CCL = FX_0CCL-FLUXXMCCL(JI,JJ,JK)
                     FLUXXMCCL(JI,JJ,JK) = 0d0
                     IF (TEST76(JI,JJ,JK)) THEN
                       FX_1 = SAVE78(JI,JJ,JK)
                       CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1,
     :                  FX_1, C_3CCL, C_2CCL, C_1CCL, FX_1CCL)
                       SD01S = SAVE77(JI,JJ,JK)
                     ELSE
                       FX_1 = SAVE79(JI,JJ,JK)
                       FX_0CCL = FX_0CCL+FX_1CCL
                       FX_1CCL = 0d0
                     END IF
                     FX_0 = SAVE74(JI,JJ,JK)
                     CALL NUMERICAL_FLUXCL(SD01S, C_3, C_2, C_1, FX_0,
     :                C_3CCL, C_2CCL, C_1CCL, FX_0CCL)
                     SD01S = SAVE73(JI,JJ,JK)
                     IF (TEST70(JI,JJ,JK)) THEN
                       C_3 = SAVE71(JI,JJ,JK)
                       ZCLXCCL(2,JJ,JK) = ZCLXCCL(2,JJ,JK)+C_3CCL
                       C_3CCL = 0d0
                     ELSE
                       C_3 = SAVE72(JI,JJ,JK)
                       CCCL(JI+1,JJ,JK) = CCCL(JI+1,JJ,JK)+C_3CCL
                       C_3CCL = 0d0
                     END IF
                   END IF
                   C_2 = SAVE55(JI,JJ,JK)
                   CCCL(JI,JJ,JK) = CCCL(JI,JJ,JK)+C_2CCL
                   C_2CCL = 0d0
                   C_1 = SAVE54(JI,JJ,JK)
                   CCCL(JI-1,JJ,JK) = CCCL(JI-1,JJ,JK)+C_1CCL
                   C_1CCL = 0d0
                 END IF
               END IF
               IF (TEST40(JI,JJ,JK)) THEN
                 NU_P = SAVE41(JI,JJ,JK)
               END IF
               IF (TEST37(JI,JJ,JK)) THEN
                 NU_M = SAVE38(JI,JJ,JK)
               END IF
               WIND = SAVE35(JI,JJ,JK)
            END DO
            JI = SAVE82(JJ,JK)
         END DO
         JJ = SAVE83(JK)
      END DO
      IF (TEST2) THEN
        DO K = NY, 1, -1
           DO J = NX, 1, -1
              IF (TEST18(J,K)) THEN
                ZCLZ(J,K) = SAVE19(J,K)
                CCCL(J,K,NZ) = CCCL(J,K,NZ)+ZCLZCCL(J,K)
                ZCLZCCL(J,K) = 0d0
              END IF
           END DO
           J = SAVE20(K)
        END DO
        K = SAVE21
        DO K = NZ, 1, -1
           DO J = NX, 1, -1
              IF (TEST14(J,K)) THEN
                ZCLY(J,2,K) = SAVE15(J,K)
                CCCL(J,NY,K) = CCCL(J,NY,K)+ZCLYCCL(J,2,K)
                ZCLYCCL(J,2,K) = 0d0
              END IF
              IF (TEST11(J,K)) THEN
                ZCLY(J,1,K) = SAVE12(J,K)
                CCCL(J,1,K) = CCCL(J,1,K)+ZCLYCCL(J,1,K)
                ZCLYCCL(J,1,K) = 0d0
              END IF
           END DO
           J = SAVE16(K)
        END DO
        K = SAVE22
        DO K = NZ, 1, -1
           DO J = NY, 1, -1
              IF (TEST7(J,K)) THEN
                ZCLX(2,J,K) = SAVE8(J,K)
                CCCL(NX,J,K) = CCCL(NX,J,K)+ZCLXCCL(2,J,K)
                ZCLXCCL(2,J,K) = 0d0
              END IF
              IF (TEST4(J,K)) THEN
                ZCLX(1,J,K) = SAVE5(J,K)
                CCCL(1,J,K) = CCCL(1,J,K)+ZCLXCCL(1,J,K)
                ZCLXCCL(1,J,K) = 0d0
              END IF
           END DO
           J = SAVE9(K)
        END DO
      ELSE
        DO K = NY, 1, -1
           DO J = NX, 1, -1
              ZCLZ(J,K) = SAVE30(J,K)
              CCCL(J,K,NZ) = CCCL(J,K,NZ)+ZCLZCCL(J,K)
              ZCLZCCL(J,K) = 0d0
           END DO
           J = SAVE31(K)
        END DO
        K = SAVE32
        DO K = NZ, 1, -1
           DO J = NX, 1, -1
              ZCLY(J,2,K) = SAVE28(J,K)
              CCCL(J,NY,K) = CCCL(J,NY,K)+ZCLYCCL(J,2,K)
              ZCLYCCL(J,2,K) = 0d0
              ZCLY(J,1,K) = SAVE27(J,K)
              CCCL(J,1,K) = CCCL(J,1,K)+ZCLYCCL(J,1,K)
              ZCLYCCL(J,1,K) = 0d0
           END DO
           J = SAVE29(K)
        END DO
        K = SAVE33
        DO K = NZ, 1, -1
           DO J = NY, 1, -1
              ZCLX(2,J,K) = SAVE25(J,K)
              CCCL(NX,J,K) = CCCL(NX,J,K)+ZCLXCCL(2,J,K)
              ZCLXCCL(2,J,K) = 0d0
              ZCLX(1,J,K) = SAVE24(J,K)
              CCCL(1,J,K) = CCCL(1,J,K)+ZCLXCCL(1,J,K)
              ZCLXCCL(1,J,K) = 0d0
           END DO
           J = SAVE26(K)
        END DO
        K = SAVE34
      END IF

C
C Restore forward results.
C
      DO JK = 1, NZ
         DO  JJ = 1, NY
            DO JI = 1, NX
               C(JI, JJ, JK) = C_RES(JI, JJ, JK)
            END DO
         END DO
      END DO
      END
