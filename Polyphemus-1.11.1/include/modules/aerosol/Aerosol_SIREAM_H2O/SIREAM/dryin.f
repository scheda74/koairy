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

      SUBROUTINE DRYIN(temp,qinti,nsize,qgbki,aai,ckvi,qgeqi,kercdi)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the aerosol surface gas-phase
C     concentration (through equilibrium) for dry aerosols.
C
C     The algorithms are detailed in Chapter 10 (section 10.1.6) of
C     the PhD work of Edouard Debry.
C     See also the reference:
C     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
C     Aerosol Science and Technology 32, 482:502, 2000.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     TEMP  : local temperature ([Kelvin]).
C     QINTI : internal  solid aerosol concentration ([�g.m-3]).
C     nsize : size of all vectors in argument following below.
C     QGBKI : bulk gas concentration ([�g.m-3]).
C     AAI   : condensation/evaporation kernel coefficient([m3.s-1]).
C     CKVI  : Kelvin effect correction ([]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     QGEQI  : gas-phase concentration at the aerosol surface ([�g.m-3].
C     KERCDI : condensation/evaporation kernel ([�g.s-1]).
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
C     2006/1/24: remove isokeq.inc and put commons directly from
C                isrpia.inc in the core of routine, because
C                we can no longer afford to put isokeq.inc in
C                ISORROPIA routines (Edouard Debry, CEREA)
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
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'vara.inc'
      INCLUDE 'varp.inc'

C ISORROPIA commons needed by this routine,
C directly taken from isrpia.inc.
      DOUBLE PRECISION XK1,XK2,XK3,XK4,XK5,XK6,XK7,XK8,
     &                 XK9,XK10,XK11,XK12,XK13,XK14,
     &                 XKW,XK21,XK22,XK31,XK32,XK41,XK42
      COMMON /EQUK/ XK1,XK2,XK3,XK4,XK5,XK6,XK7,XK8,
     &              XK9,XK10,XK11,XK12,XK13,XK14,
     &              XKW,XK21,XK22,XK31,XK32,XK41,XK42
!$OMP THREADPRIVATE(/EQUK/)


      INTEGER nsize
      DOUBLE PRECISION qinti(NINTIS),ckvi(nsize)
      DOUBLE PRECISION aai(nsize),qgbki(nsize)
      DOUBLE PRECISION qgeqi(nsize),kercdi(nsize)
      DOUBLE PRECISION temp

      LOGICAL leq1,leq2,leq3,lr47,lr56
      INTEGER icase,jesp,jj
      DOUBLE PRECISION rgas1,maa(nsize)
      DOUBLE PRECISION rk1,rk2,rk3
      DOUBLE PRECISION sat1,sat2,sat3a,sat3b
      DOUBLE PRECISION mkercd(nsize),msat
      DOUBLE PRECISION cfa,cfb,cfc

C     Initialization:
C     1 stands for nh4no3 equilibrium
C     2 stands for nh4cl equilibrium
C     3 stands for nacl/nano3 equilibrium
C     47 for nacl and nh4cl reactions
C     56 for nano3 and nh4no3 reactions

      leq1=.false.
      leq2=.false.
      leq3=.false.
      lr47=.false.
      lr56=.false.

      rgas1=ATM/(RGAS*temp)

      rk1= XK10*rgas1*rgas1
     &     *EMW(ENH3)*EMW(ENO3)
     &     *ckvi(ENH3)*ckvi(ENO3)
                                ! rk1 in (�g.m-3)2

      rk2=XK6*rgas1*rgas1
     &     *EMW(ENH3)*EMW(ECl)
     &     *ckvi(ENH3)*ckvi(ECl)
                                ! rk2 in (�g.m-3)2

      rk3=XK4*XK8/(XK3*XK9) ! adim

      sat1=qgbki(ENH3)*qgbki(ENO3) ! (�g.m-3)2
      sat2=qgbki(ENH3)*qgbki(ECl) ! (�g.m-3)2

      sat3a=qgbki(ECl)          ! �g.m-3
      sat3b=rk3*qgbki(ENO3)     ! �g.m-3

      DO jj=2,nesp_isorropia
         jesp=isorropia_species(jj)
                                ! maa(*) in m3.mol.s-1.�g-1
         maa(jesp)= aai(jesp)   ! m3.s-1
     &        /EMW(jesp)        ! �g.mol-1

                                ! mkercd(*) in mol.s-1
         mkercd(jesp)= kercdi(jesp) ! �g.s-1
     &        /EMW(jesp)        ! �g.mol-1
      END DO

      msat=2.D0*mkercd(ESO4)-mkercd(ENH3) ! mol.s-1

C     Determine which reaction is active

      IF (qinti(SNH4NO3).GT.0.D0.OR.
     &     sat1.GT.rk1) THEN

         leq1=.true.
      ENDIF

      IF (qinti(SNH4Cl).GT.0.D0.OR.
     &     sat2.GT.rk2) THEN

         leq2=.true.
      ENDIF

      IF (qinti(SNaNO3).GT.0.D0.AND.
     &     qinti(SNaCl).GT.0.D0) THEN

         leq3=.true.
      ENDIF

      IF (qinti(SNaNO3).GT.0.D0.AND.
     &     sat3a.GT.sat3b) THEN

         leq3=.true.
      ENDIF

      IF (qinti(SNaCl).GT.0.D0.AND.
     &     sat3a.LT.sat3b) THEN

         leq3=.true.
      ENDIF

      IF (qgbki(ESO4).GT.0.D0) THEN
         IF (qinti(SNaCl).GT.0.D0.OR.
     &        qinti(SNH4Cl).GT.0.D0) THEN

            lr47=.true.
         ENDIF

         IF (qinti(SNaNO3).GT.0.D0.OR.
     &        qinti(SNH4NO3).GT.0.D0) THEN

            lr56=.true.
         ENDIF
      ENDIF

      IF (leq1.AND.leq2.AND.leq3) THEN
         PRINT *,'Warning from dryin.f: << solid : leq123 >>'
      ENDIF

C     Determine which case is relevant

      IF (leq2.AND.leq3) THEN
         icase=1                ! R2 and R3 active
      ELSEIF (leq1.AND.leq2) THEN
         icase=2                ! R1 and R2 active
      ELSEIF (leq1.AND.leq3) THEN
         icase=3                ! R1 and R3 active
      ELSEIF (leq1) THEN
         icase=4                ! only R1 active
      ELSEIF (leq2) THEN
         icase=5                ! only R2 active
      ELSEIF (leq3) THEN
         icase=6                ! only R3 active

      ELSE                      ! no active equilibrium
         IF (lr47) THEN

            icase=7             ! R4 or R7 active
         ELSEIF (lr56) THEN

            icase=8             ! R5 or R6 active
         ELSE                   ! nothing active

            IF (msat.LT.0.D0) THEN
               icase=9          ! enough nh3 to neutralize so4
            ELSE
               icase=10         ! not enough nh3 to neutralize so4
                                ! in this case aerosol become acidic
            ENDIF
         ENDIF
      ENDIF

C     Solve each case

                                ! icase 3 is not physical but used
                                ! to determine the real icase
      IF (icase.EQ.3) THEN
         cfa=( maa(ENO3)*ckvi(ENO3)
     &        +maa(ECl)*ckvi(ECl)*rk3 )

         cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)
     &        +mkercd(ECl)-mkercd(ENH3)

         cfc=rk1*maa(ENH3)*ckvi(ENH3)

         IF (cfb*cfb+4.D0*cfa*cfc.le.0.D0) THEN
            WRITE(6,*)'(dryin.f): (1) sqrt(<0) '
            STOP
         ENDIF

         qgeqi(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))
     &        /(2.D0*cfa)

                                !qgeqi(ECl)=rk3*qgeqi(ENO3)
                                !qgeqi(ENH3)=rk1/qgeqi(ENO3)

         kercdi(ENO3)= aai(ENO3)
     &        *( qgbki(ENO3)
     &        -qgeqi(ENO3)
     &        *ckvi(ENO3) )

                                ! test if no3 condenses
         IF (kercdi(ENO3).GT.0.D0) THEN
            icase=2             ! if nh4no3 forms then real case=2
         ELSE
            icase=1             ! if nacl forms then real case=1
         ENDIF

      ENDIF

                                ! other cases
      IF (icase.EQ.1) THEN
         cfa=( maa(ENO3)*ckvi(ENO3)
     &        +maa(ECl)*ckvi(ECl)*rk3 )

         cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)
     &        +mkercd(ECl)-mkercd(ENH3)

         cfc=rk2/rk3*maa(ENH3)*ckvi(ENH3)

         IF (cfb*cfb+4.D0*cfa*cfc.le.0.D0) THEN
            WRITE(6,*)'(dryin.f): (2) sqrt(<0) '
            STOP
         ENDIF

         qgeqi(ENO3)= (cfb+DSQRT(cfb*cfb
     &        +4.D0*cfa*cfc))
     &        /(2.D0*cfa)

         qgeqi(ECl)=rk3*qgeqi(ENO3)

         qgeqi(ENH3)=rk2/rk3/qgeqi(ENO3)

      ELSEIF (icase.EQ.2) THEN

         cfa=( maa(ENO3)*ckvi(ENO3)
     &        +maa(ECl)*ckvi(ECl)*rk2/rk1 )

         cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)
     &        +mkercd(ECl)-mkercd(ENH3)

         cfc=rk1*maa(ENH3)*ckvi(ENH3)

         IF (cfb*cfb+4.D0*cfa*cfc.le.0.D0) THEN
            WRITE(6,*)'(dryin.f): (3) sqrt(<0) '
            STOP
         ENDIF

         qgeqi(ENO3)= (cfb+DSQRT(cfb*cfb
     &        +4.D0*cfa*cfc))
     &        /(2.D0*cfa)
         qgeqi(ENH3)=rk1/qgeqi(ENO3)
         qgeqi(ECl)=rk2/rk1*qgeqi(ENO3)

      ELSEIF (icase.EQ.4) THEN
         cfa=maa(ENO3)*ckvi(ENO3)
         cfb= 2.D0*mkercd(ESO4)
     &        +mkercd(ENO3)
     &        -mkercd(ENH3)
         cfc=rk1*maa(ENH3)*ckvi(ENH3)

         IF (cfb*cfb+4.D0*cfa*cfc.le.0.D0) THEN
            WRITE(6,*)'(dryin.f): (4) sqrt(<0) '
            STOP
         ENDIF

         qgeqi(ENO3)= (cfb+DSQRT(cfb*cfb
     &        +4.D0*cfa*cfc))
     &        /(2.D0*cfa)

         qgeqi(ENH3)=rk1/qgeqi(ENO3)

         qgeqi(ECl)=qgbki(ECl)/ckvi(ECl)

      ELSEIF (icase.EQ.5) THEN
         cfa=maa(ENH3)*ckvi(ENH3)
         cfb= 2.D0*mkercd(ESO4)+mkercd(ECl)
     &        -mkercd(ENH3)
         cfc=maa(ECl)*rk2*ckvi(ECl)

         IF (cfb*cfb+4.D0*cfa*cfc.le.0.D0) THEN
            WRITE(6,*)'(dryin.f): (5) sqrt(<0) '
            STOP
         ENDIF

         qgeqi(ENH3)= (cfb+DSQRT(cfb*cfb
     &        +4.D0*cfa*cfc))
     &        /(2.D0*cfa)

         qgeqi(ECl)=rk2/qgeqi(ENH3)

         qgeqi(ENO3)=qgbki(ENO3)/ckvi(ENO3)

      ELSEIF (icase.EQ.6) THEN
         cfa=( maa(ENO3)*ckvi(ENO3)
     &        +rk3*maa(ECl)*ckvi(ECl) )
         cfb=2.D0*mkercd(ESO4)+mkercd(ECl)
     &        +mkercd(ENO3)

         qgeqi(ENO3)=cfb/cfa

         qgeqi(ECl)=rk3*qgeqi(ENO3)

         qgeqi(ENH3)=qgbki(ECl)/ckvi(ECl)

      ELSEIF (icase.EQ.7) THEN
         qgeqi(ENH3)=qgbki(ENH3)/ckvi(ENH3)
         qgeqi(ENO3)=qgbki(ENO3)/ckvi(ENO3)

         kercdi(ENH3)=0.D0
         kercdi(ENO3)=0.D0

         kercdi(ECl)=-EMW(ECl)*2.D0
     &        *mkercd(ESO4)
         qgeqi(ECl)= ( qgbki(ECl)
     &        -kercdi(ECl)/aai(ECl) )
     &        /ckvi(ECl)

      ELSEIF (icase.EQ.8) THEN
         qgeqi(ENH3)=qgbki(ENH3)/ckvi(ENH3)
         qgeqi(ECl)=qgbki(ECl)/ckvi(ECl)

         kercdi(ENH3)=0.D0
         kercdi(ECl)=0.D0

         kercdi(ENO3)=-EMW(ECl)*2.D0
     &        *mkercd(ESO4)

         qgeqi(ENO3)= ( qgbki(ENO3)
     &        -kercdi(ENO3)/aai(ENO3) )
     &        /ckvi(ENO3)

      ELSEIF (icase.EQ.9) THEN
         qgeqi(ENO3)=qgbki(ENO3)/ckvi(ENO3)
         qgeqi(ECl)=qgbki(ECl)/ckvi(ECl)

         kercdi(ENO3)=0.D0
         kercdi(ECl)=0.D0

         kercdi(ENH3)= EMW(ENH3)*2.D0
     &        *mkercd(ESO4)

         qgeqi(ENH3)= ( qgbki(ENH3)
     &        -kercdi(ENH3)/aai(ENH3) )
     &        /ckvi(ENH3)

      ELSEIF (icase.EQ.10) THEN
         qgeqi(ENO3)=qgbki(ENO3)/ckvi(ENO3)
         qgeqi(ECl)=qgbki(ECl)/ckvi(ECl)

         kercdi(ENO3)=0.D0
         kercdi(ECl)=0.D0

C     no more electroneutrality in this case

      ENDIF

C     Giving out the kernel for case <=6

      IF (icase.LE.6) THEN
         DO jj=3,nesp_isorropia
            jesp=isorropia_species(jj)
            kercdi(jesp)= aai(jesp)
     &           *( qgbki(jesp)
     &           -qgeqi(jesp)
     &           *ckvi(jesp) )
         END DO
      ENDIF

      END

