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

      SUBROUTINE STEP(qn,qext,qti,qinti,
     &     qgeqi,vaw,dad,daw,rhoaer)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the local aerosol equilibrium in each bin.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     QN :   number aerosol concentration   ([#aero.m^-3]).
C     QEXT : external aerosol concentration ([\mu.g.m^-3]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     QTI   : total dry concentration          ([\mu.g.m^-3]).
C     QGEQI : equilibrium gas concentration    ([\mu.g.m^-3]).
C     VAW   : wet volume aerosol concentration ([\mu.m^3.m^-3]).
C     DAD   : dry aerosol dimaeter             ([\mu.m]).
C     DAW   : wet aerosol diameter             ([\mu.m]).
C     RHOAER: aerosol density                  ([\mu.g.\mu.m^-3]).
C     QINTI : internal inorganic concentration ([\mu.g.m^-3]).
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
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'pointer.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'varp.inc'

      DOUBLE PRECISION qn,qext(NEXT),qti
      DOUBLE PRECISION qinti(NINTIS),qgeqi(NEXT)
      DOUBLE PRECISION vaw,dad,daw
      DOUBLE PRECISION qtinorg,qtorg

      INTEGER jesp
      DOUBLE PRECISION vad,rhoaer

C     ******zero init
      DO jesp=1,NEXT
         qgeqi(jesp)=0.D0
      END DO

C     ******total dry mass
      qti=0.D0
      DO jesp=E1,E2
         qti=qti+qext(jesp)     ! µg.m-3
      END DO
C     ******inorganics thermodynamics

      DO jesp=1,NINTIS
         qinti(jesp)=0.d0
      ENDDO
      qext(EH2O)=0.d0

                                ! sum of inorganic mass
      qtinorg=0.D0
      DO jesp=ENa,ECl
         qtinorg=qtinorg+qext(jesp)
      END DO

                                ! if no inorg mass avoid calling thermo
      IF (qtinorg.GT.TINYM) THEN
         CALL EQINORG( qext,         ! ext inorg aero conc (µg.m-3)
     &                qinti,         ! int inorg aero conc (µg.m-3)
     &                qgeqi )        ! inorg eq gas conc (µg.m-3)
      ENDIF

C     ******organics thermodynamics
                                ! sum of organic mass
      qtorg=0.D0
      DO jesp=EARO1,EPOA
         qtorg=qtorg+qext(jesp)
      END DO

                                ! if no org mass avoid calling thermo
      IF (qtorg.GT.TINYM) THEN
         CALL EQORG( qext,   ! ext org aero conc (µg.m-3)
     &        qgeqi )        ! org pure eq gas conc (µg.m-3)
      ENDIF

C     ******aerosol volume and diameter
                                ! aerosol volume computed
                                ! with fixed or variable density
      IF (IDENS.EQ.0) THEN
                                ! fixed aerosol density µg.µm-3
         rhoaer=RHOA
      ELSE
                                ! moving aerosol density
         CALL VOLAERO(
     &        qext,             ! extern aero conc (µg.m-3)
     &        qinti,            ! int inorg conc (µg.m-3)
     &        rhoaer)           ! aerosol density (µg.µm-3)
      ENDIF

      vad=qti/rhoaer
      vaw=vad+qext(EH2O)/rhoaer
                                ! aerosol diameter
                                ! qn cannot be zero, checked in eqpart routine
      dad=(vad/qn/cst_pi6)**cst_FRAC3 ! µm
      daw=(vaw/qn/cst_pi6)**cst_FRAC3 ! µm

      END

