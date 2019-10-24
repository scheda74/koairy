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

      SUBROUTINE EQINORG(nesp_aer,qext,qinti,qgeqi)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the equilibrium between inorganic aerosols
C     and gas-phase (reverse mode).
C     It calls ISORROPIA by Nenes et al.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     QEXT : external aerosol concentration ([µg.m-3]).
C     
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     QINTI : internal inorganic aerosol concentration ([µg.m-3]).
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
      INCLUDE 'meteo.inc'
      INCLUDE 'vara.inc'
      INCLUDE 'imw.inc'
      INCLUDE 'smw.inc'
      INCLUDE 'varp.inc'

      DOUBLE PRECISION GREAT
      PARAMETER (GREAT=100.D0)

      INTEGER nesp_aer
      DOUBLE PRECISION qext(nesp_aer),qinti(NINTIS)
      DOUBLE PRECISION qgeqi(nesp_aer)

      INTEGER jesp,jj
      DOUBLE PRECISION wi(5),w(5)
      DOUBLE PRECISION aerliq(12),aersld(9)
      DOUBLE PRECISION gas(3),cntrl(2)
      DOUBLE PRECISION other(6),coef
C      CHARACTER*15 scase
      DOUBLE PRECISION organion,watorg

C no SOA in reverse mode for that moment
      organion = 0.D0
      watorg = 0.D0

C     Inputs  for Isorropia

      gas(1)=0.d0
      gas(2)=0.d0
      gas(3)=0.d0

      cntrl(1)=1.D0             ! reverse mode
      cntrl(2)=MTSBL            ! metastable option

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      coef=DMAX1(qext(ESO4),
     &     qext(ENH3),qext(ENO3),
     &     TINYM)/GREAT
#else
      coef=DMAX1(qext(ENa),qext(ESO4),
     &     qext(ENH3),qext(ENO3),qext(ECl),
     &     TINYM)/GREAT
#endif
                                ! convert µg to moles
                                ! and scaling for reverse mode
      DO jj=1,nesp_isorropia
         jesp=isorropia_species(jj)
         wi(jj)= qext(jesp) ! µg.m-3
     &        /EMW(jesp)  ! µg.mol-1
     &        /coef             ! adim
      END DO

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      wi(1) = 0.D0              !Do not consider sea salt in isoropia
      wi(5) = 0.D0
#endif

      CALL ISOROPIA(wi,RH,TEMP,cntrl,w,gas,
     &     aerliq,aersld,other,organion,watorg)

C     clipping to tinym

      IF (gas(1).lt.0.d0) gas(1)=tinym
      IF (gas(2).lt.0.d0) gas(2)=tinym
      IF (gas(3).lt.0.d0) gas(3)=tinym
      
C     Outputs isorropia

                                ! sulfate surf conc always 0. µg.m-3
      qgeqi(ESO4)=0.D0

                                ! convert moles.m-3 to µg.m-3 
      qgeqi(ENH3)=gas(1)*EMW(ENH3)
      qgeqi(ENO3)=gas(2)*EMW(ENO3)
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
      qgeqi(ECl) =gas(3)*EMW(ECl)
#endif


                                ! liquid inorg aerosol
      DO jesp=IH,IOH
         qinti(jesp)= DMAX1(aerliq(jesp),0.D0)
     &        *coef             ! scaling back
     &        *IMW(jesp)        ! moles to µg
      END DO

                                ! liquid water content
      qext(EH2O)= qinti(IH2O)
     &     +qinti(IOH)*1.05882352941D0 ! mwh2o/mwioh

                                ! solid inorg aerosol
      DO jesp=SNaNO3,SLC
         qinti(jesp)= DMAX1(aersld(jesp-12),0.D0)
     &        *coef             ! scaling back
     &        *SMW(jesp)        ! moles to µg
      END DO

      END
