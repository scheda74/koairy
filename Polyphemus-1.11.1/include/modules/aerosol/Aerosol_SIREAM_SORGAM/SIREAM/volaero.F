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

      SUBROUTINE VOLAERO(qext,qinti,rhoaer)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the dry and wet aerosol volumes and the
C     aerosol density according to the internal composition.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C
C     QEXT :  external aerosol concentration   ([\mu.g.m^-3]).
C     QINTI:  internal inorganic concentration ([\mu.g.m^-3]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     RHOAER : aerosol density ([\mu.g.\mu.m^-3]).
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
      INCLUDE 'dynaero.inc'
      INCLUDE 'liqmd.inc'
      INCLUDE 'solmd.inc'
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'

      DOUBLE PRECISION qext(NEXT),qinti(NINTIS),rhoaer

      INTEGER jesp
      DOUBLE PRECISION vid,viw,vis,vil,vod,vad,vaw,sumint

C     ******inorganic volume
      vid=0.D0
      vis=0.D0
      vil=0.D0
      viw=0.D0
      vod=0.D0
      sumint=0.d0

C     Mineral Dust
      vis = vis + qext(EMD)/LMD(EMD)
      sumint = sumint + qext(EMD)

C     Black Carbon
      vis = vis + qext(EBC)/LMD(EBC)
      sumint = sumint + qext(EBC)

C     Inorganic
                                ! compute solid aerosol volume
      DO jesp=SNaNO3,SLC
         vis=vis+qinti(jesp)
     &        /SMD(jesp)
      END DO

                                ! compute liquid aerosol volume

                                ! sodium volume
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      vil=vil + qext(ENa)/LMD(ENa)
#endif

#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
      vil=vil + qinti(INa)/LMD(ENa)
#endif

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      vil = vil + qext(ECl)/LMD(ECl) ! HCl volume
#endif
                                ! ammonium volume
      vil=vil+( qinti(INH3)
     &     +qinti(INH4)
     &     *0.944444444444D0    ! mwnh3/mwinh4
     &     )/LMD(ENH3)          ! µg.µm-3

                                ! nitric acid volume
      vil=vil+( qinti(IHNO3)
     &     +qinti(INO3)
     &     *1.01612903226D0     ! mwhno3/mwino3
     &     )/LMD(ENO3)          ! µg.µm-3

                                ! chlorhydric acid volume
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
         vil=vil+( qinti(IHCl)
     &        +qinti(ICl)
     &        *1.02816901408D0  ! mwhcl/mwicl
     &        )/LMD(ECl)        ! µg.µm-3
#endif
                                ! sulfuric acid volume
      vil=vil+( qinti(IHSO4)
     &     *1.01030927835D0     ! mwh2so4/mwihso4
     &     +qinti(ISO4)
     &     *1.02083333333D0     ! mwh2so4/mwiso4
     &     )/LMD(ESO4)          ! µg.µm-3

      vid=vil+vis               ! dry inorg vol µm3.m-3

                                ! water volume
      viw=qext(EH2O)/LMD(EH2O)

                                ! total inorg internal mass

      DO jesp=1,NINTIS
         sumint=sumint+qinti(jesp)
      ENDDO

                                ! correction for water
      sumint=sumint-qinti(IH2O)-qinti(IOH)+qext(EH2O)

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                                ! correction when no NaCl in internal composition
      sumint = sumint - qinti(INa) + qext(ENa)
      sumint = sumint - qinti(IHCl) - qinti(ICl)*1.02816901408D0
     &     + qext(ECl)
#endif

C     dry organic volume and total internal mass
      DO jesp=EARO1,EPOA
         vod=vod+qext(jesp)/LMD(jesp)
         sumint=sumint+qext(jesp)
      END DO

C     ******tot dry and wet vol, µm3.m-3
      vad=vod+vid
      vaw=vad+viw

C     Notice that the density is based on the internal composition
C     and as such bigger than the minimal density and less that the maximal
C     one

      IF ((vaw.GT.0.d0).AND.(sumint.GT.0.D0)) THEN
         rhoaer=sumint/vaw
      ELSE
         rhoaer=RHOA
      ENDIF

      END
