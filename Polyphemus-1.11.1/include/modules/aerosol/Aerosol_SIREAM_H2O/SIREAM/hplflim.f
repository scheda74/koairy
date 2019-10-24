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

      SUBROUTINE HPLFLIM(alfa,qih,nsize,qgbki,aai,ckvi,qgeqi,kercdi)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the flux limitation for the
C     condensation/evaporation flux. The algorithm is based
C     on the limitation of the aerosol acidity rate.
C
C     The details may be found in the PhD Work of Edouard Debry,
C     Chapter 10 (section 10.1.5) or in the reference:
C     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
C     Aerosol Science and Technology 32, 482:502, 2000.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C     alfa  : percentage (0.0< <1.0) of H+ aerosol mass allowed to c/e.
C     QIH   : internal H+ aerosol mass ([µg]).
C     nsize : size of vectors following below.
C     QGBKI : bulk gas concentration ([µg.m-3]).
C     AAI   : condensation/evaporation kernel coefficient([m3.s-1]).
C
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     CKVI  : kelvin effect coefficient ([]).
C     QGEQI : gas-phase equilibrium concentration ([µg.m-3]).
C     KERCDI: condensation/evaporation kernel ([µg.s-1]).
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
      INCLUDE 'imw.inc'
      INCLUDE 'vara.inc'
      INCLUDE 'varp.inc'

      INTEGER nsize
      DOUBLE PRECISION qgbki(nsize),aai(nsize)
      DOUBLE PRECISION ckvi(nsize),qgeqi(nsize)
      DOUBLE PRECISION kercdi(nsize),qih,alfa

      INTEGER jesp,jj
      DOUBLE PRECISION maa(nsize),mkercd(nsize)
      DOUBLE PRECISION cfa,cfb,cfc,cc
      DOUBLE PRECISION mih,melec,mlim,q

C     Compute mol fluxes

      DO jj=1,nesp_isorropia
         jesp=isorropia_species(jj)
                                ! maa(*) in m3.mol.s-1.µg-1
         maa(jesp)= aai(jesp)   ! m3.s-1
     &        /EMW(jesp)        ! µg.mol-1

                                ! mkercd(*) in mol.s-1
         mkercd(jesp)= kercdi(jesp) ! µg.s-1
     &        /EMW(jesp)        ! µg.mol-1
      END DO

C     H+ limitation

      mih=qih/IMW(IH)           ! mol of H+ in aerosol

                                ! maximum of mih variation tolerated
      mlim=mih*alfa             ! mol.s-1

                                ! electroneutrality  ! mol.s-1
      melec= 2.D0*mkercd(ESO4)+mkercd(ENO3)
     &     +mkercd(ECl)-mkercd(ENH3)

                                ! correction factor default value
      cc=0.D0

                                ! correction calculation
      IF (DABS(melec).GT.mlim) THEN
                                ! we give to mlim the sign of melec
         mlim=DSIGN(mlim,melec)

                                ! cfa,cfb,cfc are coefficients of
                                ! 2nd order eq : cfa*cc^2+cfb*cc+cfc=0
                                ! satisfied by the correction factor cc

         cfa= maa(ENO3)*qgeqi(ENO3)*ckvi(ENO3)
     &        +maa(ECl)*qgeqi(ECl)*ckvi(ECl)

         cfb=  mlim
     &        -2.D0*maa(ESO4)*qgbki(ESO4) ! mol.s-1
     &        -maa(ENO3)*qgbki(ENO3)
     &        -maa(ECl) *qgbki(ECl)
     &        +maa(ENH3)*qgbki(ENH3)

         cfc=-maa(ENH3)*qgeqi(ENH3)*ckvi(ENH3)

                                ! one can note cfa>=0 and cfc<=0
                                ! then there always exist a + but
                                ! possibly zero root

                                ! root computation
         IF (cfa.GT.0.D0) THEN
            IF (cfb*cfb-4.D0*cfa*cfc.le.0.D0) THEN
               WRITE(6,*)'(hplflim.f): sqrt(<0)'
               STOP
            ENDIF

            q=-5.D-01*(cfb+DSIGN(1.D0,cfb)
     &           *DSQRT(cfb*cfb-4.D0*cfa*cfc))
            cc=DMAX1(q/cfa,cfc/q) ! we select the + root
         ELSE
            IF (cfb.NE.0.D0) cc=-cfc/cfb
         ENDIF
      ENDIF

C     A correction is done if only upper calculation
C     root has changed cc to a strictly positive
C     value, otherwise it is considered as non stiff cases

      IF (cc.GT.0.D0) THEN
         qgeqi(ENH3)=qgeqi(ENH3)/cc
         qgeqi(ENO3)=qgeqi(ENO3)*cc
         qgeqi(ECl) =qgeqi(ECl)*cc

         DO jj=3,nesp_isorropia
            jesp=isorropia_species(jj)
            kercdi(jesp)= aai(jesp)*(qgbki(jesp)-qgeqi(jesp)*ckvi(jesp))
         END DO
      ENDIF

      END
