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

      SUBROUTINE KERCOND(nesp_aer,qn,qext,qgbki,qinti,
     &     vawi,dawi,qgeqi,aai,ckvi,kercdi)

C----------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the condensation/evaporation kernels for
C     all semi-volatile species in a given bin.
C     
C----------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     QN   : number aerosol concentration     ([#aero.m-3]).
C     QEXT : external aerosol concentration   ([µg.m-3]).
C     QGBKI: bulk gas concentration           ([µg.m-3]).
C     QINTI: internal inorganic concentration ([µg.m-3]).
C     VAWI : volume aerosol concentration     ([µm3.m-3]).
C     DAWI : wet aerosol diameter             ([µm]).
C     
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     QGEQI: equilibrium gas concentration     ([µg.m-3]).
C     
C     
C     -- OUTPUT VARIABLES
C     
C     AAI    : c/e kernel coefficient          ([m3.s-1]).
C     CKVI   : Kelvin effect coefficient       ([]).
C     KERCDI : c/e kernel                      ([µg.s-1]).
C     
C----------------------------------------------------------------------
C     
C     -- REMARKS
C     
C----------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C     
C----------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2004: Edouard Debry, CEREA.
C     
C----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'vara.inc'

      INTEGER nesp_aer
      DOUBLE PRECISION qn,qext(nesp_aer),qgbki(nesp_aer)
      DOUBLE PRECISION qinti(NINTIS),aai(nesp_aer)
      DOUBLE PRECISION qgeqi(nesp_aer),kercdi(nesp_aer)
      DOUBLE PRECISION ckvi(nesp_aer),vawi,dawi,dawic,rhop

      INTEGER jesp,nsize
      DOUBLE PRECISION qih,emw_tmp,rhop_tmp

C     ******Initialization to zero
      DO jesp=1,nesp_aer
         aai(jesp)=0.D0
         ckvi(jesp)=1.D0        ! kelv coef init is 1.D0 (KS) 
         kercdi(jesp)=0.D0
      END DO

C     ******c/e kernel coefficient
      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
            CALL COMPUTE_CONDENSATION_TRANSFER_RATE(
     $           DIFFG(jesp),      ! diffusion coef (m2.s-1)
     $           VQMG(jesp),       ! quadratic mean speed (m.s-1)
     $           STICK(jesp),      ! accomadation coef (adim)
     $           dawi,             ! wet aero diameter (µm)
     $           aai(jesp) )       ! c/e kernel coef (m3.s-1)
         ENDIF
      END DO

C     ******Kelvin effect
      IF (IKV2.EQ.1) THEN
         dawic =DMAX1(dawi,Dmin)
C     Aerosol density in µ g.µ m -3
         rhop = 0.d0
         DO jesp= 1,nesp_aer
            rhop=rhop+qext(jesp)
         enddo         
         
         rhop =rhop/vawi        ! µg.µm-3
         rhop_tmp = rhop * 1.D12 ! kg/m3

         DO jesp=E1,E2
            IF (aerosol_species_interact(jesp).GT.0) THEN
               emw_tmp = EMW(jesp) * 1.D-6 ! g/mol
               CALL COMPUTE_KELVIN_COEFFICIENT(
     $              TEMP,          ! temperature (Kelvin)
     $              emw_tmp,       ! ext mol weight (g.mol-1)
     $              SIGMA(jesp),   ! surface tension (N.m-1)
     $              dawic,         ! wet aero diameter (µm)
     $              rhop_tmp,      ! aerosol density (kg.m-3)
     $              ckvi(jesp) )   ! kelvin effect coef (adim)
            ENDIF
         ENDDO
      ENDIF
      
C     ******Not limited c/e kernels
      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
            kercdi(jesp)= aai(jesp) ! kernel coef (m3.s-1)
     &            *( qgbki(jesp)    ! bulk gas conc (µg.m-3)
     &            -qgeqi(jesp)      ! equi gas conc (µg.m-3)
     &            *ckvi(jesp) )     ! kelvin coef (adim) 
         ENDIF
      ENDDO

C     ******Limited inorganic kernels
      nsize = nesp_aer
                                ! aerosol liquid water content
      IF (qext(EH2O).EQ.0.D0) THEN ! solid
         
         CALL DRYIN( TEMP,      ! local temperature (Kelvin)
     &        qinti,            ! int sld inorg conc (µg.m-3)
     &        nsize,            ! size of vectors following below
     &        qgbki,            ! bulk gas conc (µg.m-3)
     &        aai,              ! kernel coef (m3.s-1)
     &        ckvi,             ! kelvin coef (adim)
     &        qgeqi,            ! equi gas conc (µg.m-3)
     &        kercdi)           ! modified c/e kernel
                                ! liq or mix : H+ limitation flux
      ELSE
         qih=qinti(IH)/qn       ! µg
         
         CALL HPLFLIM( ALFHP,   ! percentage of H+ allowed to c/e
     &        qih,              ! int H+ conc (µg)
     &        nsize,            ! size of vectors following below
     &        qgbki,            ! bulk gas conc (µg.m-3)
     &        aai,              ! kernel coef (m3.s-1)
     &        ckvi,             ! kelvin coef (adim)
     &        qgeqi,            ! equi gas conc (µg.m-3)
     &        kercdi)           ! modified c/e kernel
      ENDIF
      
      END


