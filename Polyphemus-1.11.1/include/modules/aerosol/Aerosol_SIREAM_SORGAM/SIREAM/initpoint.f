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

      SUBROUTINE INITPOINT(nbin_aer,iq)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine performs the initialization of all the pointers
C     related to the Q(*) vector at current time.
C     All the variables are in common blocks.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nbin_aer: number of aerosol bins.
C     iq: index of aerosol species in q(*) vector.
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
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
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'

      INTEGER nbin_aer,iq(NEXT,nbin_aer)
      INTEGER jesp,js,icpt

C     Species to be computed

      E1=EMD
      E2=EPOA
      G1=ESO4
      G2=EBiBmP

C     Pointer for concentrations

      icpt=nbin_aer

                                ! pointer of species at equilibrium
      DO jesp=E1,E2
         DO js=1,ICUT
            icpt=icpt+1
            IQ(jesp,js)=icpt
         END DO
      END DO

                                ! pointer of dynamic species
      DO jesp=E1,E2
         DO js=ICUT+1,nbin_aer
            icpt=icpt+1
            IQ(jesp,js)=icpt
         END DO
      END DO

                                ! H2O pointers
      DO js=1,nbin_aer
         icpt=icpt+1
         IQ(EH2O,js)=icpt
      END DO

C     Pointer for gas-phase concentration

      DO jesp=1,NEXT
         icpt=icpt+1
         IG(jesp)=icpt
      END DO

      IG1=IG(1)

      END
