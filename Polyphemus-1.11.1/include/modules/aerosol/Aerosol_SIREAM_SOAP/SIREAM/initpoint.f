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

      SUBROUTINE INITPOINT(nbin_aer,nesp_aer,
     $     aerosol_species_interact_loc,iq,
     $     nesp_isorropia_loc,isorropia_species_loc, 
     $     nesp_aec_loc,aec_species_loc,
     $     nesp_pankow_loc,pankow_species_loc,
     $     nesp_poa_loc,poa_species_loc,
     $     md_species_loc, bc_species_loc,
     $     nesp_cloud_interact, cloud_species_interact)

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
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varp_cloud.inc'


      INTEGER nesp_aer,nbin_aer,iq(nesp_aer,nbin_aer)
      INTEGER jesp,js,icpt
      INTEGER aerosol_species_interact_loc(nesp_aer)
      INTEGER nesp_isorropia_loc
      INTEGER isorropia_species_loc(nesp_isorropia_loc)
      INTEGER nesp_aec_loc,aec_species_loc(nesp_aec_loc)
      INTEGER nesp_pankow_loc,pankow_species_loc(nesp_pankow_loc)
      INTEGER nesp_poa_loc,poa_species_loc(nesp_poa_loc)
      INTEGER md_species_loc,bc_species_loc
      INTEGER nesp_cloud_interact
      INTEGER cloud_species_interact(nesp_cloud_interact)

C     Species to be computed

      E1=1
      E2=nesp_aer-1 ! to avoid water

      nesp_isorropia=nesp_isorropia_loc
      nesp_aec=nesp_aec_loc
      nesp_pankow=nesp_pankow_loc
      nesp_pom=nesp_poa_loc

      IF (nesp_isorropia.GT.NEXT) THEN
         STOP 'initpoint.f : nesp_isorropia > NEXT'
      ENDIF
      IF (nesp_aec.GT.NEXT) THEN
         STOP 'initpoint.f : nesp_aec > NEXT'
      ENDIF
      IF (nesp_pankow.GT.NEXT) THEN
         STOP 'initpoint.f : nesp_pankow > NEXT'
      ENDIF
      IF (nesp_pom.GT.NEXT) THEN
         STOP 'initpoint.f : nesp_pankow > NEXT'
      ENDIF

      DO jesp=1,NEXT
         isorropia_species(jesp) = -1
         aec_species(jesp) = -1
         pankow_species(jesp) = -1
         poa_species(jesp) = -1
         aerosol_species_interact(jesp) = -1
      ENDDO

      DO jesp=1,nesp_isorropia
         isorropia_species(jesp) = isorropia_species_loc(jesp)
      ENDDO

      DO jesp=1,nesp_aec
         aec_species(jesp) = aec_species_loc(jesp)
      ENDDO

      DO jesp=1,nesp_pankow
         pankow_species(jesp) = pankow_species_loc(jesp)
      ENDDO

      DO jesp=1,nesp_pom
         poa_species(jesp) = poa_species_loc(jesp)
      ENDDO

      ENa=isorropia_species(1)
      ESO4=isorropia_species(2)
      ENH3=isorropia_species(3)
      ENO3=isorropia_species(4)
      ECl=isorropia_species(5)
      EMD=md_species_loc
      EBC=bc_species_loc
      EH2O=nesp_aer ! water always at the end

      DO jesp=1,nesp_aer
         aerosol_species_interact(jesp)=
     &      aerosol_species_interact_loc(jesp)
      ENDDO

C     Pointer for concentrations

      icpt=nbin_aer

      DO jesp=1,nesp_aer
         DO js=1,nbin_aer
            icpt=icpt+1
            IQ(jesp,js)=icpt
         END DO
      END DO

C     Pointer for gas-phase concentration

      DO jesp=1,nesp_aer
         icpt=icpt+1
         IG(jesp)=icpt
      END DO

      IG1=IG(1)

C     Pointers for cloud species.

      ictmNH3=cloud_species_interact(1)
      ictmHNO3=cloud_species_interact(2) 
      ictmHCl=cloud_species_interact(3)
      ictmSO2=cloud_species_interact(4) 
      ictmH2O2=cloud_species_interact(5) 
      ictmHCHO=cloud_species_interact(6) 
      ictmHNO2=cloud_species_interact(7)
      ictmO3=cloud_species_interact(8) 
      ictmOH=cloud_species_interact(9) 
      ictmHO2=cloud_species_interact(10) 
      ictmNO3=cloud_species_interact(11) 
      ictmNO=cloud_species_interact(12) 
      ictmNO2=cloud_species_interact(13)
      ictmPAN=cloud_species_interact(14)
      ictmH2SO4=cloud_species_interact(15)

      END
