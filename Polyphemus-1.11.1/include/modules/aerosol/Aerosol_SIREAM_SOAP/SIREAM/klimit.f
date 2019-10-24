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

      SUBROUTINE KLIMIT(neq,nesp_aer,nbin_aer,q,iq,k,AA)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine limits the condensation/evaporation rates for
C     aerosol and gases in order to avoid clippings.     
C     Two kinds of limitation are performed:
C     - the 1st is aerosol clipping : as it may reduce evaporation then 
C     enlarge condensation;, it is done before condensation limitation.
C     - the 2nd is gas clipping : in practice it may reduce, per species,
C     aerosol condensation only in bins that leads to gas clipping.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     iq: aerosol pointers.
C     Q   : gas/aerosol concentration ([µg.m-3]).
C     AA: condensation coefficients ([s^-1]).
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     K   :  c/e time derivatives ([µg.m-3.s-1]).
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
      INCLUDE 'dynaero.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varq.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq),k(neq)
      INTEGER iq(nesp_aer,nbin_aer)
      DOUBLE PRECISION AA(nesp_aer,nbin_aer)

      INTEGER js,jesp,jj
      DOUBLE PRECISION ksum,klim,ktlim
      DOUBLE PRECISION qsum,aatot
      DOUBLE PRECISION frac,qnew

C     ****** prevent aerosol clipping
      DO jesp=E1,E2
         DO js=ICUT+1,nbin_aer
            jj=IQ(jesp,js)      

                                ! only when evaporation
            IF (k(jj).LT.0.D0) THEN

                                ! if q(jj) is <=TINYM or =0
                                ! then k should be >=0, but
                                ! due to bad matrix inversion
                                ! this case may occur
               IF (q(jj).LE.TINYM) k(jj)=0.D0

                                ! test clipping in other cases

               qnew=q(jj)+k(jj)          
               IF (qnew.LT.0.D0) THEN
                                ! we are sure that q>=tinym
                                ! otherwise k would be = 0
                                ! from previous case
                  k(jj)=(TINYM-q(jj)) !/dt ! <=0 µg.m-3.s-1

                                ! we force q to be a
                                ! 'little' more than TINYM
                  k(jj)=0.99D0*k(jj)

               ENDIF
            ENDIF
         END DO
      END DO

C     ****** prevent gas clipping
      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
                                ! compute total mass rate per species
         ksum=0.D0
         DO js=ICUT+1,nbin_aer
            jj=IQ(jesp,js)
            ksum=ksum+k(jj)     ! µg.m-3.s-1
         END DO

                                ! we perform limiting in
                                ! case of condensation only
         IF (ksum.GT.0.D0) THEN
            
                                ! this is the total lumped mass
                                ! to perserve from clipping

            jj=IG(jesp)
            qsum=q(jj)
            DO js=1,ICUT
               jj=IQ(jesp,js)
               qsum=qsum+q(jj)  ! µg.m-3
            END DO

                                ! test if clipping occurs
                                ! then perform the limitation

            IF (ksum.GT.qsum) THEN
                                ! sum of aa(*) c/e coefficient
               aatot=0.D0
               DO js=1,nbin_aer
                  aatot=aatot+AA(jesp,js) 
               END DO

                                ! tot max rate, µg.m-3.s-1
               ktlim=qsum       !/dt

               DO js=ICUT+1,nbin_aer
                                ! fraction, adim
                  frac=AA(jesp,js)/aatot ! aatot != 0

                                ! we allow only a given fraction
                                ! of ktlim to condense on given bin
                  klim=ktlim*frac

                                ! apply the limit
                                ! only if necessary
                  jj=IQ(jesp,js)

                  IF (k(jj).GT.klim) k(jj)=klim

               END DO


            ENDIF
         ENDIF
         ENDIF
      END DO

      END
