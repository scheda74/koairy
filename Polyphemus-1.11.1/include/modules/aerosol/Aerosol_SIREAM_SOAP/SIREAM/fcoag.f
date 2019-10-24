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

      SUBROUTINE FCOAG(neq,nesp_aer,nbin_aer,q,iq,dqdt,couples_coag,
     s      first_index_coag,second_index_coag,coefficient_coag,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD, bin_density)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the source terms for the system of
C     Ordinary Differential Equations defined by the coagulation.    
C     
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C     
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     Q   : gas/aerosol concentrations ([\mu.g.m^-3]).
C     ZSC: volumic emissions for aerosol ([\mu g.m^-3.s^-1]).
C     iq: index of aerosol species in q(*) vector.
C     couples_coag: coagulation couples for each bin.
C     first_index_coag: first index of coagulation couples.
C     second_index_coag: second index of coagulation couples.
C     coefficient_coag: coagulation partition coefficient ([adim]).
C     XSF: neperian logarithm of fixed aerosol bin mass ([adim]).
C     MSF: fixed aerosol bin dry mass ([\mu g]).
C     DSF: fixed aerosol bin dry diameter ([\mu m]).
C     
C     -- INPUT/OUTPUT VARIABLES
C
C     Q   : gas/aerosol concentrations ([\mu.g.m^-3]).
C     XSD: neperian logarithm of moving aerosol bin mass ([adim]).
C     MSD: moving aerosol bin dry mass ([\mu g]).
C     DSD: moving aerosol bin dry diameter ([\mu m]).
C     QT: total aerosol concentration per bin ([\mu g.m^-3]).
C     
C     -- OUTPUT VARIABLES
C     
C     DQDT : time derivative ([µg.m-3.s-1]).    
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
C     2005/10: Add ITHRM hook to compute diameter (Edouard Debry, CEREA).
C     2005/10: Put KERCOAG routine (Edouard Debry, CEREA).
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
      INCLUDE 'meteo.inc'
      INCLUDE 'varg.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq),dqdt(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER js,jj,jesp,j1,j2
      DOUBLE PRECISION gain(neq),loss(nbin_aer)
      DOUBLE PRECISION qnum(nbin_aer),qesp(nbin_aer)
      DOUBLE PRECISION gain2(nbin_aer)

      DOUBLE PRECISION kercg(nbin_aer,nbin_aer)
      DOUBLE PRECISION fcg(nbin_aer,nbin_aer,nbin_aer)

      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(nesp_aer,nbin_aer),RHOW(nbin_aer)
      DOUBLE PRECISION bin_density(nbin_aer)
C     Initialization to zero.

      DO jj=1,neq
         gain(jj)=0.D0
      END DO

C     Local gas/aerosol equilibrium needed for computing the
C     aerosol diameter.

      IF (ITHRM.EQ.0) THEN
C     always recompute diam with isorropia
         CALL EQPART(nesp_aer,nbin_aer,1,nbin_aer,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,
     s      QINT,QGEQ,VSW,MSW,DSW,RHOW)
      ELSEIF (ITHRM.EQ.1) THEN
C     always recompute diam with fastdiam (Gerber)
         CALL FASTDIAM(nesp_aer,nbin_aer,1,nbin_aer,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,VSW,MSW,DSW,RHOW,bin_density)
      ENDIF

C     Coagulation kernels
      DO j2=1,nbin_aer
         DO j1=1,j2
            CALL COMPUTE_BIDISPERSE_COAGULATION_KERNEL(TEMP,AIRFMP,
     &           DSW(j1),DSW(j2),
     &           MSW(j1),MSW(j2),
     &           kercg(j1,j2))
            
                                ! symmetric kernels
            kercg(j2,j1)=kercg(j1,j2)
         END DO
      END DO
      
C     Compute coagulation coefficients
      DO js=1,nbin_aer
         DO jj=1,couples_coag(js)
            j1=first_index_coag(js,jj)
            j2=second_index_coag(js,jj)
            
            fcg(j1,j2,js)= coefficient_coag(j1,j2,js)
     &           *kercg(j1,j2)

                                ! symmetric coef
            fcg(j2,j1,js)=fcg(j1,j2,js)
         END DO
      END DO

C     Compute coagulation gain

                                ! for number distribution
      DO js=1,nbin_aer
         qnum(js)=q(js)
      END DO

      CALL CALCGAIN(nbin_aer,qnum,qnum,fcg,gain2,couples_coag,
     s      first_index_coag,second_index_coag)

      DO js=1,nbin_aer
         gain(js)=gain2(js)
      END DO

                                ! for each species
      DO jesp=E1,E2
         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            qesp(js)=q(jj)
         END DO

         CALL CALCGAIN(nbin_aer,qnum,qesp,fcg,gain2,couples_coag,
     s      first_index_coag,second_index_coag)

         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            gain(jj)=gain2(js)
         END DO
      END DO

C     Compute coagulation loss

      CALL CALCLOSS(nbin_aer,qnum,kercg,loss)

C     Update time derivatives      

                                ! for number distribution
      DO js=1,nbin_aer
         dqdt(js)= dqdt(js)
     &        +5.D-01*gain(js)
     &        -q(js)*loss(js)
      END DO
      
                                ! for each species
      DO jesp=E1,E2
         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            dqdt(jj)= dqdt(jj)
     &           +gain(jj)
     &           -q(jj)*loss(js)
         END DO
      END DO

      END

