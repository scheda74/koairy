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

      SUBROUTINE FGDE(neq,nesp_aer,nbin_aer,q,iq,dqdt,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA, bin_density)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the source terms for the system of
C     Ordinary Differential Equations defined by the GDE.    
C     
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C     
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     Q   : gas/aerosol concentrations ([\mu.g.m^-3]).
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
C     AA: condensation coefficient             ([m^3.s^-1]).
C     
C     -- OUTPUT VARIABLES
C     
C     DQDT : time derivative ([�g.m-3.s-1]).
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
      INCLUDE 'varp.inc'
      
      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq),dqdt(neq)
      INTEGER iq(nesp_aer,nbin_aer)
      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)

      DOUBLE PRECISION AA(nesp_aer,nbin_aer)

      INTEGER jj,     js, jesp

      DOUBLE PRECISION bin_density(nbin_aer)

C     Initialization to zero.

      DO jj=1,neq
         dqdt(jj)=0.D0
      END DO

C     Compute the physical processes

      IF (ICG2.EQ.1) then
c$$$         DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN avant fcoag"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO

         CALL FCOAG(neq,nesp_aer,nbin_aer,q,iq,dqdt,
     s     couples_coag,first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD, bin_density)
c$$$         DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN apres fcoag"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO
      ENDIF

      IF (ICE2.EQ.1) then
c$$$          DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN avant fcond"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO
         
         CALL FCOND(neq,nesp_aer,nbin_aer,q,iq,dqdt,QT,
     s     XSF,MSF,DSF,XSD,MSD,DSD,AA)
c$$$          DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN apres fcond"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO
      ENDIF

      IF (INU2.EQ.1) then
c$$$         DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN avant fnucl"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO

         CALL FNUCL(neq,nesp_aer,nbin_aer,q,iq,dqdt)
c$$$         DO js = 1,nbin_aer
c$$$            DO jesp = 1,nesp_aer-1
c$$$               if (isNaN(q(IQ(jesp,js)))) then
c$$$                  print *, js, jesp, "q(IQ) NaN apres fnucl"
c$$$               endif
c$$$            ENDDO
c$$$         ENDDO
      ENDIF

      END
