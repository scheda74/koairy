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

      SUBROUTINE FCOND(neq,nesp_aer,nbin_aer,q,iq,dqdt,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,AA)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the source terms for the system of
C     Ordinary Differential Equations defined by condensation/evaporation
C     along a moving grid.
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
C     AA: condensation coefficient             ([m^3.s^-1]).
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
      INCLUDE 'pointer.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq),dqdt(neq)

      INTEGER js,jj,jesp,GLOC,iq(nesp_aer,nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(nesp_aer,nbin_aer),RHOW(nbin_aer)

      DOUBLE PRECISION KERCD(nesp_aer,nbin_aer),CKV(nesp_aer,nbin_aer)
      DOUBLE PRECISION AA(nesp_aer,nbin_aer)

C     Compute gas mass conservation.

      CALL MASSCNSRV(neq,nesp_aer,nbin_aer,q,iq)

C     Compute c/e kernels in dynamic section.
      CALL EQPART(nesp_aer,nbin_aer,ICUT2+1,nbin_aer,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,
     s      QINT,QGEQ,VSW,MSW,DSW,RHOW)
      IKV2=IKELV
      CALL KERCD2(nesp_aer,nbin_aer,ICUT2+1,nbin_aer,neq,q,iq,QT,
     s      QINT,QGEQ,VSW,DSW,AA,CKV,KERCD)

C     Update time derivatives.

      IF(ISULFCOND.EQ.1) THEN
         GLOC = ESO4           ! Do not take sulfate into account
      ELSE
         GLOC = -1
      ENDIF

      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
            IF (jesp.NE.GLOC) THEN
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         IF (jesp.NE.ECl) THEN
#endif
            DO js=ICUT2+1,nbin_aer
               jj=IQ(jesp,js)
               dqdt(jj)=dqdt(jj)+KERCD(jesp,js)*q(js)
            END DO
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         ENDIF
#endif
            ENDIF
         ENDIF
      END DO

      END
