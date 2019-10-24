C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Karine Sartelet
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

      SUBROUTINE SULFDYN(neq,nesp_aer,nbin_aer,q,iq,
     &      QT,XSF,MSF,DSF,XSD,MSD,DSD, bin_density)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes equilibrium of small bins and
C     integrates sulfate for all bins
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     NBIN_AER: number of aerosol bins.
C     IQ: aerosol pointers.
C     QT: 1D total aerosol mass in each bin ([\mu.g.m^-3])
C     XSF: logarithm of fixed mean aerosol mass in bins. ([])
C     MSF: fixed mean aerosol mass in bins. ([\mu.g])
C     DSF: fixed mean aerosol diameter in bins. ([\mu.m])
C     XSD: logarithm of moving mean aerosol mass in bins. ([])
C     MSD: moving mean aerosol mass in bins. ([\mu.g])
C     DSD: moving mean aerosol diameter in bins. ([\mu.m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q: gas/aerosol concentrations ([\mu.g.m^-3]).
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
C     2005/10: Add ITHRM hook to compute diameter (Edouard Debry, CEREA)
C	  2013/11/27: Added bin_density (Stephanie Deschamps, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2005: Karine Sartelet, CEREA.
C
C------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'time.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER jesp,js,jj
      DOUBLE PRECISION aatot,dexploc

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(nesp_aer,nbin_aer),RHOW(nbin_aer)

      DOUBLE PRECISION AA_tmp(nbin_aer)
      double precision bin_density(nbin_aer)

C     Compute diameter with thermo or gerber formula
      IF (ITHRM.EQ.0) THEN
         CALL EQPART(nesp_aer,nbin_aer,1,nbin_aer,neq,q,iq,QT,
     s        XSF,MSF,DSF,XSD,MSD,DSD,
     s        QINT,QGEQ,VSW,MSW,DSW,RHOW)
      ELSEIF (ITHRM.EQ.1) THEN
         CALL FASTDIAM(nesp_aer,nbin_aer,1,nbin_aer,neq,q,iq,QT,
     s        XSF,MSF,DSF,XSD,MSD,DSD,VSW,MSW,DSW,RHOW, bin_density)

      ENDIF

C     Integrate sulfate now for all bins
      jesp=ESO4
      aatot=0.D0

      DO js=1,nbin_aer
         CALL COMPUTE_CONDENSATION_TRANSFER_RATE(DIFFG(jesp), ! diffusion coef (m2.s-1)
     $        VQMG(jesp),       ! quadratic mean speed (m.s-1)
     $        STICK(jesp),      ! accomadation coef (adim)
     $        DSW(js),          ! wet aero diameter (µm)
     $        AA_tmp(js) )      ! c/e kernel coef (m3.s-1)
         aatot = aatot + AA_tmp(js)*q(js)
      END DO

      If(aatot.GE.1.d-25) then
         dexploc=DEXP(-aatot*dt2)

         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            q(jj)=q(jj)+(AA_tmp(js)*q(js)/aatot)*(1.D0-dexploc)
     &           *q(IG(ESO4))
         END DO
         q(IG(ESO4))=DMAX1(q(IG(ESO4))*dexploc,0.d0)
      ENDIF

      END
