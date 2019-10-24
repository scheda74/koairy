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

      SUBROUTINE EQPART(nesp_aer,nbin_aer,is1,is2,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,QINT,QGEQ,VSW,MSW,DSW,RHOW)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the gas/aerosol equilibrium.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nbin_aer: number of aerosol bins.
C     IS1 : first bin to compute.
C     IS2 : last bin to compute.
C     NEQ : number of equations.
C     Q : aerosol concentration ([µg.m-3]).
C     iq: index of aerosol species in q(*) vector.
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
C     VSW: wet volume aerosol concentration ([\mu.m^3.m^-3]).
C     MSW: single wet aerosol mass             ([\mu g]).
C     DSW: wet aerosol diameter             ([\mu.m]).
C     RHOW: wet aerosol density       ([\mu g.\mu.m^-3]).
C     QINT : internal inorganic concentration ([\mu.g.m^-3]).
C     QGEQ : equilibrium gas concentration    ([\mu.g.m^-3]).
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
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'varp.inc'

      INTEGER is1,is2,neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq)

      INTEGER js,jj,jesp,iq(nesp_aer,nbin_aer)
      DOUBLE PRECISION qext(nesp_aer)

      DOUBLE PRECISION QGEQTMPA(nesp_aer),QINTMPA(NINTIS)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(nesp_aer,nbin_aer),RHOW(nbin_aer)

      DO js=is1,is2
                                ! zero init
         DO jesp=1,nesp_aer
            qext(jesp)=0.D0
         END DO

                                ! compute total dry mass
         QT(js)=0.D0
         DO jesp=E1,E2
            jj=IQ(jesp,js)
            qext(jesp)=q(jj)
            QT(js)=QT(js)+qext(jesp) ! µg.m-3
         END DO

         IF ( q(js).GT.TINYN.AND.
     &        QT(js).GT.TINYM ) THEN

            CALL STEP( nesp_aer,
     &           q(js),   ! num aero conc (#part.m-3)
     &           qext,          ! extern aero conc (µg.m-3)
     &           QT(js),        ! tot dry conc (µg.m-3)
     &           QINTMPA,       ! int inorg conc (µg.m-3)
     &           QGEQTMPA,      ! equi gas conc (µg.m-3)
     &           VSW(js),       ! vol aero conc (µm3.m-3)
     &           DSD(js),       ! dry aero diam (µm)
     &           DSW(js),       ! wet aero diam (µm)
     &           RHOW(js))      ! aerosol density (µg/µm-3)

            DO jj=1,nesp_aer
               QGEQ(jj,js) = QGEQTMPA(jj)
            ENDDO

            DO jj=1,NINTIS
               QINT(jj,js) = QINTMPA(jj)
            ENDDO

            jj=IQ(EH2O,js)
            q(jj)=qext(EH2O)
                                ! qn cannot be less than tinyn
            MSD(js)=QT(js)/q(js) ! single dry aero mass (µg)
            XSD(js)=DLOG(MSD(js)) ! dry log mass (adim)
            MSW(js)=MSD(js)+q(jj)/q(js) ! single wet mass (µg)

         ELSE
                                ! if too few aerosols or too few mass
                                ! we set variables of given bins as
                                ! its initial fixed ones,
                                ! thus avoiding zero values
            MSD(js)=MSF(js)
            DSD(js)=DSF(js)
            XSD(js)=XSF(js)
            MSW(js)=MSF(js)
            DSW(js)=DSF(js)
            RHOW(js)=RHOA
         ENDIF
      END DO

      END
