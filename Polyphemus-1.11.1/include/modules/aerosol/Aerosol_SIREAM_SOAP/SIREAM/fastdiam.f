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

      SUBROUTINE FASTDIAM(nesp_aer,nbin_aer,is1,is2,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,VSW,MSW,DSW,RHOW,rho_dry)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine fastly computes dry and wet aerosol diameters  
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     nbin_aer: number of aerosol bins.
C     IS1 : first bin to compute.
C     IS2 : last bin to compute.
C     NEQ : number of equations.
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
C     
C     -- OUTPUT VARIABLES
C
C     VSW: wet volume aerosol concentration ([\mu.m^3.m^-3]).
C     MSW: single wet aerosol mass             ([\mu g]).
C     DSW: wet aerosol diameter             ([\mu.m]).
C     RHOW: wet aerosol density       ([\mu g.\mu.m^-3]).      
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2005: Edouard Debry, CEREA.
C     
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'dynaero.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'vara.inc'
      
      INTEGER is1,is2,neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq)

      INTEGER js,jj,jesp,iq(nesp_aer,nbin_aer)
      DOUBLE PRECISION qext(nesp_aer),vad,qn

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION DSW(nbin_aer),RHOW(nbin_aer)
      double precision rho_dry(nbin_aer)

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
     &        QT(js).GT.TINYM.AND.
     &        rho_dry(js).GT.TINYM ) THEN

            qn=q(js)            ! num aero conc (#part.m-3)
            RHOW(js)=rho_dry(js)       ! fixed aerosol density (µg.µm-3)

            vad=QT(js)/RHOW(js)     ! dry aerosol volume (µm3.m-3)
            DSD(js)=(vad/qn/cst_pi6)**cst_FRAC3 ! µm, cannot be less than tinyn

            CALL gerber_wet_diameter(RH, ! relative humidity (0<adim<1)
     &           TEMP,          ! temperature (Kelvin)
     &           DSD(js),      ! dry aero diam (m)
     &           DSW(js))      ! wet aero diam (m)

            VSW(js)=qn*cst_pi6*DSW(js)*DSW(js)*DSW(js) ! wet aerosol volume (µm3.m-3)
            qext(EH2O)=DMAX1(VSW(js)-vad,0.D0)*LMD(EH2O) ! water conc (µg.m-3)

            jj=IQ(EH2O,js)
            q(jj)=qext(EH2O)

            MSD(js)=QT(js)/q(js) ! single dry aero mass (µg)
            XSD(js)=DLOG(MSD(js)) ! dry log mass (adim)
            MSW(js)=MSD(js)+q(jj)/q(js) ! single wet mass (µg)

         ELSE
                                ! if too few aerosols or too few mass
                                ! we set variables of given bins as
                                ! its initial fixed ones, 
                                ! thus avoiding zero values
            RHOW(js)=RHOA      ! fixed aerosol density (g.m-3)
            MSD(js)=MSF(js)
            DSD(js)=DSF(js)
            XSD(js)=XSF(js)
            MSW(js)=MSF(js)
            DSW(js)=DSF(js)

         ENDIF
      END DO

      END
