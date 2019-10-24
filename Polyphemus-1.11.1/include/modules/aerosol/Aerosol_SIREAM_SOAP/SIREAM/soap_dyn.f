C-----------------------------------------------------------------------
C     Copyright (C) 2017, CEREA - ENPC - EDF R&D
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
C     component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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

      SUBROUTINE SOAP_DYN(nesp_aer, nbin_aer, rh,
     &     ionic, proton, lwc, liquid,
     &     temp, psoap_config, psurrogate,deltat,
     &     DSD, neq, q, iq)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes ...
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     FLAG: whether to solved hydrophilic (=0) or hydrophobic (=1) species.
C     PROTON: hydronium ion concentration ([\mu g.m^-3]).
C     LWC: total liquid water content ([\mu g.m^-3]).7
C     RH: relative humidity 0< <1 ([]).
C     TEMP: temperature ([Kelvin]).
C     IOLIGO: flag for oligomerization (true if =1)
C
C     -- INPUT/OUTPUT VARIABLES
C
C     AERO: aerosol bulk concentration ([\mu g.m^-3]).
C     GAS: gas concentration ([\mu g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C     ORGANION: organic ions ([\mu mol.m^-3]).
C     LWCORG: organic liquid water content ([\mu g.m^-3]).
C     CHP: hydronium ion concentration in water ([mol.L^-1]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2017: Youngseob Kim, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varq.inc'
     
      INTEGER nesp_aer,nbin_aer,neq
      double precision q(neq)
      INTEGER iq(nesp_aer,nbin_aer)
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
      DOUBLE PRECISION liquid(12)

      INTEGER jj,jesp,js
      double precision qaero(nesp_aer), qgas(nesp_aer)
      integer psoap_config, psurrogate
      double precision deltat
      double precision DSD(nbin_aer)
      double precision csol(nbin_aer) !! Concentration of solid particles (PBC + PMD)

      lwcorg = 0.D0

C     ******compute total aerosol mass
      DO jesp=E1,E2
         qaero(jesp)=0.D0
C         DO js=(ICUT2+1),nbin_aer
         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            qaero(jesp)=qaero(jesp)+q(jj)
         END DO
         qgas(jesp) = q(IG(jesp))
      enddo

      if(lwc.GT.1.d-19) then
        chp = proton / lwc * 1.0e3
      else
        chp = 1.d-19
        lwc = 1.d-19
      endif

      do js=1,nbin_aer
         csol(js) = q(IQ(EMD,js)) + q(IQ(EBC,js))
      enddo

      CALL soap_main(lwc, rh, temp, ionic, chp, lwcorg,
     &     psoap_config, psurrogate,deltat,DSD,csol,liquid,
     &     nesp_aer, neq, q, qaero, qgas)

      END
