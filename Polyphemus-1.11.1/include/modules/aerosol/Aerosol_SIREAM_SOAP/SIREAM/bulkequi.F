C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry and Kathleen Fahey
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

      SUBROUTINE BULKEQUI(neq,nesp_aer,nbin_aer,q,iq,QT,
     s     XSF,MSF,DSF,XSD,MSD,DSD,bin_density,psoap_config,psurrogate,
     s     ISOAPDYN,ionic,proton,lwc,liquid)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the gas/aerosol bulk equilibrium.
C     The partitioning among sections (bins) is made with the method
C     advocated by Spyros Pandis (flux ratio).
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ: number of aerosol variables.
C     nbin_aer: number of aerosol bins.
C     q: aerosol concentrations ([\mu g.m^-3]).
C     iq: index of aerosol species in q(*) vector.
C     QT: total aerosol concentration per bin ([\mu g.m^-3]).
C     XSF: neperian logarithm of fixed aerosol bin mass ([adim]).
C     MSF: fixed aerosol bin dry mass ([\mu g]).
C     DSF: fixed aerosol bin dry diameter ([\mu m]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q: the vector of aerosol variables.
C     XSD: neperian logarithm of moving aerosol bin mass ([adim]).
C     MSD: moving aerosol bin dry mass ([\mu g]).
C     DSD: moving aerosol bin dry diameter ([\mu m]).
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
C     2005/10: implement ability to not compute diameters, or
C     compute fastly with gerber (Edouard Debry, CEREA)
C     2017/09: adapt for implementation of SOAP (Youngseob Kim)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Edouard Debry and Kathleen Fahey, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varg.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER jj,js,jesp,iclip,ibegin,iter
      DOUBLE PRECISION aatot(nesp_aer),dq(nesp_aer),qold(nesp_aer)
      DOUBLE PRECISION frac(nesp_aer,nbin_aer)
      DOUBLE PRECISION aatoteq,qgasa,qgasi

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(nesp_aer,nbin_aer),RHOW(nbin_aer)

      DOUBLE PRECISION AA(nesp_aer,nbin_aer)

      DOUBLE PRECISION qaero(nesp_aer),qgas(nesp_aer)
      DOUBLE PRECISION organion, lwcorg, proton, lwc, ionic
      DOUBLE PRECISION liquid(12)

      double precision bin_density(nbin_aer)
      integer psoap_config,psurrogate
      integer ISOAPDYN
      double precision fractot(nbin_aer), bintot(nbin_aer), aerotot
      double precision qaerona, qaerocl, qgascl !! YK

C     ******zero init
      DO jesp=1,nesp_aer
         dq(jesp)=0.D0
         qold(jesp)=0.D0
         aatot(jesp)=TINYA
         qaero(jesp)=0.D0
         qgas(jesp)=0.D0
      END DO

      DO jesp=1,nesp_aer
         DO js=1,nbin_aer
            frac(jesp,js)=0.D0
         END DO
      END DO

C     ******compute local equilibrium

C     : compute local equi only if asked
      IF (ITHRM.EQ.0) THEN
         CALL EQPART(nesp_aer,nbin_aer,1,ICUT,neq,q,iq,QT,
     s        XSF,MSF,DSF,XSD,MSD,DSD,QINT,QGEQ,VSW,MSW,DSW,RHOW)
      ELSEIF (ITHRM.EQ.1) THEN
         CALL FASTDIAM(nesp_aer,nbin_aer,1,ICUT,neq,q,iq,QT,
     s        XSF,MSF,DSF,XSD,MSD,DSD,VSW,MSW,DSW,RHOW,bin_density)
      ENDIF

C     ******compute only c/e coefficients
!     compute total aa coef (s-1)
      DO jesp=1,nesp_aer        ! avoid sulfate
         IF (aerosol_species_interact(jesp).GT.0) THEN
            IF (jesp.NE.ESO4) THEN
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         IF (jesp.NE.ECl) THEN
#endif
            DO js=1,ICUT
               CALL COMPUTE_CONDENSATION_TRANSFER_RATE(DIFFG(jesp), ! diffusion coef (m2.s-1)
     $              VQMG(jesp), ! quadratic mean speed (m.s-1)
     $              STICK(jesp), ! accomadation coef (adim)
     $              DSW(js), ! wet aero diameter (µm)
     $              AA(jesp,js) ) ! c/e kernel coef (m3.s-1)
               aatot(jesp)= aatot(jesp)
     &              +AA(jesp,js)*q(js)
            END DO
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         ENDIF
#endif
         ENDIF
         ENDIF
      END DO

C     for sulfate
      IF (ISULFCOND.EQ.0) then
         jesp=ESO4
         DO js=1,nbin_aer
            CALL COMPUTE_CONDENSATION_TRANSFER_RATE(DIFFG(jesp), ! diffusion coef (m2.s-1)
     $           VQMG(jesp),    ! quadratic mean speed (m.s-1)
     $           STICK(jesp),   ! accomadation coef (adim)
     $           DSW(js),       ! wet aero diameter (µm)
     $           AA(jesp,js) )  ! c/e kernel coef (m3.s-1)
            aatot(jesp)= aatot(jesp)
     &           +AA(jesp,js)*q(js)
         END DO

C     for equilibrium bins
         aatoteq=0.d0
         DO js=1,ICUT
            aatoteq=aatoteq+aa(jesp,js)*q(js)
         END DO
      ELSE
         aatoteq=0.D0
      ENDIF

C     ******compute total aerosol mass
      DO jesp=E1,E2
         qaero(jesp)=0.D0
         DO js=1,ICUT
            jj=IQ(jesp,js)
            qaero(jesp)=qaero(jesp)+q(jj)
         END DO

         qgas(jesp) = q(IG(jesp))

!     store initial aero conc
         qold(jesp)=qaero(jesp)
      END DO

C     ******compute inorg gas/aero bulk equilibrium
      IF (ISULFCOND.EQ.0) then
         jesp=ESO4
         qgasi=qgas(jesp)
!     compute 'apparent' SO4 gas conc
!     i.e. the SO4 gas conc actually
!     seen by equilibrium aerosols
         qgasa= aatoteq/aatot(jesp) ! fraction < 1, adim
     &        *qgasi            ! gas conc, µg.m-3

!     give back into working array
         qgas(jesp)=qgasa
      ELSE
         jesp=ESO4
         qgasi = 0.D0
         qgasa = 0.D0
         qgas(jesp)=0.D0
      ENDIF

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      qaerona = qaero(ENa)
      qaerocl = qaero(ECl)
      qgascl = qgas(ECl)
      qaero(ENa) = 0.D0
      qaero(ECl) = 0.D0
      qgas(ECl) = 0.D0
#endif

      organion = 0.D0
      lwcorg = 0.D0
      proton = 0.D0
      ionic = 0.D0

c      DO iter = 1,NITER_AEC_AQ
         CALL isoropia_drv(nesp_aer,
     &        qaero,qgas,organion, lwcorg, ionic, proton, lwc, rh, temp,
     &        liquid)

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      qaero(ENa) = qaerona
      qaero(ECl) = qaerocl
      qgas(ECl) = qgascl
#endif 
!! YK as SCRAM

         if (ISOAPDYN.EQ.0) CALL soap_eq(nesp_aer, nbin_aer,
     &        neq, q, lwcorg, lwc, rh, ionic, proton, 
     &        temp, qaero, qgas, liquid, psoap_config,psurrogate)

c      ENDDO

c      qaero(EH2O) = lwc   + lwcorg

!     update gas conc µg.m-3
      IF (ISULFCOND.EQ.0) q(IG(ESO4))=qgasi-qgasa
!     give back initial SO4 gas conc
!     minus that consumed by equi bins
      q(IG(ENH3))=qgas(ENH3)
      q(IG(ENO3))=qgas(ENO3)
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
      q(IG(ECl)) =qgas(ECl)
#endif
      DO jesp=1,nesp_aec
         jj = aec_species(jesp)
         IF (aerosol_species_interact(jj).GT.0) q(IG(jj))=qgas(jj)
      ENDDO

      DO jesp=1,nesp_pankow
         jj = pankow_species(jesp)
         IF (aerosol_species_interact(jj).GT.0) q(IG(jj))=qgas(jj)
      ENDDO

      DO jesp=E1,E2             ! compute delta aero conc
         dq(jesp)=qaero(jesp)-qold(jesp)
      ENDDO

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      dq(ENa)=0.D0              ! do not treat sodium and chloride
      dq(ECl)=0.D0
#endif

C     ******redistribute on each size according to rates
      IF (ISULFCOND.EQ.0) then
         ibegin=-1              ! if sulfate equilibrium
      ELSE
         ibegin=ESO4            ! if sulfate dynamic, then avoid it
      ENDIF

      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
            IF (jesp.NE.ibegin) THEN
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         IF (jesp.NE.ECl) THEN
#endif
            iclip=0
            DO js=1,ICUT
                                ! fraction coef, adim

               frac(jesp,js)= AA(jesp,js)*q(js)
     &              /aatot(jesp)

               jj=IQ(jesp,js)
               q(jj)=q(jj)+dq(jesp)*frac(jesp,js)
               if (q(jj).le.0.d0) iclip=1
            END DO

C     warning : can lead to clipping for small sizes

            IF (iclip.EQ.1) THEN
               DO js=1,ICUT
                  jj=IQ(jesp,js)
                  q(jj)=qaero(jesp)*frac(jesp,js)
               END DO
            ENDIF
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         ENDIF
#endif
      ENDIF
      ENDIF
      END DO

    ! Water redistribution: assume the same distribution as the distribution of total aerosols in each bin. YK
      bintot = 0.D0
      fractot = 0.D0
      aerotot = 0.D0
      DO js=1,nbin_aer
         DO jesp=1,nesp_aer
            jj=IQ(jesp,js)
            bintot(js) = bintot(js) + q(jj)
            aerotot = aerotot + q(jj)
         END DO
      END DO

      DO js=1,nbin_aer
         if (aerotot .gt. 0.d0) then
            frac(EH2O,js) = bintot(js) / aerotot
            jj=IQ(EH2O,js)
            q(jj)=qaero(EH2O) * frac(EH2O,js)
         else
            q(jj) = 0.d0
         end if
      END DO

      END
