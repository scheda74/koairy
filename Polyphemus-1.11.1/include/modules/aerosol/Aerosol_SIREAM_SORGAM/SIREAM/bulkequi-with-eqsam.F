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

      SUBROUTINE BULKEQUI(neq,nbin_aer,q,iq,QT,
     s            XSF,MSF,DSF,XSD,MSD,DSD)

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
      INCLUDE 'pointer.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'emw.inc'
      INCLUDE 'stick.inc'
      INCLUDE 'eqsam.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(NEXT,nbin_aer)

      INTEGER jj,js,jesp,iclip,ibegin
      DOUBLE PRECISION aatot(NEXT),qext2(NEXT)
      DOUBLE PRECISION qtot2(NEXT),dq(NEXT)
      DOUBLE PRECISION wi(NEINORG),w(NEINORG)
      DOUBLE PRECISION aerliq(NILIQ),aersld(NISLD)
      DOUBLE PRECISION gas(3),cntrl(2),other(6)
      DOUBLE PRECISION frac(NEXT,nbin_aer),qold2(NEXT)
      DOUBLE PRECISION aatoteq,qgasa,qgasi
      CHARACTER*15 scase

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(NEXT,nbin_aer),RHOW(nbin_aer)

      DOUBLE PRECISION AA(NEXT,nbin_aer)

C     For eqsam
      real eq_yi(eq_imax,eq_nca)
      real eq_yo(eq_imax,eq_nco)
      INTEGER iopt


C     ******zero init
      DO jesp=1,NEXT
         dq(jesp)=0.D0
         qtot2(jesp)=0.D0
         qext2(jesp)=0.D0
         qold2(jesp)=0.D0
         aatot(jesp)=TINYA
      END DO

      DO jesp=1,NEXT
         DO js=1,nbin_aer
            frac(jesp,js)=0.D0
         END DO
      END DO

C     ****** if sulfate computed dynamically (sulfdyn.f) avoid it
      IF (ISULFCOND.EQ.0) THEN
         ibegin=G1
      ELSE
         ibegin=ENH3
      ENDIF


C     : compute local equi only if asked
      IF (ITHRM.EQ.0) THEN
         CALL EQPART(nbin_aer,1,ICUT,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,QINT,QGEQ,VSW,MSW,DSW,RHOW)
      ELSEIF (ITHRM.EQ.1) THEN
         CALL FASTDIAM(nbin_aer,1,ICUT,neq,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD,VSW,MSW,DSW,RHOW)
      ENDIF

C     ******compute only c/e coefficients
                                ! compute total aa coef (s-1)
      DO jesp=ENH3,G2           ! avoid sulfate



            DO js=1,ICUT
               CALL COMPUTE_CONDENSATION_TRANSFER_RATE(DIFFG(jesp), ! diffusion coef (m2.s-1)
     $              VQMG(jesp), ! quadratic mean speed (m.s-1)
     $              STICK(jesp), ! accomadation coef (adim)
     $              DSW(js),    ! wet aero diameter (µm)
     $              AA(jesp,js) ) ! c/e kernel coef (m3.s-1)
               aatot(jesp)= aatot(jesp)
     &              +AA(jesp,js)*q(js)
            END DO



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
      ENDIF

C     ******compute total aerosol mass
      DO jesp=E1,E2
         DO js=1,ICUT
            jj=IQ(jesp,js)
            qext2(jesp)=qext2(jesp)+q(jj)
         END DO

                                ! store initial aero conc
         qold2(jesp)=qext2(jesp)
      END DO

C     ******compute inorg gas/aero bulk equilibrium
      IF (ISULFCOND.EQ.0) then
                             ! save initial real SO4 gas conc
         qgasi=q(IG(ESO4))

                             ! compute 'apparent' SO4 gas conc
                             ! i.e. the SO4 gas conc actually
                             ! seen by equilibrium aerosols
         qgasa= aatoteq/aatot(ESO4) ! fraction < 1, adim
     &        *qgasi         ! gas conc, µg.m-3

                                ! give back into working array
         q(IG(ESO4))=qgasa
      endif

      DO jesp=ENa,ECl
         qtot2(jesp)=qext2(jesp)+q(IG(jesp))
      ENDDO






      IF (ISULFCOND.EQ.1) qtot2(ESO4)=qext2(ESO4)

      DO jesp=1,ECl-ENa+1
         wi(jesp)= qtot2(ENa-1+jesp) ! µg.m-3
     &        /EMW(ENa-1+jesp)  ! µg.mol-1
      END DO


      IF (MTHRM==0) then !isoropia
         cntrl(1)=0.D0          ! forward mode
         cntrl(2)=MTSBL         ! metastable option

                                ! convert µg to moles
         CALL ISOROPIA(wi,RH,TEMP,cntrl,w,gas,
     &        aerliq,aersld,scase,other)

                                ! update gas conc µg.m-3
         q(IG(ENH3))=gas(1)*EMW(ENH3)
         q(IG(ENO3))=gas(2)*EMW(ENO3)

         q(IG(ECl)) = gas(3)*EMW(ECl)

      ELSE ! eqsam
         iopt=1

         eq_yi(1,1)=TEMP
         eq_yi(1,2)=RH
         eq_yi(1,11)=PRES*1D-3 ! Pa->hPa
                                ! convert µg to moles
C     convert µg to µmol
         eq_yi(1,3)= wi(3)      !NH3  (g) + NH4+  (p)
     $        * 1d6             ! µmol.m^3
         eq_yi(1,4)=wi(2)       ! H2SO4    + SO4-- (p)
     $        * 1d6             ! µmol.m^3
         eq_yi(1,5)=wi(4)       ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
     $        * 1d6             ! µmol.m^3
         eq_yi(1,6)=wi(1)       ! Na+ (ss  + xsod) (a)   [umol/m^3 air]
     $        * 1d6             ! µmol.m^3
         eq_yi(1,7)=wi(5)       ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
     $        * 1d6             ! µmol.m^3
c         eq_yi(1,8)=5D-10 * q(EMD) ! K+   (p) from Dust
c     $        * 1d6             ! µmol.m^3
c         eq_yi(1,9)=5D-10 * q(EMD) ! K+   (p) from Dust
c     $        * 1d6             ! µmol.m^3
c         eq_yi(1,10)=5D-10 * q(EMD) ! Mg++ (p) from Dust
c     $        * 1d6             ! µmol.m^3
         eq_yi(1,8)= 2.D-03
         eq_yi(1,9)= 2.D-03
         eq_yi(1,10)= 2.D-03

        CALL eqsam_v03d(eq_yi,eq_yo,eq_nca,eq_nco,
     $        iopt,eq_loop,eq_imax)
                                ! update gas conc µg.m-3
         q(IG(ENH3))=eq_yo(1,10)*EMW(ENH3)*1d-6
         q(IG(ENO3))=eq_yo(1,9)*EMW(ENO3)*1d-6
         q(IG(ECl)) =eq_yo(1,8)*EMW(ECl)*1d-6
      ENDIF
      IF (ISULFCOND.EQ.0) q(IG(ESO4))=0.D0    ! all SO4 condense

                                ! update aero conc
                                ! compute delta aero conc
      DO jesp=ibegin,ECl        ! avoid sodium
         qext2(jesp)=qtot2(jesp)-q(IG(jesp))
         dq(jesp)=qext2(jesp)-qold2(jesp)
      END DO





                                ! give back initial SO4 gas conc
                                ! minus that consumed by equi bins
      IF (ISULFCOND.EQ.0) q(IG(ESO4))=qgasi-qgasa

C     ******compute org gas/aero bulk equilibrium
      DO jesp=EARO1,EBiBmP
         qtot2(jesp)=qext2(jesp)+q(IG(jesp)) ! µg.m-3
      END DO
      qtot2(EPOA)=qext2(EPOA)

                                ! solve non-linear system
                                ! qtot2(*) is gas+aerosol
      CALL EQORG2(qtot2,qext2)
                                ! qext2(*) is aero equi conc

                                ! update gas conc µg.m-3
                                ! compute delta aero conc
      DO jesp=EARO1,EBiBmP
         q(IG(jesp))=qtot2(jesp)-qext2(jesp)
         dq(jesp)=qext2(jesp)-qold2(jesp)
      END DO

C     ******redistribute on each size according to rates
      DO jesp=ibegin,G2



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
                  q(jj)=qext2(jesp)*frac(jesp,js)
               END DO
            ENDIF



      END DO

      END
