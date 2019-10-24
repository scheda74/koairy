C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University.  It is a
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

      subroutine vsrmchem(NGAS,NAER,NS,RHO_AERO,fixed_rho_aero,
     &     DBF_AERO,dsf_aero,xbf_aero,ZA,HUMID,press2,temp,lwc_c,
     &     t0,t1,rain_rate,pH,qscav)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves the aqueous-phase model.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     HUMID : specific humidity    ([kg.kg^{-1}]).
C     PRESS2: pressure             ([Pa]).
C     TEMP  : temperate            ([T]).
C     LWC_C : liquid water content ([g.m^{-3}]).
C     T0/T1 : initial/final time   ([s]).
C     rain_rate : rain precipitation rate ([mm/hr]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     ZA: gas-phase and aerosol concentration ([\mu.g/m^3]).
C
C     -- OUTPUT VARIABLES
C
C     PH : pH.
C     qscav : scavenged quantity by in-cloud scavenging
C     for dissolved aerosols species ([\mu.g/m^3])
C     and dissolves gas      species ([ppm])
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     1) Optimize conversion and define molar mass as parameters.
C     2) Introduce pointers for aerosols (E*) and for gas (ICTM*).
C     3) Clipping at the end (TINYAQ).
C     4) Rewrite mass balance for S and N.
C     5) Add inclusion of data_mass.inc and pointer_ctm.inc.
C     6) Remove call to VSRM (now included as such).
C     7) Define NITSUBAQ (number of subcycling timesteps = 10 here).
C     8) Remove include 'dropcom.inc'.
C     9) Move the initialization (before in AQOPERATOR.f).
C     10)Remove IAQ.
C     11)Clarify the computation of N2O5 and update the computation
C     of nitrate mass.
C     12)Remove the GOTO (with NITVSRM possible restart).
C     13)Replace NSECT by NS.
C     14)Remove tfin.
C     15)Conversion factors f* as input parameters (not computed) for
C     AQOPERATOR.F.
C     16)Remove call to DROPINIT.F (included as such).
C     17)Remove computations of DAER (replaced by DSF).
C     18)Update computation of total mass.
C     19)Add 3bis to avoid numerical difficulties for low SO2 (before:
C     in DECISIONS.F).
C     20)Comment call to DECISIONS.F
C     21)Uniform distribution of N2O5 in bins.
C     22)Remove include 'dropcom.inc'.
C     23) 2005/11/25: Implement in-cloud scavenging (Edouard Debry).
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Kathleen Fahey, CEREA, , on the basis of the VSRM model
C     (Carneggie Mellon University).
C     2005/10/3, cleaning and update, Bruno Sportisse, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'CONST_A.INC'
      include 'data_mass.inc'
      include 'pointer_ctm.inc'
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'pointer.inc'
      include 'num_aq.inc'

      INTEGER NGAS,NAER,NS

      double precision coefloc,tinit,coefloc2
      double precision dh2SO4,dso2
      double precision gas(ngas_aq), aerosol(NS,naers)
      double precision gasav(ngas_aq), aerosav(NS,naers)
      double precision rh,temp,press,lwc_c,t0,t1,initso2
      double precision dnew(NS),deltat,press2
      double precision numberconc(NS),totmass(NS)
      double precision ZA(NGAS+NS*NAER)
      double precision HUMID
      double precision nitbef,nitaf
      double precision sulfbef,sulfaf,sbal
      double precision fNH3,fHNO3,fHCl,fSO2
      double precision fHCHO,fO3,fOH,fHCOOH
      double precision fNO,fNO3,fNO2,fPAN,fHO2,fH2O2,fHNO2
      double precision pH
      double precision qscav(NGAS+NS*NAER),rain_rate
      double precision fdist(NS), fdist2(NS)

      double precision collision_eff,drop_diam,scav_coef
      double precision scav_frac,aqSO2,aqH2O2
      double precision DBF_AERO(NS+1),RHO_AERO(NS)
      double precision dsf_aero(NS),xbf_aero(NS+1)
      double precision fixed_rho_aero,totmasstot

      integer ifirstact
      integer isect,isp,i,istep,ii,nistep,ind_ok,j

      data collision_eff /0.9d0/

      integer nesp

      nesp = NGAS+NS*NAER

c     1) Initialisation
c     ------------------
      call initactiv(NS,dsf_aero,ifirstact)

      nistep = NITSUBAQ

      do i=1,ngas_aq
         gas(i)=0.d0
      enddo

      do isect=1,NS
         do isp=1,naers
            aerosol(isect,isp)=0.d0
         enddo
      enddo

      do isp=1,NESP
         qscav(isp)=0.d0
      enddo

      deltat = (t1-t0)/60.d0/NITSUBAQ !timestep in min

      press  = press2/101325.d0 ! pressure in atm

C     2) Compute sectional diameters and bimodal distribution
C     -------------------------------------------------------

C      call newdist(NS,DBF_AERO,dsf_aero,fdist,fdist2,ifirstact)

C     3) Conversion from \mug/m^3 to ppm
C     -----------------------------------

      coefloc = 8.314d-5*temp/press

      fNH3 = coefloc/mmNH3      !NH3
      fHNO3= coefloc/mmHNO3     !HNO3
      fHCl = coefloc/mmHCl      !HCl
      fSO2 = coefloc/mmSO2      !SO2
      fH2O2= coefloc/mmH2O2     !H2O2
      fHCHO= coefloc/mmHCHO     !HCHO
      fHNO2= coefloc/mmHNO2     !HNO2
      fO3  = coefloc/mmO3       !O3
      fOH  = coefloc/mmOH       !OH
      fHO2 = coefloc/mmHO2      !HO2
      fNO3 = coefloc/mmNO3      !NO3
      fNO  = coefloc/mmNO       !NO
      fNO2 = coefloc/mmNO2      !NO2
      fPAN = coefloc/mmPAN      !PAN
      fHCOOH = coefloc/mmHCOOH  !HCOOH

      coefloc2 = 8.314d0*temp/press2

      gas(nga)    = ZA(ictmNH3)  * coefloc2 / mmNH3  !NH3
      gas(ngn)    = ZA(ictmHNO3) * coefloc2 / mmHNO3 !HNO3
      gas(ngc)    = ZA(ictmHCl)  * coefloc2 / mmHCl  !HCl
      gas(ngso2)  = ZA(ictmSO2)  * coefloc2 / mmSO2  !SO2
      gas(ngh2o2) = ZA(ictmH2O2) * coefloc2 / mmH2O2 !H2O2
      gas(nghcho) = ZA(ictmHCHO) * coefloc2 / mmHCHO !HCHO
      gas(nghno2) = ZA(ictmHNO2) * coefloc2 / mmHNO2 !HNO2
      gas(ngo3)   = ZA(ictmO3)   * coefloc2 / mmO3   !O3
      gas(ngoh)   = ZA(ictmOH)   * coefloc2 / mmOH   !OH
      gas(ngho2)  = ZA(ictmHO2)  * coefloc2 / mmHO2  !HO2
      gas(ngno3)  = ZA(ictmNO3)  * coefloc2 / mmNO3  !NO3
      gas(ngno)   = ZA(ictmNO)   * coefloc2 / mmNO   !NO
      gas(ngno2)  = ZA(ictmNO2)  * coefloc2 / mmNO2  !NO2
      gas(ngpan)  = ZA(ictmPAN)  * coefloc2 / mmPAN  !PAN

      gas(nghcooh)   = 0.1d0*gas(ngh2o2) !HCOOH
      gas(ngch3o2h)  = 0.2d0*gas(ngh2o2) !CH3OOH(g)  ppm = 0.2*H2O2
      gas(ngch3o2)   = 1.0d-6            !CH3O2(g)   ppm
      gas(ngch3oh)   = 1.0d-3            !CH3OH(g)   ppm = 1 ppb
      gas(ngch3co3h) = 0.05d0*gas(ngh2o2)!CH3C(O)OOH ppm = 0.05*H2O2

      initso2 = gas(ngso2)

      DO i = 1,NS
         ii = NGAS + i
         aerosol(i,nas) = ZA(ii+(ENa-1)*NS) ! Na
         aerosol(i,na4) = ZA(ii+(ESO4-1)*NS) ! SO4
         aerosol(i,naa) = ZA(ii+(ENH3-1)*NS) ! NH4
         aerosol(i,nan) = ZA(ii+(ENO3-1)*NS) ! NO3
         aerosol(i,nac) = ZA(ii+(ECl-1)*NS) ! Cl
         aerosol(i,nar) = ZA(ii+(EMD-1)*NS) ! DUST
         aerosol(i,nae) = ZA(ii+(EBC-1)*NS) ! EC
         aerosol(i,nao) = ZA(ii+(EPOA-1)*NS) ! POA
         aerosol(i,naw) = ZA(ii+(EH2O-1)*NS) ! H2O
      ENDDO

C     3bis) At low SO2 concentration, transfer all SO2 to
C     SO4-- in order to avoid numerical difficulties
C     ---------------------------------------------------
      IF (gas(ngso2).LE.minso2) THEN
         IF (gas(ngh2o2).GE.gas(ngso2)) THEN
            gas(ngh2o2) = gas(ngh2o2) - gas(ngso2)
         ELSE
            gas(ngh2o2) = 0.d0
         ENDIF
         dso2 = (ZA(ictmSO2)/(NS-IFIRSTACT+1))*mmSO4/mmSO2
         gas(ngso2) = 0.d0

         DO i=IFIRSTACT,NS
            aerosol(i,na4) = aerosol(i,na4)+dso2
         ENDDO
      ENDIF

C     4) Transfer all H2SO4(g) to aerosol phase
C     ------------------------------------------

      dh2SO4 = (ZA(ictmH2SO4)/NS)*mmSO4/mmH2SO4
      DO i=1,NS
         aerosol(i,na4) = aerosol(i,na4)+dh2SO4
      ENDDO
      ZA(ictmH2SO4) = 0.d0

C     5) Save initial concentrations (for restart)

C     6) Save initial concentrations (for restart)
C     ----------------------------------------------

      do isect=1,NS
         do isp=1,naers
            aerosav(isect,isp) = aerosol(isect,isp)
         enddo
      enddo

      do i=1,ngas_aq
         gasav(i) = gas(i)
      enddo

C     6) CALCULATION OF TOTAL SULFUR/NITROGEN  MASS BEFORE THE CALL
C     --------------------------------------------------------------

      sulfbef = gas(ngso2)/coefloc
      nitbef  = (gas(nga)+gas(ngn)+gas(ngpan)+gas(ngno3)
     &     +gas(ngno)+gas(ngno2)+gas(nghno2))/coefloc

      do i=1, NS
         nitbef  = nitbef  + aerosol(i,naa)   /mmNH4
     &        + aerosol(i,nan)   /mmNO3
         sulfbef = sulfbef + aerosol(i,na4)   /mmSO4
     &        + aerosol(i,nahso5)/mmHSO5
     &        + aerosol(i,nahmsa)/mmHMSA
      enddo
      nitbef = nitbef  * mmN
      sulfbef= sulfbef * mmS

C     7) Calculation of aerosol number
C     --------------------------------

      do i = 1,NS
         totmass(i)= aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)
         numberconc(i) = totmass(i)/(dsf_aero(i)**3.d0)
     &        /cst_pi6/RHO_AERO(i)
      enddo


C     7bis) Calculation of the bulk distribution factors fdist
      totmasstot = 0.D0
      do i=ifirstact,NS
         totmasstot = totmass(i) + totmasstot
      enddo
      do i=1,ifirstact-1
         fdist(i) = 0.d0
      enddo
      do i=ifirstact,NS
         fdist(i) = totmass(i)/totmasstot
      enddo

C     8) Compute relative humidity
C     ----------------------------

      call COMPUTE_RELATIVE_HUMIDITY(HUMID,temp,press2,rh)
      rh = DMIN1(DMAX1(rh, Threshold_RH_inf), Threshold_RH_sup)

C     ***************************************************
C     Begin computations
C     ***************************************************

      ind_ok = 0

      DO j=1,NITVSRM
         IF (ind_ok.eq.0) THEN

c     10) Simulate cloud/fog
c     ----------------------

            do istep = 1, nistep
               tinit= (istep-1)*deltat

               call aqoperator(NS, tinit, deltat, gas, aerosol,
     &              lwc_c, temp, press,
     &              fHCHO,fhcooh,fso2,fh2o2,fNH3,fhno3,fhcl,
     &              ifirstact,
     &              pH,
     &              aqSO2,aqH2O2,
     &              fdist)

            enddo

c     11) Compute TOTAL SULFUR/NITROGEN MASS AFTER THE CALL
C     --------------------------------------------------------------

            sulfaf = gas(ngso2)/coefloc
            nitaf  = (gas(nga)+gas(ngn)+gas(ngpan)+gas(ngno3)
     &           +gas(ngno)+gas(ngno2)+gas(nghno2))/coefloc

            do i=1, NS
               nitaf  = nitaf  + aerosol(i,naa)   /mmNH4
     &              + aerosol(i,nan)   /mmNO3
               sulfaf = sulfaf + aerosol(i,na4)   /mmSO4
     &              + aerosol(i,nahso5)/mmHSO5
     &              + aerosol(i,nahmsa)/mmHMSA
            enddo
            nitaf = nitaf  * mmN
            sulfaf= sulfaf * mmS

            sbal    = DABS(sulfaf/sulfbef-1.D0)

C     12) Check if mass balance is OK
C     Restart otherwise (from step 9)
C     with a timestep divided by 2.
C     -------------------------------

            if ((sulfbef.gt.0.1d0) .and.
     &           (sbal.gt.RTOLSULF) .and.
     &           (initso2.gt.minso2)) then
C               write(*,*) 'VSRM: Sulfate balance not OK: restart'

C     computation has failed
C     ----------------------
               do isect=1,NS
                  do isp=1,naers
                     aerosol(isect,isp) = aerosav(isect,isp)
                  enddo
               enddo

               do i=1,ngas_aq
                  gas(i) = gasav(i)
               enddo

               deltat = deltat/2.d0
               nistep = 2*nistep

            else

C     computation is OK
C     -----------------
               ind_ok = 1
            endif

         ENDIF
      ENDDO

C     13) Compute new distribution by projection to
C     the aerosol population
C     ----------------------------------------------

      do i = 1,NS
         totmass(i)=aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)
         dnew(i)=(totmass(i)/numberconc(i)
     &        /rho_aero(i)/cst_pi6)**cst_FRAC3

      enddo

      call redistaq(NS,DBF_AERO,xbf_aero,fixed_rho_aero,dnew,aerosol)

C     14) Update gas & aerosol concentrations
C     ---------------------------------------

      ZA(ictmNH3) = gas(nga)   /coefloc2 * mmNH3          !NH3
      ZA(ictmHNO3)= gas(ngn)   /coefloc2 * mmHNO3         !HNO3
      ZA(ictmHCl) = gas(ngc)   /coefloc2 * mmHCl          !HCl
      ZA(ictmSO2) = gas(ngso2) /coefloc2 * mmSO2          !SO2
      ZA(ictmH2O2)= gas(ngh2o2)/coefloc2 * mmH2O2         !H2O2
      ZA(ictmHCHO)= gas(nghcho)/coefloc2 * mmHCHO         !HCHO
      ZA(ictmHNO2)= gas(nghno2)/coefloc2 * mmHNO2         !HNO2
      ZA(ictmO3)  = gas(ngo3)  /coefloc2 * mmO3           !O3
      ZA(ictmOH)  = gas(ngoh)  /coefloc2 * mmOH           !OH
      ZA(ictmHO2) = gas(ngho2) /coefloc2 * mmHO2          !HO2
      ZA(ictmNO3) = gas(ngno3) /coefloc2 * mmNO3          !NO3
      ZA(ictmNO)  = gas(ngno)  /coefloc2 * mmNO           !NO
      ZA(ictmNO2) = gas(ngno2) /coefloc2 * mmNO2          !NO2
      ZA(ictmPAN) = gas(ngpan) /coefloc2 * mmPAN          !PAN

      DO i = 1,NS
         ii = NGAS + i
         ZA(ii+(ENa-1)*NS) = aerosol(i,nas) ! Na
         ZA(ii+(ESO4-1)*NS)= aerosol(i,na4) ! SO4
         ZA(ii+(ENH3-1)*NS)= aerosol(i,naa) ! NH4
         ZA(ii+(ENO3-1)*NS)= aerosol(i,nan) ! NO3
         ZA(ii+(ECl-1)*NS) = aerosol(i,nac) ! Cl
         ZA(ii+(EMD-1)*NS) = aerosol(i,nar) ! DUST
         ZA(ii+(EBC-1)*NS) = aerosol(i,nae) ! EC
         ZA(ii+(EPOA-1)*NS)= aerosol(i,nao) ! POA
         ZA(ii+(EH2O-1)*NS)= aerosol(i,naw) ! H2O
      ENDDO
C     15) Compute in-cloud scavenging
C     -------------------------------

      IF (rain_rate.GT.0.d0) THEN
                                ! compute scavenging coefficient, same for all species
         drop_diam = 9.76d-4 * (rain_rate**0.21d0 ) ! m
         scav_coef = 4.17d-7 * collision_eff * rain_rate / drop_diam ! s-1

                                ! compute scavenged fraction of droplet
         scav_frac = 1.D0 - dexp(scav_coef * (t0 - t1) ) ! adim

                                ! compute scavenged quantiy
         DO i = ifirstact,NS
            ii = NGAS + i
            DO isp = EMD,EPOA
               qscav(ii+(isp-1)*NS) = ZA(ii+(isp-1)*NS) * scav_frac
            ENDDO
         ENDDO

         qscav(ictmSO2) = aqSO2 * scav_frac * 0.7901d0
         qscav(ictmH2O2)= aqH2O2 * scav_frac

                                ! remove scavenged quantity
         DO i = ifirstact,NS
            ii = NGAS + i
            DO isp = EMD,EPOA
               ZA(ii+(isp-1)*NS) = ZA(ii+(isp-1)*NS)
     &              - qscav(ii+(isp-1)*NS)
            ENDDO
         ENDDO

         ZA(ictmSO2) = ZA(ictmSO2) - qscav(ictmSO2)
         ZA(ictmH2O2)= ZA(ictmH2O2) - qscav(ictmH2O2)
      ENDIF

C     16)
C----
      DO i=1,NESP
         ZA(i)=dmax1(ZA(i),TINYAQ)
      ENDDO

      return
      end
