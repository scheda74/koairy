C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
C     Author(s): Maryline Tombette, Yelva Roustan
C
C     This file is part of the Simple Aqueous model (SIMPLE_AQUEOUS), a
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

      subroutine simple_aqueous_module(ns, ns_aer, nbin_aer,
     $     density_aer, fixed_density_aer, DBF_AERO, dsf_aero, xbf_aero,
     $     ZA, press_in_pa, temp, lwc_c, t0, t1, rain_rate, pH, qscav)

C-----------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves the aqueous-phase model.
C
C-----------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     ns : number of gas phase species in the chemical mechanism.
C     ns_aer : number of particle phase species in the chemical mechanism.
C     nbin_aer : number of size section in the aerosol model.
C     density_aer : size variable aerosol density ([microg / microm^3]).
C     fixed_density_aer: fixed aerosol density ([microg / microm^3]).
C     press_in_pa : pressure ([Pa]).
C     temp : temperate ([K]).
C     LWC_C : liquid water content ([g.m^{-3}]).
C     T0, T1 : initial and final time ([s]).
C     rain_rate : rain precipitation rate ([mm/hr]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     ZA: gas-phase and aerosol concentration ([\mu.g/m^3]).
C
C     -- OUTPUT VARIABLES
C
C     PH : pH.
C     qscav : scavenged quantity by in-cloud scavenging.
C     for dissolved aerosols species ([\mu.g/m^3]).
C     and dissolves gas species ([ppm]).
C
C-----------------------------------------------------------------------
C
C     -- REMARKS
C
C     Implemented on the basis of Pandis/Seinfeld book.
C
C-----------------------------------------------------------------------
C
C     -- MODIFICATIONS
C     Modified to improve the numerical behaviour and to correct several
C     errors (Yelva Roustan, 2015).
C
C-----------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Marilyne Tombette, CEREA.
C
C-----------------------------------------------------------------------

      implicit none

      include 'CONST_A.INC'
      include 'pointer_ctm.inc'
      include 'pointer.inc'
      include 'data_mass.inc'
      include 'aerpar_simple.inc'
      include 'droppar_simple.inc'
      include 'num_aq.inc'

      integer ns, ns_aer, nbin_aer

      double precision coef_cte
      double precision gas(ns_aq), aerosol(nbin_aer, naers)
      double precision henrys(ns_aq)
      double precision cst_dissoc(Nreact_dissoc)
      double precision cst_oxydation(Nreact_oxydation)
      double precision temp, press_in_atm, lwc_c, lwc, t0, t1
      double precision dnew(nbin_aer), deltat, deltatmax, press_in_pa
      double precision numberconc(nbin_aer), totmass(nbin_aer)
      double precision ZA(ns + nbin_aer * ns_aer)
      double precision pH
      double precision qscav(ns + nbin_aer * ns_aer), rain_rate
      double precision fdist(nbin_aer)

C     Concentration needed to compute the pH.
      double precision conc_in(5)

C     Concentration updated with the pH.
      double precision conc_out(10)

      double precision collision_eff, drop_diam, scav_coef
      double precision so2g, henrys_eff_so2
      double precision scav_frac, aqSO2, aqH2O2, aqO3
      double precision DBF_AERO(nbin_aer + 1), density_aer(nbin_aer)
      double precision dsf_aero(nbin_aer), xbf_aero(nbin_aer + 1)
      double precision fixed_density_aer, totmasstot

      double precision dsivdt_o3, dsivdt_h2o2, dsivdt
      double precision fc_o3, fc_h2o2, fc_siv
      double precision tot_aer(naers), aa
      double precision sulfate, nitrate, ammonium
      double precision nitrate_init, ammonium_init, sulfate_init
      double precision report, fo3, fh2o2
      integer ifirstact
      integer isect, isp, i, istep, ii, nistep

      data collision_eff /0.9d0/

      double precision cst_pr, cst_pr2, cst_RT, cst_RTLWC
      integer Nesp

C     Initialisation.
C     ---------------

      cst_pr = 1.987d-3         ! perfect gas constant in kcal.K-1.mol-1
      cst_pr2 = 8.206d-2        ! perfect gas constant in atm.L_air.K-1.mol-1
      cst_RT = cst_pr2 * temp   ! in atm.L_air.mol-1

      lwc = lwc_c * 1.d-6       ! in L_water/L_air
      cst_RTLWC = cst_pr2 * temp * lwc ! in atm.mol-1.L_water

      Nesp = ns + nbin_aer * ns_aer

      press_in_atm  = press_in_pa / 101325.d0 ! pressure in atm

      call initactiv(nbin_aer, dsf_aero, ifirstact)

      nistep = NITSUBAQ
      deltat = (t1 - t0) / nistep !timestep in min

      do i = 1, ns_aq
         gas(i) = 0.d0
      enddo

      do isect = 1, nbin_aer
         totmass(isect) = 0.d0
      enddo

      do isp = 1, Nesp
         qscav(isp) = 0.d0
      enddo

      do isect = 1, nbin_aer
         do isp = 1, naers
            aerosol(isect, isp) = 0.d0
         enddo
      enddo

C     Compute henry's and equilibrium constants.
C     ------------------------------------------

      coef_cte= (1.d0/temp - 1.d0/temp_ref)

      do i = 1, ns_aq
         henrys(i) = chenry(i) * exp(-dhhenry(i) * coef_cte / cst_pr)
      enddo

      do i = 1, nreact_dissoc
         cst_dissoc(i) = ckdissoc(i)
     $        * exp(dhkdissoc(i) * coef_cte)
      enddo

      do i = 1, nreact_oxydation
         cst_oxydation(i) = ckoxydation(i)
     $        * exp(-dhkoxydation(i) *coef_cte
     $        / cst_pr)
      enddo

C     Initialize local variable.
C     --------------------------

C     Gas phase concentrations with conversion \mu g.m-3 -> atm.
      gas(igso2)  = ZA(ictmSO2) * 1.d-9 * cst_RT / mmSO2  !SO2
      gas(ignh3)  = ZA(ictmNH3)  * 1.d-9 * cst_RT / mmNH3  !NH3
      gas(ighno3) = ZA(ictmHNO3) * 1.d-9 * cst_RT / mmHNO3 !HNO3
      gas(igh2o2) = ZA(ictmH2O2) * 1.d-9 * cst_RT / mmH2O2 !H2O2
      gas(igo3)   = ZA(ictmO3)   * 1.d-9 * cst_RT / mmO3   !O3
      gas(igCO2) = 350.d0 * 1.d-6 * press_in_atm ! 350 ppm in averaged.

C     Particle phase concentration in \mu g.m-3.
      do i = 1, nbin_aer
         ii = ns + i
         aerosol(i, nas) = ZA(ii + (ENa-1)  * nbin_aer) ! Na
         aerosol(i, na4) = ZA(ii + (ESO4-1) * nbin_aer) ! SO4
         aerosol(i, naa) = ZA(ii + (ENH3-1) * nbin_aer) ! NH4
         aerosol(i, nan) = ZA(ii + (ENO3-1) * nbin_aer) ! NO3
         aerosol(i, nac) = ZA(ii + (ECl-1)  * nbin_aer) ! Cl
         aerosol(i, nar) = ZA(ii + (EMD-1)  * nbin_aer) ! DUST
         aerosol(i, nae) = ZA(ii + (EBC-1)  * nbin_aer) ! EC
         aerosol(i, nao) = ZA(ii + (EPOA-1) * nbin_aer) ! POA
         aerosol(i, naw) = ZA(ii + (EH2O-1) * nbin_aer) ! H2O
      enddo

C     Transfer all H2SO4(g) to particle phase.
C     ----------------------------------------
C     Could be improved to be proportional to the particle surface.
      do i = 1, nbin_aer
         aerosol(i, na4) = aerosol(i, na4) +
     $        ZA(ictmH2SO4) / nbin_aer * mmSO4 / mmH2SO4
      enddo
      ZA(ictmH2SO4) = 0.d0

C     Compute total mass and number concentration for each size section.
C     ------------------------------------------------------------------
      do i = 1, nbin_aer
         totmass(i) = aerosol(i,naa) + aerosol(i,na4) +
     $        aerosol(i,nan) + aerosol(i,nac) + aerosol(i,nas) +
     $        aerosol(i,nao) + aerosol(i,nae) + aerosol(i,nar)
         numberconc(i) = totmass(i) / (dsf_aero(i)**3.d0)
     $        / cst_pi6 / density_aer(i)
      enddo

C     Calculation of the bulk distribution factors fdist.
C     ---------------------------------------------------
      totmasstot = 0.D0
      do i = ifirstact, nbin_aer
         totmasstot = totmass(i) + totmasstot
      enddo

      do i = 1, ifirstact - 1
         fdist(i) = 0.d0
      enddo

      do i = ifirstact, nbin_aer
         fdist(i) = totmass(i) / totmasstot
      enddo

C     Compute total sulfate, ammonium and nitrate.
C     --------------------------------------------
      sulfate = 0.d0
      nitrate = 0.d0
      ammonium = 0.d0

      do i = ifirstact, nbin_aer
         sulfate = sulfate + aerosol(i, na4) ! in \mu g.m-3
         nitrate = nitrate + aerosol(i, nan) ! in \mu g.m-3
         ammonium = ammonium + aerosol(i, naa) ! in \mu g.m-3
      enddo

C     Save the initial values of total sulfate, ammonium and nitrate.
      sulfate_init = sulfate
      nitrate_init = nitrate
      ammonium_init = ammonium

C     Initialize the cloud/fog conditions.
C     ------------------------------------
      conc_in(1) = gas(igso2)   ! S(IV) in atm
      conc_in(2) = sulfate * 1.d-9 * cst_RT / mmSO4 ! S(VI) in atm
      conc_in(3) = ammonium * 1.d-9 * cst_RT / mmNH4 + gas(ignh3) ! in atm
      conc_in(4) = nitrate * 1.d-9 * cst_RT / mmNO3 + gas(ighno3) ! in atm
      conc_in(5) = gas(igco2)   ! in atm

C     First computation of pH.
      call compute_ph(conc_in, henrys, cst_dissoc, temp, lwc, ph,
     $     conc_out)

C     First computation of O3 and H2O2 aqueous concentration.
C     An instantaneous equilibrium between the gas phase concentration
C     and the aqueous phase concentration is assumed.
            aqo3 = gas(igo3) * henrys(igo3) /
     $           (1.d0 + cst_RTLWC * henrys(igo3)) ! in mol.L-1_water

            aqh2o2 = gas(igh2o2) * henrys(igh2o2) /
     $           (1.d0 + cst_RTLWC * henrys(igh2o2)) ! in mol.L-1_water

C     Time integration with subcycling time step.
C     -------------------------------------------
      do istep = 1, nistep

C     No aqueous chemistry without S(IV).
         if (conc_in(1).gt.0.d0) then

C     Compute reaction rates of S(IV) with O3 and H2O2.
C     -------------------------------------------------

C     Oxidation by O3 with the reaction rate from Seinfel & Pandis.
            dsivdt_o3 = - (cst_oxydation(1) * conc_out(2) !SO2 aq
     $           + cst_oxydation(2) * conc_out(3) !HSO3-
     $           + cst_oxydation(3) * conc_out(4)) !SO3--
     $           * aqo3

C     Oxidation by H2O2 with the reaction rate from VSRM.
            dsivdt_h2o2 = - 1300000.d0 * exp( -4430.d0 * coef_cte)
     $           * conc_out(2) * aqh2o2
     $           / (1.d0 + 16.d0 * conc_out(1))

C     Total rate.
            dsivdt = dsivdt_o3 + dsivdt_h2o2

C     Compute concentrations.
C     -----------------------

C     Forecast of S(IV), O3 and H2O2 concentrations.
            fc_siv = conc_in(1) + dsivdt * deltat * cst_RTLWC ! in atm
            fc_o3 = gas(igo3) + dsivdt_o3 * deltat * cst_RTLWC ! in atm
            fc_h2o2 = gas(igh2o2) + dsivdt_h2o2 * deltat * cst_RTLWC ! in atm

            if ( (fc_siv.lt.0.d0) .or.
     $           ( (fc_o3.lt.0.d0) .or. (fc_h2o2.lt.0.d0) ) ) then

C     Determine the maximum time step for the most limiting reactant among
C     SO2, O3 and H2O2.
               deltatmax = DMIN1( (conc_in(1) / cst_RTLWC / (-dsivdt)),
     $              DMIN1((gas(igo3) / cst_RTLWC / (-dsivdt_o3)),
     $              (gas(igh2o2) / cst_RTLWC / (-dsivdt_h2o2))))

               conc_in(1) = DMAX1( 0.d0,
     $              (conc_in(1) + dsivdt * deltatmax * cst_RTLWC ))
               gas(igo3) = DMAX1( 0.d0,
     $              (gas(igo3) + dsivdt_o3 * deltatmax * cst_RTLWC ))
               gas(igh2o2) = DMAX1( 0.d0,
     $              (gas(igh2o2) + dsivdt_h2o2 * deltatmax * cst_RTLWC))
               conc_in(2) = conc_in(2) - dsivdt * deltatmax * cst_RTLWC
            else
               conc_in(1) = fc_siv
               gas(igo3) = fc_o3
               gas(igh2o2) = fc_h2o2
               conc_in(2) = conc_in(2) - dsivdt * deltat * cst_RTLWC
            endif

C     Update aqueous concentration and pH.
C     ------------------------------------

C     An instantaneous equilibrium between the gas phase concentration
C     and the aqueous phase concentration is assumed.
            aqh2o2 = gas(igh2o2) * henrys(igh2o2) /
     $           (1.d0 + cst_RTLWC * henrys(igh2o2)) ! in mol.L-1_water

            aqo3 = gas(igo3) * henrys(igo3) /
     $           (1.d0 + cst_RTLWC * henrys(igo3)) ! in mol.L-1_water

C     Update pH.
            call compute_ph(conc_in, henrys, cst_dissoc,
     $           temp, lwc, ph, conc_out)

         endif
      enddo

C     New total sulfate, nitrate and ammonium.
C     ----------------------------------------

C     Unlike nitrate and ammonium, sulfate is not stored in conc_out.
      sulfate = conc_in(2) * mmSO4 * 1.d9 / cst_RT ! in \mu g.m-3
      nitrate = conc_out(6) * mmNO3 * 1.d9 * lwc ! in \mu g.m-3
      ammonium = conc_out(8) * mmNH4 * 1.d9 * lwc ! in \mu g.m-3

C     Dissolved S(IV) is added to gas phase SO2.
      gas(igso2) = conc_in(1) ! in atm
C     Portion of S(IV) concerned by in-cloud scavenging.
      aqso2 = conc_out(2) + conc_out(3) + conc_out(4) ! in mol.L-1

C     Dissolved HNO3 is added to gas phase HNO3.
      gas(ighno3) = conc_out(9) + conc_out(5) * cst_RTLWC ! in atm

C     Dissolved NH3 is added to gas phase NH3.
      gas(ignh3) = conc_out(10) + conc_out(7) * cst_RTLWC ! in atm

C     Compute new distribution by projection to
C     the aerosol population (as in VSRM).
C     -----------------------------------------

      do isect = 1, nbin_aer
         aerosol(isect, nan) = aerosol(isect, nan)
     $        + fdist(isect) * (nitrate - nitrate_init)
         aerosol(isect, naa) = aerosol(isect, naa)
     $        + fdist(isect) * (ammonium - ammonium_init)
         aerosol(isect, na4) = aerosol(isect, na4)
     $        + fdist(isect) * (sulfate - sulfate_init)
      enddo

      do isp = 1, naers
         do isect = 1, nbin_aer
            aerosol(isect, isp) = dmax1(aerosol(isect, isp), tinyaq2)
         enddo
      enddo

      do i = 1, nbin_aer
         totmass(i) = aerosol(i, naa) + aerosol(i, na4)
     $        + aerosol(i, nan) + aerosol(i, nac) + aerosol(i, nas)
     $        + aerosol(i, nao) + aerosol(i, nae) + aerosol(i, nar)
         dnew(i) = (totmass(i) / numberconc(i) / density_aer(i)
     $        / cst_pi6)**cst_FRAC3
      enddo

      call redistaq_VSRM(nbin_aer, dbf_aero, xbf_aero,
     $     fixed_density_aer, dnew, aerosol)

C     14) Update gas & aerosol concentrations in model.
C     -------------------------------------------------

      ZA(ictmHNO3) = gas(ighno3) * 1.d9 * mmHNO3 / cst_RT !HNO3
      ZA(ictmSO2)  = gas(igso2)  * 1.d9 * mmSO2  / cst_RT !SO2
      ZA(ictmH2O2) = gas(igh2o2) * 1.d9 * mmH2O2 / cst_RT !H2O2
      ZA(ictmO3)   = gas(igo3)   * 1.d9 * mmO3   / cst_RT !O3
      do i = 1, nbin_aer
         ii = ns + i
         ZA(ii + (ENa -1) * nbin_aer) = aerosol(i, nas) ! Na
         ZA(ii + (ESO4-1) * nbin_aer) = aerosol(i, na4) ! SO4
         ZA(ii + (ENH3-1) * nbin_aer) = aerosol(i, naa) ! NH4
         ZA(ii + (ENO3-1) * nbin_aer) = aerosol(i, nan) ! NO3
         ZA(ii + (ECl -1) * nbin_aer) = aerosol(i, nac) ! Cl
         ZA(ii + (EMD -1) * nbin_aer) = aerosol(i, nar) ! DUST
         ZA(ii + (EBC -1) * nbin_aer) = aerosol(i, nae) ! EC
         ZA(ii + (EPOA-1) * nbin_aer) = aerosol(i, nao) ! POA
         ZA(ii + (EH2O-1) * nbin_aer) = aerosol(i, naw) ! H2O
      enddo

C     Computes in-cloud scavenging.
C     -----------------------------

      if (rain_rate.gt.0.d0) then
C     Compute the in-cloud scavenging coefficient (same for all species).
         drop_diam = 9.76d-4 * (rain_rate**0.21d0 ) ! m
         scav_coef = 4.17d-7 * collision_eff * rain_rate / drop_diam ! s-1

C     Compute scavenged fraction of droplet.
         scav_frac = 1.d0 - dexp(scav_coef * (t0 - t1) ) ! adim

C     Compute scavenged quantity.
         do i = ifirstact, nbin_aer
            ii = ns + i
            do isp = EMD, EPOA
               qscav(ii + (isp-1) * nbin_aer) =
     $              ZA(ii + (isp-1) * nbin_aer) * scav_frac
            enddo
         enddo

C     Conversion from mol.L-1 to \mu g.m-3.
         qscav(ictmSO2) = aqSO2 * lwc * 1.d9 * mmSO2
     $        * scav_frac * 0.7901d0 ! Remarks: what is the meaning of the coefficient 0.7901d0?
         qscav(ictmH2O2)= aqH2O2 * lwc * 1.d9 * mmH2O2 * scav_frac

C     Remove scavenged quantity.
         do i = ifirstact, nbin_aer
            ii = ns + i
            do isp = EMD, EPOA
               ZA(ii + (isp-1) * nbin_aer) = ZA(ii + (isp-1) * nbin_aer)
     $              - qscav(ii + (isp-1) * nbin_aer)
            enddo
         enddo
         ZA(ictmSO2) = ZA(ictmSO2) - qscav(ictmSO2)
         ZA(ictmH2O2) = ZA(ictmH2O2) - qscav(ictmH2O2)
      endif

C     Treshold TINYAQ.
      do i = 1, Nesp
         ZA(i) = dmax1(ZA(i), TINYAQ)
      enddo

      return
      end
