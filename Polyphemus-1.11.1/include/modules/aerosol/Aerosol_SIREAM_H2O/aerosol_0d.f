C-----------------------------------------------------------------------
C     Copyright (C) 2003-2012, ENPC - INRIA - EDF R&D
C     Author(s): Debry Edouard, Youngseob Kim
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

      subroutine aerosol_0d (nesp, nesp_isorropia_loc, nesp_aec_loc,
     $     nesp_pankow_loc, nesp_poa_loc, LWCmin, ts, DLhumid,
     $     DLtemp, DLpress, delta_t, DLconc, noptions_aer,
     $     options_aer, nesp_aer, nbin_aer, ncycle_aer, DLLWC, DLrain,
     $     bin_bound_aer, fixed_density_aer, density_aer,
     $     couples_coag, first_index_coag, second_index_coag,
     $     coefficient_coag, vertical_interface,
     $     DLconc_aer, Wet_Deposition, Wet_Deposition_aer,
     $     pH, saturation_pressure, partition_coefficient,
     $     vaporization_enthalpy, accomodation_coefficient,
     $     surface_tension, saturation_pressure_mass,
     $     saturation_pressure_torr, deliquescence_relative_humidity,
     $     molecular_weight_aer, molecular_diameter_aer,
     $     collision_factor_aer, mass_density_aer,
     $     aerosol_species_interact_loc, isorropia_species_loc,
     $     aec_species_loc, pankow_species_loc,
     $     poa_species_loc, md_species_loc, bc_species_loc,
     $     nesp_cloud_interact,cloud_species_interact,
     $     lwcavg, heightfog, ifog, DLnum_conc_aer, 
     $     Wet_Deposition_Number_aer, DQLIMIT, p_nucl_fact, k_nucl_fact)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes one timestep for aerosol chemistry SIREAM.
C     Chemical kinetics is solved in each grid cell.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NESP: number of gas species.
C     LWCmin: air liquid water content threshold.
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLHUMID: specific humidity field at initial time ([%]).
C     DLTEMP: temperature field at initial time ([K]).
C     DLPRESS: pressure field at initial time ([Pa]).
C     DELTA_T: time step ([s]).
C     NOPTIONS_AER: number of aerosol module options.
C     OPTIONS_AER: list of aerosol module options.
C     NESP_AER: number of aerosol species.
C     NBIN_AER: number of aerosol bins.
C     NCYCLE_AER: number of cycle in aerosol computation.
C     DLLWC: air liquid water content ([fraction]).
C     DLRAIN: field of rain rate at initial time.
C     BIN_BOUND_AER: aerosol diameters at bin bounds.
C     FIXED_DENSITY_AER: fixed aerosol density ([kg/m^3]).
C     DENSITY_AER: size variable aerosol density ([kg/m^3]).
C     COUPLES_COAG: coagulation couples for each bin.
C     FIRST_INDEX_COAG: first bin index of coagulation couples.
C     SECOND_INDEX_COAG: second bin index of coagulation couples.
C     COEFFICIENT_COAG: coagulation partition coefficient.
C     VERTICAL_INTERFACE: 1D field of vertical interface height ([m]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of gas concentrations ([\mu.g/m^3]).
C     DLCONC_AER: array of aerosol concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
C     DLNUM_CONC_AER: array of aerosol number concentration ([m^(-3)])
C     
C     -- OUTPUT VARIABLES
C
C     Wet_Deposition: 1D wet fluxes of gaseous species due to in-cloud
C     scavenging ([\mu.g/m^2/s]).
C     Wet_Deposition_aer: 2D wet fluxes of particulate species due to
C     in-cloud scavenging ([\mu.g/m^2/s]).
C     Wet_Deposition_Number_aer: 2D wet number fluxes of particulate species due to 
C     in-cloud scavenging ([m^-2/s]).
C     
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2002/02/26: new treatment of sources (Jaouad Boutahar, CEREA).
C     2006/09/29: updated header (Edouard Debry).
C     2010/03/04: SIMPLE_AQUEOUS option included (Youngseob KIM)
C     2013/11/27: Added number concentrations (Stephanie Deschamps, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Edouard Debry, CEREA, September 2006.
C
C------------------------------------------------------------------------

      implicit none

      include 'CONST.INC'
      include 'CONST_A.INC'
      include 'param.inc'
      include 'dynaero.inc'
      include 'varg.inc'
      include 'vara.inc'

      double precision ts, delta_t

      integer nesp

      double precision DLconc(nesp), DLtemp
      double precision DLhumid, DLpress

      integer Jt, Jsp, i, Jb

      integer nesp_aer, nbin_aer, ncycle_aer
      integer noptions_aer, options_aer(noptions_aer)
      double precision DLLWC, DLrain
      double precision bin_bound_aer(nbin_aer + 1)
      double precision density_aer(nbin_aer), fixed_density_aer
      integer couples_coag(nbin_aer)
      integer first_index_coag(nbin_aer, 4 * nbin_aer)
      integer second_index_coag(nbin_aer, 4 * nbin_aer)
      double precision coefficient_coag(nbin_aer, nbin_aer, nbin_aer)
      double precision vertical_interface(2)
      double precision layer_height
      double precision DLconc_aer(nbin_aer, nesp_aer)
      DOUBLE PRECISION DLnum_conc_aer(nbin_aer)

      double precision Wet_Deposition(nesp)
      double precision Wet_Deposition_aer(nbin_aer, nesp_aer)
      DOUBLE PRECISION Wet_Deposition_Number_aer(nbin_aer)

      integer ICLD, IINCLD
      double precision rho_aero(nbin_aer+1)

      double precision lwca, pH, lwc_surf
      double precision tschem_aer, tfchem_aer, dtchem_aer
      double precision ZA(nesp + nbin_aer * nesp_aer)
      DOUBLE PRECISION ZNA(nbin_aer)
      double precision rain_rate
      double precision qscav(nesp + nbin_aer * nesp_aer)
      DOUBLE PRECISION qscav_num(nbin_aer)
      double precision LWCmin

      double precision XSF(nbin_aer), MSF(nbin_aer)
      double precision DSF(nbin_aer), XBF(nbin_aer+1)
      double precision MBF(nbin_aer+1), DBF(nbin_aer+1)
      double precision HSF(nbin_aer)

      integer iq(nesp_aer,nbin_aer)

      double precision saturation_pressure(nesp_aer)
      double precision partition_coefficient(nesp_aer)
      double precision vaporization_enthalpy(nesp_aer)
      double precision accomodation_coefficient(nesp_aer)
      double precision surface_tension(nesp_aer)
      double precision saturation_pressure_mass(nesp_aer)
      double precision saturation_pressure_torr(nesp_aer)
      double precision deliquescence_relative_humidity(nesp_aer)

      double precision molecular_weight_aer(nesp_aer)
      double precision molecular_diameter_aer(nesp_aer)
      double precision collision_factor_aer(nesp_aer)
      double precision mass_density_aer(nesp_aer)
      integer aerosol_species_interact_loc(nesp_aer)

      integer nesp_isorropia_loc, nesp_aec_loc, nesp_pankow_loc
      integer nesp_poa_loc
      integer isorropia_species_loc(nesp_isorropia_loc)
      integer aec_species_loc(nesp_aec_loc)
      integer pankow_species_loc(nesp_pankow_loc)
      integer poa_species_loc(nesp_poa_loc)
      integer md_species_loc, bc_species_loc
      integer nesp_cloud_interact
      integer cloud_species_interact(nesp_cloud_interact)

      double precision lwcavg, heightfog
      integer ifog

      INTEGER section_pass
      double precision conc_tot
      double precision DQLIMIT
      double precision k_nucl_fact, p_nucl_fact
      double precision diam, rho

      integer type
      DOUBLE PRECISION DSF2(nbin_aer), dold(nbin_aer)
      DOUBLE PRECISION conc(nbin_aer, nesp_aer)

      if (nesp_aer.gt.NEXT) then
         stop 'Number of aerosol species exceeds NEXT in param.inc'
      endif

      k_fact = k_nucl_fact
      p_fact = p_nucl_fact

C     Set parameters
      ICOAG = options_aer(1)
      ICOND = options_aer(2)
      INUCL = options_aer(3)
      ICLD = options_aer(4)
      IINCLD = options_aer(5)
      IKELV = options_aer(6)
      IDENS = options_aer(7)
      ICUT = options_aer(8)
      ISULFCOND = options_aer(9)
      KDSLV = options_aer(10)
      IREDIST = options_aer(11)
      ITERN = options_aer(12)
      ITHRM = options_aer(13)
      IOLIGO = options_aer(14)
      ITHERMO = options_aer(15)
      INUM = options_aer(16)

      do Jsp = 1, NEXT
         PSATREF(Jsp) = 0.d0
         KPARTREF(Jsp) = 0.d0
         DHVAP(Jsp) = 0.d0
         SIGMA(Jsp) = 0.d0
         STICK(Jsp) = 0.d0
         QSATREF(Jsp) = 0.d0
         TSATREF(Jsp) = 0.d0
         DRH(Jsp) = 0.d0
         EMW(Jsp) = 0.d0
         SIGM(Jsp) = 0.d0
         PARM(Jsp) = 0.d0
         LMD(Jsp) = 0.d0
      enddo

      do Jsp = 1, nesp_aer
         PSATREF(Jsp) = saturation_pressure(Jsp)
         KPARTREF(Jsp) = partition_coefficient(Jsp)
         DHVAP(Jsp) = vaporization_enthalpy(Jsp)
         SIGMA(Jsp) = surface_tension(Jsp)
         STICK(Jsp) = accomodation_coefficient(Jsp)
         QSATREF(Jsp) = saturation_pressure_mass(Jsp)
         TSATREF(Jsp) = saturation_pressure_torr(Jsp)
         DRH(Jsp) = deliquescence_relative_humidity(Jsp)
         EMW(Jsp) = molecular_weight_aer(Jsp)
         SIGM(Jsp) = molecular_diameter_aer(Jsp)
         PARM(Jsp) = collision_factor_aer(Jsp)
         LMD(Jsp) = mass_density_aer(Jsp)
      enddo

C     Set pointers.
      call initpoint(nbin_aer, nesp_aer, aerosol_species_interact_loc,
     $     iq, nesp_isorropia_loc, isorropia_species_loc, nesp_aec_loc,
     $     aec_species_loc, nesp_pankow_loc, pankow_species_loc,
     $     nesp_poa_loc, poa_species_loc, md_species_loc,
     $     bc_species_loc, nesp_cloud_interact, cloud_species_interact)

C     Aerosol density converted from kg / m^3 to microg / microm^3.
      RHOA = fixed_density_aer * 1.d-09
      do Jb = 1, nbin_aer
         rho_aero(Jb) = density_aer(Jb) * 1.d-9
      enddo

C     Aerosol discretization converted in microm.
      do Jb = 1, nbin_aer + 1
         DBF(Jb) = bin_bound_aer(Jb) * 1.d6
         MBF(Jb) = RHOA * cst_pi6 * DBF(Jb)**3
         XBF(Jb) = dlog(MBF(Jb))
      enddo

      do Jb = 1, nbin_aer
         XSF(Jb) = (XBF(Jb) + XBF(Jb+1)) * 5.d-1
         MSF(Jb) = dsqrt(MBF(Jb) * MBF(Jb + 1))
         DSF(Jb) = dsqrt(DBF(Jb) * DBF(Jb + 1))
      enddo

C     Width of each fixed bin.
      do Jb = 1, nbin_aer
         HSF(Jb) = XBF(Jb + 1) - XBF(Jb)
      enddo

!     Section_pass for the repartition_euler module
      section_pass = 1;
      DO WHILE(diam_pass.GT.DBF(section_pass+1) 
     &     .AND. diam_pass.LE.DBF(nbin_aer+1))
         section_pass = section_pass + 1
      ENDDO
      IF (section_pass .EQ. 1) THEN
         PRINT * , "These scheme is not adapt, 
     & please change redistribution method"
      ENDIF

C     Cloud liquid water content (g/m^3).
      lwca = DLLWC * 1.d3 * DLpress / 101325.d0 * 28.97d0 / Pr / DLtemp

C     Initialize layer_height.
      layer_height = vertical_interface(2) - vertical_interface(1)

C     Initialize pH.
      pH = 0.d0

C     Initialize in-cloud wet fluxes.
      do Jsp = 1, nesp
         Wet_Deposition(Jsp) = 0.d0
      enddo

      do Jsp = 1, nesp_aer
         do Jb = 1, nbin_aer
            Wet_Deposition_aer(Jb, Jsp) = 0.d0
         enddo
      enddo

      IF (INUM.EQ.1) THEN
         DO Jb=1,Nbin_aer
            Wet_Deposition_Number_aer(Jb) = 0.d0
         ENDDO
      ENDIF

      DO Jsp=1,NESP
         ZA(Jsp) = DLconc(Jsp)
         IF (ZA(Jsp).LT.0.D0) then
            write(*,*) "negative concentrations", ZA(Jsp)
            ZA(Jsp) = 0.D0
         ENDIF
      ENDDO
      
      DO Jb=1,nbin_aer
         conc_tot = 0.d0
         DO Jsp=1,nesp_aer
            i = nesp+(Jsp-1)*nbin_aer+Jb
            if (DLconc_aer(Jb,Jsp) .GT. TINYM) then
               ZA(i) = DLconc_aer(Jb,Jsp)
               conc_tot = conc_tot + ZA(i)
            else
               ZA(i) = 0.d0
            endif
         ENDDO
         IF ((conc_tot .LE. TINYM).AND.(INUM.EQ.1)) then
            ZNA(Jb) = 0.d0
         ENDIF
      ENDDO

      IF (INUM.EQ.1) THEN
         DO Jb=1, nbin_aer
            if (DLnum_conc_aer(Jb) .GT. TINYN) then
               ZNA(Jb) = DLnum_conc_aer(Jb)   
            else
               ZNA(Jb) = 0.d0
               DO Jsp=1,nesp_aer
                  i = nesp + (Jsp - 1)*nbin_aer + Jb
                  ZA(i) = 0.d0
               ENDDO
            endif
         ENDDO
      ENDIF

      IF (IDENS .EQ. 1) THEN
         DO Jb=1,nbin_aer 
            CALL compute_density(nbin_aer,nesp_aer, nesp_aer,TINYM,
     &                            DLconc_aer,LMD,Jb,rho_aero(Jb))
         ENDDO
      ENDIF
      
      DO Jb=1,nbin_aer
         conc_tot = 0.d0
         DO Jsp = 1, Nesp_aer-1
            i = nesp+(Jsp-1)*nbin_aer+Jb
            conc_tot = conc_tot + ZA(i)
         ENDDO
         IF (INUM.EQ.1 .AND. ZNA(Jb).GT.TINYN
     s        .AND. IDENS .EQ. 1) then
            DSF(Jb) = (conc_tot/ZNA(Jb)
     &           /cst_PI6/rho_aero(Jb))**cst_FRAC3
         ELSE
            DSF(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
         ENDIF
      ENDDO

      IF(ncycle_aer.NE.1) then
         write(*,*) 'Pb ncycle_aer > 1', ncycle_aer
         stop
      ENDIF

C     Loop on grid cells.
      do Jt = 1, ncycle_aer
         tschem_aer = ts + (Jt - 1) * delta_t / ncycle_aer
         tfchem_aer = tschem_aer + delta_t / ncycle_aer
         dtchem_aer = tfchem_aer - tschem_aer


         DO Jb=1,nbin_aer
            conc_tot = 0.d0
            DO Jsp = 1, Nesp_aer-1
               i = nesp+(Jsp-1)*nbin_aer+Jb
               conc_tot = conc_tot + ZA(i)
               conc(Jb, Jsp) = ZA(i)

            ENDDO
            IF (IDENS .EQ. 1) THEN
               CALL compute_density(nbin_aer,nesp_aer, nesp_aer,TINYM,
     &              conc,
     &              LMD,Jb,rho_aero(Jb))
            ENDIF
            IF (DLnum_conc_aer(Jb).GT.0.D0) THEN
               dold(Jb) = (conc_tot/DLnum_conc_aer(Jb)
     &              /cst_PI6/rho_aero(Jb))**cst_FRAC3
            ELSE
               dold(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
            ENDIF
         ENDDO

         if (lwca.ge.LWCmin) then

            if (IINCLD.eq.1) then
               rain_rate = DLrain
            else
               rain_rate = 0.d0
            endif

            if (ICLD.eq.1) then ! Use VSRM chemical mechanism.

               DO Jb=1,nbin_aer
                  conc_tot = 0.d0
                  DO Jsp = 1, Nesp_aer-1
                     i = nesp+(Jsp - 1)*nbin_aer+Jb
                     conc_tot = conc_tot + ZA(i)
                     conc(Jb, Jsp) = ZA(i)
                  ENDDO
                  IF (IDENS .EQ. 1) THEN
                     CALL compute_density(nbin_aer,nesp_aer,nesp_aer,
     &                    TINYM,conc,
     &                    LMD,Jb,rho_aero(Jb))
                  ENDIF
                  IF (DLnum_conc_aer(Jb).GT.0.D0) THEN
                     DSF2(Jb) = (conc_tot/ZNA(Jb)
     &                    /cst_PI6/rho_aero(Jb))**cst_FRAC3
                  ELSE
                     DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
                  ENDIF
               ENDDO

               call vsrmchem(nesp, nesp_aer, nbin_aer, rho_aero, RHOA,
     $              DBF, DSF, XBF, ZA, DLhumid, DLpress, DLtemp, lwca,
     $              tschem_aer, tfchem_aer, rain_rate, pH, qscav,
     $              ZNA, qscav_num,
     $              section_pass,DQLIMIT)

               DO Jb=1,nbin_aer
                  conc_tot = 0.d0
                  DO Jsp = 1, Nesp_aer-1
                     i = nesp+(Jsp - 1)*nbin_aer+Jb
                     conc_tot = conc_tot + ZA(i)
                     conc(Jb, Jsp) = ZA(i)
                  ENDDO
                  IF (IDENS .EQ. 1) THEN
                     CALL compute_density(nbin_aer,nesp_aer,nesp_aer,
     &                 TINYM,conc,
     &                 LMD,Jb,rho_aero(Jb))
                  ENDIF

                  IF (ZNA(Jb).GT.TINYN) THEN
                     DSF2(Jb) = (conc_tot/ZNA(Jb)
     &                    /cst_PI6/rho_aero(Jb))**cst_FRAC3
                  ELSE
                     DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
                  ENDIF
               ENDDO

               do Jsp = 1, nesp
                  Wet_Deposition(Jsp) = Wet_Deposition(Jsp) +
     $                 qscav(Jsp) * layer_height / dtchem_aer
               enddo

               do Jsp = 1, nesp_aer
                  do Jb = 1, nbin_aer
                     Wet_Deposition_aer(Jb, Jsp) =
     $                    Wet_Deposition_aer(Jb, Jsp) +
     $                    qscav(nesp + nbin_aer*(Jsp - 1) + Jb) *
     $                    layer_height / dtchem_aer
                  enddo
               enddo

               IF (INUM.EQ.1) THEN
                  
                  DO Jb=1,Nbin_aer
                     Wet_Deposition_Number_aer(Jb) =  
     $                    Wet_Deposition_Number_aer(Jb) +
     $                    qscav_num(Jb) * 
     $                    layer_height / dtchem_aer
                  ENDDO
                  
               ENDIF

            else if(ICLD.EQ.2) then ! Use "simple aqueous" chemical mechanism.
               call simple_aqueous_module(nesp, nesp_aer, nbin_aer,
     $              rho_aero, RHOA, DBF, DSF, XBF, ZA, DLpress,
     $              DLtemp, lwca, tschem_aer, tfchem_aer, rain_rate, pH,
     $              qscav, section_pass, ZNA, DQLIMIT)

               do Jsp = 1, nesp
                  Wet_Deposition(Jsp) = Wet_Deposition(Jsp)
     $                 + qscav(Jsp) * layer_height / dtchem_aer
               enddo

               do Jsp = 1, nesp_aer
                  do Jb = 1, nbin_aer
                     Wet_Deposition_aer(Jb, Jsp) =
     $                    Wet_Deposition_aer(Jb,Jsp)
     $                    + qscav(nesp + nbin_aer * (Jsp - 1) + Jb)
     $                    * layer_height / dtchem_aer
                  enddo
               enddo

               IF (INUM.EQ.1) THEN
                  DO Jb=1,Nbin_aer
                     Wet_Deposition_Number_aer(Jb) =  
     $                    Wet_Deposition_Number_aer(Jb) +
     $                    qscav_num(Jb) * 
     $                    layer_height / dtchem_aer
                  ENDDO
               ENDIF


            else
               call aerodyn(nesp, nesp_aer, nbin_aer, DLtemp, DLpress,
     $              DLhumid, tschem_aer, tfchem_aer, ZA, couples_coag,
     $              first_index_coag, second_index_coag,
     $              coefficient_coag, XSF, MSF, DSF, XBF, MBF, DBF,
     $              HSF,iq, ZNA, section_pass,DQLIMIT)

               DO Jb=1,nbin_aer
                  conc_tot = 0.d0
                  DO Jsp = 1, Nesp_aer-1
                     i = nesp+(Jsp - 1)*nbin_aer+Jb
                     conc_tot = conc_tot + ZA(i)
                     conc(Jb, Jsp) = ZA(i)
                  ENDDO
                  IF (IDENS .EQ. 1) THEN
                     CALL compute_density(nbin_aer,nesp_aer,nesp_aer,
     &                 TINYM,conc,
     &                 LMD,Jb,rho_aero(Jb))
                  ENDIF
                  IF (ZNA(Jb).GT.TINYN) THEN
                     DSF2(Jb) = (conc_tot/ZNA(Jb)
     &                    /cst_PI6/rho_aero(Jb))**cst_FRAC3
                  ELSE
                     DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
                  ENDIF
               ENDDO

            endif

         else
            call aerodyn(nesp, nesp_aer, nbin_aer, DLtemp, DLpress,
     $           DLhumid, tschem_aer, tfchem_aer, ZA, couples_coag,
     $           first_index_coag, second_index_coag, coefficient_coag,
     $           XSF, MSF, DSF, XBF, MBF, DBF, HSF, iq, ZNA, 
     $           section_pass, DQLIMIT)

            DO Jb=1,nbin_aer
               conc_tot = 0.d0
               DO Jsp = 1, Nesp_aer-1
                  i = nesp+(Jsp - 1)*nbin_aer+Jb
                  conc_tot = conc_tot + ZA(i)
                  conc(Jb, Jsp) = ZA(i)
               ENDDO
               IF (IDENS .EQ. 1) THEN
                  CALL compute_density(nbin_aer,nesp_aer,nesp_aer,TINYM,
     &                 conc,
     &                 LMD,Jb,rho_aero(Jb))
               ENDIF
               IF (ZNA(Jb).GT.TINYN) THEN
                  DSF2(Jb) = (conc_tot/ZNA(Jb)
     &                 /cst_PI6/rho_aero(Jb))**cst_FRAC3
               ELSE
                  DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
               ENDIF
            ENDDO

         ENDIF
         DO Jb=1,nbin_aer
            conc_tot = 0.d0
            DO Jsp = 1, Nesp_aer-1
               i = nesp+(Jsp - 1)*nbin_aer+Jb
               conc_tot = conc_tot + ZA(i)
               conc(Jb, Jsp) = ZA(i)
            ENDDO
            IF (IDENS .EQ. 1) THEN
            CALL compute_density(nbin_aer,nesp_aer,nesp_aer,TINYM,
     &              conc,
     &              LMD,Jb,rho_aero(Jb))
            ENDIF
            IF (ZNA(Jb).GT.TINYN) THEN
               DSF2(Jb) = (conc_tot/ZNA(Jb)
     &              /cst_PI6/rho_aero(Jb))**cst_FRAC3
            ELSE
               DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
            ENDIF
         ENDDO
      ENDDO                     !End of loop on Ncycle_aer
      DO Jb=1,nbin_aer
         conc_tot = 0.d0
         DO Jsp = 1, Nesp_aer-1
            i = nesp+(Jsp - 1)*nbin_aer+Jb
            conc_tot = conc_tot + ZA(i)
            conc(Jb, Jsp) = ZA(i)
         ENDDO
         IF (IDENS .EQ. 1) THEN
         CALL compute_density(nbin_aer,nesp_aer,nesp_aer,TINYM,
     &        conc,
     &        LMD,Jb,rho_aero(Jb))
         ENDIF

         IF (ZNA(Jb).GT.TINYN) THEN
            DSF2(Jb) = (conc_tot/ZNA(Jb)
     &           /cst_PI6/rho_aero(Jb))**cst_FRAC3
         ELSE
            DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
         ENDIF
      ENDDO
      DO Jb=1,nbin_aer
         DO Jsp=1,nesp_aer
            i = nesp+(Jsp - 1)*nbin_aer+Jb
            DLconc_aer(Jb,Jsp) = ZA(i)
            IF (isNaN(DLconc_aer(Jb,Jsp))) then
               write(*,*) jb,jsp,dlconc_aer(jb,jsp)
               STOP "aerosol.f :  concentration is NaN"
            endif
         ENDDO
      ENDDO

C     Number concentration 
      IF(INUM.EQ.1) THEN
         DO Jb=1,nbin_aer
            DLnum_conc_aer(Jb) = ZNA(Jb)
            if (isNaN(DLnum_conc_aer(Jb))) then 
               write(*,*), "aerosol.f for number concentration :",
     &               Jb, DLnum_conc_aer(Jb)
                stop
             endif
          ENDDO
       ENDIF
       
       DO Jsp=1,nesp
          DLconc(Jsp) = ZA(Jsp)
          IF (isNaN(dlconc(Jsp))) THEN 
             WRITE(*,*) Jsp, dlconc(Jsp)
             STOP "aerosol.f :  concentration is NaN"
          ENDIF
       ENDDO

       do Jb = 1, nbin_aer
          conc_tot = 0.d0
          DO Jsp = 1, Nesp_aer-1
             conc_tot = conc_tot + DLconc_aer(Jb, Jsp)
          ENDDO
          IF (IDENS .EQ. 1) THEN
             CALL compute_density(nbin_aer,nesp_aer, nesp_aer,TINYM,
     &            DLconc_aer,
     &            LMD,Jb,rho_aero(Jb))
          ENDIF
          IF (DLnum_conc_aer(Jb).GT.0.D0) THEN
             DSF2(Jb) = (conc_tot/DLnum_conc_aer(Jb)
     &            /cst_PI6/rho_aero(Jb))**cst_FRAC3
          ELSE
             DSF2(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
          ENDIF
         
       enddo

         IF (ICLD.GE.1.AND.IINCLD.EQ.1.and.ifog.eq.1) THEN
            call settling_0d(nesp_aer,nbin_aer,
     &           DSF,DLconc_aer,
     &           lwcavg,tschem_aer,tfchem_aer,heightfog)
         ENDIF


      END
