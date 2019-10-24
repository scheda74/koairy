Module zAerosolSCRAM
  use aInitialization
  use jAdaptstep
  use omp_lib
  implicit none
!!-----------------------------------------------------------------------
!!     Copyright (C) 2003-2018, CERA (ENPC - EDF R&D)
!!     Author(s): Debry Edouard, Shupeng Zhu, Youngseob Kim
!!
!!     This file is part of the SCRAM, a
!!     component of the air quality modeling system Polyphemus.
!!
!!     Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
!!
!!     Polyphemus is free software; you can redistribute it and/or modify
!!     it under the terms of the GNU General Public License as published
!!     by the Free Software Foundation; either version 2 of the License,
!!     or (at your option) any later version.
!!
!!     Polyphemus is distributed in the hope that it will be useful, but
!!     WITHOUT ANY WARRANTY; without even the implied warranty of
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!!     General Public License for more details.
!!
!!     For more information, visit the Polyphemus web site:
!!     http://cerea.enpc.fr/polyphemus/
!!-----------------------------------------------------------------------
  contains
      SUBROUTINE aerosol(nesp_loc,&
            nesp_isorropia_loc,nesp_aec_loc,nesp_pankow_loc,&
            nesp_poa_loc,&
            LWCmin,initial_time,DLhumid,&
            DLtemp,DLpress,delta_t,DLconc,noptions_aer,&
            options_aer,nesp_aer,&
       	    nsize_section_aer,nbin_aer,ncomposition_aer,ngroup_aer,&
            ncycle_aer,DLLWC,DLrain,&
            bin_bound_aer,fixed_density_aer,&
            aerosol_species_group_relation,composition_bounds,&
            Ncoefficient,index_first,index_second,coefficient,&
            coef_size, vertical_interface,&
            DLconc_aer,Wet_Deposition,Wet_Deposition_aer,&
            pH,saturation_pressure_loc, partition_coefficient_loc,&
            vaporization_enthalpy_loc, accomodation_coefficient_loc,&
            surface_tension_loc, saturation_pressure_mass_loc,&
            saturation_pressure_torr_loc, deliquescence_relative_humidity_loc,&
            molecular_weight_aer_loc, molecular_diameter_aer_loc,&
            collision_factor_aer_loc, mass_density_aer_loc,&
            aerosol_species_interact_loc, isorropia_species_loc,&
            aec_species_loc, pankow_species_loc,&
            poa_species_loc, md_species_loc, bc_species_loc,&
            nesp_cloud_interact,cloud_species_interact,&
            lwcavg, heightfog, ifog, DLnum_conc_aer, &
            Wet_Deposition_Number_aer, DQLIMIT_loc, p_nucl_fact, k_nucl_fact, &
            psoap_config, psurrogate, dtaeromin_loc, relative_error)

!!------------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!
!!     This routine computes one timestep for gas-phase chemistry RADM.
!!     Chemical kinetics is solved in each grid cell.
!!
!!------------------------------------------------------------------------
!!
!!     -- INPUT VARIABLES
!!
!!     nesp_loc: number of gas species.
!!     nesp_isorropia_loc: number of isorropia species
!!     nesp_aec_loc: number of aec species
!!     nesp_pankow_loc: number of pankow species
!!     nesp_poa_loc: number of poa species
!!     LWCmin: air liquid water content threshold.
!!      initial time: initial time (GMT, computed from January 1st, [s]).
!!     DLHUMID: 3D specific humidity field at initial time ([%]).
!!     DLTEMP: 3D temperature field at initial time ([K]).
!!     DLPRESS: 3D pressure field at initial time ([Pa]).
!!     DELTA_T: time step ([s]).
!!     NOPTIONS_AER: number of aerosol module options.
!!     OPTIONS_AER: 1D list of aerosol module options.
!!     NESP_AER: number of aerosol species.
!!     nsize_section_aer: number of aerosol size sections.
!!     ncomposition_aer: number of aerosol composition sections.
!!     NBIN_AER: number of aerosol bins=nsize_section_aer*ncomposition_aer.
!!     NCYCLE_AER: number of cycle in aerosol computation.
!!     DLLWC: 3D air liquid water content ([fraction]).
!!     DLRAIN: 2D field of rain rate at initial time.
!!     BIN_BOUND_AER: aerosol diameters at bin bounds.
!!     FIXED_DENSITY_AER: fixed aerosol density ([g/m^3]).
!!     DENSITY_AER: size variable aerosol density ([g/m^3]).
!!     aerosol_species_group_relation: relations between aerosol species index and groups index
!!     composition_bounds: List of aerosol compositions bounds
!!     coagulation_coefficient_file: location of coagulation repartition database
!!
!!     -- INPUT/OUTPUT VARIABLES
!!
!!     DLCONC: array of 3D gas concentrations ([\mu.g/m^3]).
!!     DLCONC_AER: array of 3D aerosol concentrations ([\mu.g/m^3]).
!!     # Before entry, it is given at initial time of the timestep.
!!     # On exit, it is computed at final time of the timestep.
!!     DLNUM_CONC_AER: array of aerosol number concentration ([m^(-3)])
!!
!!     -- OUTPUT VARIABLES
!!
!!     Wet_Deposition: 2D wet fluxes of gaseous species due to in-cloud
!!     scavenging ([\mu.g/m^2/s]).
!!     Wet_Deposition_aer: 2D wet fluxes of particulate species due to
!!     in-cloud scavenging ([\mu.g/m^2/s]).
!!     Wet_Deposition_Number_aer: 2D wet number fluxes of particulate species due to
!!     in-cloud scavenging ([m^-2/s]).
!!
!!------------------------------------------------------------------------
!!
!!
!!     -- AUTHOR(S)
!!
!!     Shupeng ZHU, CEREA, August 2014.
!!     Youngseob Kim, CEREA, May 2018.
!!
!!------------------------------------------------------------------------

      DOUBLE PRECISION initial_time,delta_t, dtaeromin_loc

      INTEGER nesp_loc
      INTEGER size_cfg

      DOUBLE PRECISION DLtemp
      DOUBLE PRECISION DLhumid,DLpress
      DOUBLE PRECISION LWCmin
      INTEGER jesp,j,f,k,t,Jt,b,g,i
      DOUBLE PRECISION tschem_aer, tfchem_aer, dtchem_aer!time for chemistry
      
      INTEGER nesp_aer,nbin_aer,ncycle_aer
      INTEGER nsize_section_aer,ncomposition_aer,ngroup_aer
      INTEGER noptions_aer,options_aer(noptions_aer)
      INTEGER aerosol_species_group_relation(nesp_aer)
      INTEGER coef_size
      DOUBLE PRECISION DLLWC,DLrain
      DOUBLE PRECISION pH,lwc_surf
      DOUBLE PRECISION rain_rate
      DOUBLE PRECISION concentration_gas_loc(nesp_loc)      
      DOUBLE PRECISION qscav_gas(nesp_loc)
      DOUBLE PRECISION qscav_aer(nbin_aer,nesp_aer)
      DOUBLE PRECISION qscav_num(nbin_aer)      
      DOUBLE PRECISION bin_bound_aer(nsize_section_aer + 1)
      DOUBLE PRECISION composition_bounds(ncomposition_aer*ngroup_aer*2)
      DOUBLE PRECISION fixed_density_aer
      DOUBLE PRECISION vertical_interface(2)
      INTEGER, intent(in) :: Ncoefficient(nbin_aer)!SZ
      INTEGER, intent(in) :: index_first(coef_size)!SZ
      INTEGER, intent(in) :: index_second(coef_size)!SZ
      DOUBLE PRECISION, intent(in) :: coefficient(coef_size)!SZ
      DOUBLE PRECISION, intent(inout) :: DLconc(nesp_loc)
      DOUBLE PRECISION, intent(inout) :: DLconc_aer(nbin_aer,nesp_aer)
      DOUBLE PRECISION, intent(inout) :: DLnum_conc_aer(nbin_aer)

      DOUBLE PRECISION saturation_pressure_loc(nesp_aer)
      DOUBLE PRECISION partition_coefficient_loc(nesp_aer)
      DOUBLE PRECISION vaporization_enthalpy_loc(nesp_aer)
      DOUBLE PRECISION accomodation_coefficient_loc(nesp_aer)
      DOUBLE PRECISION surface_tension_loc(nesp_aer)
      DOUBLE PRECISION saturation_pressure_mass_loc(nesp_aer)
      DOUBLE PRECISION saturation_pressure_torr_loc(nesp_aer)
      DOUBLE PRECISION deliquescence_relative_humidity_loc(nesp_aer)
      DOUBLE PRECISION molecular_weight_aer_loc(nesp_aer)
      DOUBLE PRECISION molecular_diameter_aer_loc(nesp_aer)
      DOUBLE PRECISION collision_factor_aer_loc(nesp_aer)
      DOUBLE PRECISION mass_density_aer_loc(nesp_aer)
    
      INTEGER aerosol_species_interact_loc(nesp_aer)
      
      INTEGER nesp_isorropia_loc,nesp_aec_loc,nesp_pankow_loc
      INTEGER nesp_poa_loc
      INTEGER isorropia_species_loc(nesp_isorropia_loc)
      INTEGER aec_species_loc(nesp_aec_loc)
      INTEGER pankow_species_loc(nesp_pankow_loc)
      INTEGER poa_species_loc(nesp_poa_loc)
      INTEGER md_species_loc,bc_species_loc
      INTEGER nesp_cloud_interact
      INTEGER cloud_species_interact(nesp_cloud_interact)
      !!out put value
      DOUBLE PRECISION, intent(out) :: Wet_Deposition(nesp_loc)
      DOUBLE PRECISION, intent(out) :: Wet_Deposition_aer(nbin_aer,nesp_aer)
      DOUBLE PRECISION, intent(out) :: Wet_Deposition_Number_aer(nbin_aer)

      DOUBLE PRECISION lwcavg, heightfog,cloud_water
      double precision DQLIMIT_loc,layer_height
      INTEGER ifog,tag_file,id
      double precision k_nucl_fact, p_nucl_fact
      double precision total_nb,total_ms,tot_tot
      double precision total_ms_old(nesp_aer)
      double precision relative_error
      
      ! pointers for SOAP 
      integer psoap_config, psurrogate

     
      !!read the coefficient repartition for coagulation

      k_fact = k_nucl_fact
      p_fact = p_nucl_fact
      
      Humidity=DLhumid
      Temperature=DLtemp
      Pressure=DLpress
      dtaeromin = dtaeromin_loc
      epser = relative_error
      
      call COMPUTE_RELATIVE_HUMIDITY(Humidity,Temperature,Pressure,Relative_Humidity)
      Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
      DQLIMIT=DQLIMIT_loc
!!    aInitialization     
    CALL Init_global_parameters(nesp_loc,nbin_aer,nsize_section_aer,nesp_aer,&
				     ncomposition_aer,ngroup_aer,&
				     aerosol_species_interact_loc,&
				     nesp_isorropia_loc,isorropia_species_loc,&
				     nesp_aec_loc,aec_species_loc,&
				     nesp_pankow_loc,pankow_species_loc,&
				     nesp_poa_loc,poa_species_loc,&
				     md_species_loc, bc_species_loc,&
				     nesp_cloud_interact,cloud_species_interact,&
				     saturation_pressure_loc,partition_coefficient_loc,&
				     vaporization_enthalpy_loc,accomodation_coefficient_loc,&
				     surface_tension_loc,saturation_pressure_mass_loc,&
				     saturation_pressure_torr_loc,&
				     deliquescence_relative_humidity_loc,molecular_weight_aer_loc,&
				     molecular_diameter_aer_loc,collision_factor_aer_loc,&
				     mass_density_aer_loc,noptions_aer,options_aer,&
				     aerosol_species_group_relation)


      !!read coagulation_coefficient data
      IF (with_coag.EQ.1) THEN

	call ReadCoefficient(coef_size,Ncoefficient,index_first,index_second,coefficient)
	! defined in ModuleCoefficientRepartition

	  ! Check the quality of coagulation repartition coefficients
	call check_repart_coeff()

      ENDIF

      !!read composition_bounds
      allocate(discretization_composition(N_fracmax, N_groups, 2))
      id=1
      do f=1, N_fracmax
	do g=1,N_groups
	  do k=1,2
	    discretization_composition(f,g,k)=composition_bounds(id)
	    id=id+1
	  enddo
	enddo
      enddo

!!     Cloud liquid water content (g/m^3)
      cloud_water=DLLWC*1000.d0*Pressure/101325.d0*28.97d0/Pr/Temperature

!!     Initialize layer_height
      layer_height = vertical_interface(2) -vertical_interface(1)

!!     Initialize pH
      pH = 0.d0

!!     Initialize in-cloud wet fluxes. (out put)
      do jesp=1, nesp_loc
	Wet_Deposition(jesp)=0.d0
      enddo

      do j=1,nbin_aer
	do jesp=1, nesp_aer
	  Wet_Deposition_aer(j,jesp)=0.d0
	enddo
      enddo
      IF (with_number.EQ.1) THEN
	do j=1,nbin_aer
	  Wet_Deposition_Number_aer(j) = 0.d0
	enddo
      ENDIF

!!     Loop on grid cells.
      DO Jt=1,ncycle_aer
      
      !print*,"start time loop",Jt,ncycle_aer
      
         tschem_aer=initial_time+(Jt-1)*delta_t/ncycle_aer
         tfchem_aer=tschem_aer+delta_t/ncycle_aer
         dtchem_aer=tfchem_aer-tschem_aer
	!initialize gas
         DO jesp=1,nesp_loc ! N_gas
	    if(IsNaN(DLconc(jesp)*0.d0)) then
	      print*,"Error of infinity/NaN initial",DLconc(jesp)
	      stop
	    endif
	    concentration_gas_loc(jesp) = DLconc(jesp)
            IF (concentration_gas_loc(jesp).LT.0.D0) THEN
	      concentration_gas_loc(jesp) = 0.D0
	    ENDIF
         ENDDO
         
	 total_nb=0.d0
!     Number concentration: loop on bins
         IF (with_number.EQ.1) THEN
            DO j=1, N_size
	      if(IsNaN(DLnum_conc_aer(j)*0.d0)) then
		print*,"Error of infinity/NaN initial",DLnum_conc_aer(j)
		stop
	      endif
               concentration_number(j) = DLnum_conc_aer(j)
               total_nb=total_nb+concentration_number(j)
! 	      if(concentration_number(j).gt.0.d0) then
! 	      print*,j,concentration_number(j)
! 	      endif
               !concentration_number(j) = DMAX1(concentration_number(j), TINYN)
               !why risk to add number into the system!?
            ENDDO
         ENDIF
         
	 total_ms=0.d0
	 total_aero_mass=0.d0
         DO j=1,N_size
            DO jesp=1,N_aerosol-1
	      if(IsNaN(DLconc_aer(j,jesp)*0.d0)) then
		print*,"Error of infinity/NaN initial",DLconc_aer(j,jesp)
		stop
	      endif
               concentration_mass(j,jesp) = DLconc_aer(j,jesp)
               total_aero_mass(jesp)=total_aero_mass(jesp)+concentration_mass(j,jesp)
               total_ms=total_ms+concentration_mass(j,jesp)
               !concentration_mass(j) = DMAX1(concentration_mass(j), TINYM)
               !why risk to add mass into the system!?
            ENDDO
         ENDDO

! 	 if(total_ms*total_nb.eq.0.d0.and.total_ms.ne.total_nb) then
! 	   Print*,"Bad input data mass=",total_ms,"number=",total_nb
! 	   stop
! 	 endif

         DO jesp=1,N_aerosol
	    IF(aerosol_species_interact(jesp).GT.0) THEN
	      concentration_gas(jesp)=concentration_gas_loc(aerosol_species_interact(jesp))
	    ENDIF
	    total_mass(jesp)=total_aero_mass(jesp)+concentration_gas(jesp)
! 	    if(total_mass(jesp).gt.1.d3.and.jesp.ne.N_aerosol) print*,"total_mass>100",&
! 	      jesp,total_mass(jesp),total_aero_mass(jesp),concentration_gas(jesp)

	    total_ms_old(jesp)= total_mass(jesp)
         ENDDO

        !print*,"finished concentration initialization"
  !!     Aerosol density converted in microg / microm^3.
        fixed_density= fixed_density_aer * 1.D-09!µg/µm3
        fixed_density_l=fixed_density_aer*1.D+09!µg/m3

        density_aer_bin=fixed_density
        density_aer_size=fixed_density
        !print*,"finished density"
        if(with_fixed_density.ne.1) then! modified
          !fixed density become overall averaged density
          !print*,"compute_all_density"
          call compute_all_density()
        endif
        call check_nan_inf(1)
  !!     Aerosol discretization converted into microm.
        do b= 1, (N_sizebin+1)
          diam_bound(b)=bin_bound_aer(b)* 1.D06
          if(b.eq.N_sizebin+1) then
            mass_bound(b)= density_aer_size(b-1)* cst_pi6 * diam_bound(b)**3
          else
            mass_bound(b)= density_aer_size(b)* cst_pi6 * diam_bound(b)**3
          endif
          log_bound(b)=dlog10(mass_bound(b))
        end do
        !print*,"problem of diameters"

        do b =1,N_sizebin
          size_diam_av(b)=dsqrt(diam_bound(b)*diam_bound(b+1))
          size_log_av(b)=(log_bound(b) + log_bound(b+1))*5.D-01
          size_sect(b)=log_bound(b+1) - log_bound(b)
        enddo

!	call mass_conservation(concentration_mass,concentration_number,&
!		concentration_gas, total_mass)	

	IF (with_number.NE.1) THEN
	  call compute_number()
        ENDIF    

	call compute_average_bin_diameter()
! 	call check_nan_inf(2)
!	IF (with_number.NE.1) THEN
!	  call computer_number()
!	ENDIF
! 	call check_nan_inf(3)
	do j =1,N_size
	  b=concentration_index(j,1)
	  cell_log_av(j)=size_log_av(b)
	enddo

	call compute_average_diameter()
! 	call check_nan_inf(4)

	call compute_wet_mass_diameter(1,N_size,concentration_mass,concentration_number,&
		concentration_inti,wet_mass,wet_diameter,wet_volume)

! 	call check_nan_inf(5)		

         IF (cloud_water.GE.LWCmin) THEN
            IF (with_incloud_scav.EQ.1) THEN
               rain_rate = DLrain
            ELSE
               rain_rate = 0.d0
            ENDIF           
            
            IF (aqueous_module.EQ.1) THEN
	       CALL VSRMCHEM(nesp_loc,N_aerosol,&
               N_size,N_sizebin,nesp_cloud_interact,&
               density_aer_bin,fixed_density,&
               diam_bound,cell_diam_av,log_bound,&
               concentration_gas_loc,concentration_mass,&
               Humidity,Pressure,Temperature,cloud_water,&
               tschem_aer,tfchem_aer,rain_rate,pH,&
               qscav_gas,qscav_aer,qscav_num,&
               with_number,concentration_number,DQLIMIT,&
               cloud_species_interact,concentration_index,&
               List_species,nesp_isorropia,nesp_aec,nesp_pankow,&
               nesp_pom,section_pass,redistribution_method,&
               with_fixed_density)

               DO jesp=1,N_aerosol
                  IF(aerosol_species_interact(jesp).GT.0) THEN
                     concentration_gas(jesp)=concentration_gas_loc(aerosol_species_interact(jesp))
                  ENDIF
               ENDDO !! added by YK: update of concentration_gas

	       !here solves the aqueous-phase model.
	       !compute qscav_gas ; qscav_aer ; qscav_num
               DO jesp=1,nesp_loc ! N_gas
                  Wet_Deposition(jesp)=Wet_Deposition(jesp)+&
                  qscav_gas(jesp) * layer_height /dtchem_aer
               ENDDO

               DO jesp=1,N_aerosol
                  DO j=1,N_size
                     Wet_Deposition_aer(j,jesp) =Wet_Deposition_aer(j,jesp) +&
                     qscav_aer(j,jesp) * layer_height / dtchem_aer
                  ENDDO
               ENDDO
               
               IF (with_number.EQ.1) THEN
                  
                  DO j=1,N_size
                     Wet_Deposition_Number_aer(j) =Wet_Deposition_Number_aer(j) +&
                     qscav_num(j) *layer_height / dtchem_aer
                  ENDDO
                
               ENDIF


            ELSE IF(aqueous_module.EQ.2) THEN ! Use simple aqueous module.

               CALL simple_aqueous_module(nesp_loc,N_aerosol,&
               N_size,N_sizebin,nesp_cloud_interact,&
               density_aer_bin,fixed_density,&
               diam_bound,cell_diam_av,log_bound,&
               concentration_gas_loc,concentration_mass,&
               Humidity,Pressure,Temperature,cloud_water,&
               tschem_aer,tfchem_aer,rain_rate,pH,&
               qscav_gas,qscav_aer,qscav_num,&
               with_number,concentration_number,DQLIMIT,&
               cloud_species_interact,concentration_index,&
               List_species,nesp_isorropia,nesp_aec,nesp_pankow,&
               nesp_pom,section_pass,redistribution_method,&
               with_fixed_density)

               DO jesp=1,N_aerosol
                  IF(aerosol_species_interact(jesp).GT.0) THEN
                     concentration_gas(jesp)=concentration_gas_loc(aerosol_species_interact(jesp))
                  ENDIF
               ENDDO !! added by YK: update of concentration_gas

               !print*,"finished simple_aqueous_module"
	       !here solves the aqueous-phase model.
	       !compute qscav_gas ; qscav_aer ; qscav_num
               DO jesp=1,nesp_loc ! N_gas
                  Wet_Deposition(jesp)=Wet_Deposition(jesp)+&
                  qscav_gas(jesp) * layer_height /dtchem_aer
               ENDDO
               
               DO jesp=1,N_aerosol
                  DO j=1,N_size
                     Wet_Deposition_aer(j,jesp) =Wet_Deposition_aer(j,jesp) +&
                     qscav_aer(j,jesp) * layer_height / dtchem_aer
                  ENDDO
               ENDDO

               IF (with_number.EQ.1) THEN

                  DO j=1,N_size
                     Wet_Deposition_Number_aer(j) =Wet_Deposition_Number_aer(j) +&
                     qscav_num(j) *layer_height / dtchem_aer
                  ENDDO
	      ENDIF

            ELSE

		CALL AERODYN(tschem_aer,tfchem_aer, &
       psoap_config, psurrogate)

               !start condensation coagulation and nucleation
            ENDIF

         ELSE

		CALL AERODYN(tschem_aer,tfchem_aer, &
       psoap_config, psurrogate)
               !start condensation coagulation and nucleation

         ENDIF

	! call mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass) ! YK

	  tot_tot=0.d0

! 	 OPEN(UNIT=10,ACCESS='APPEND',FILE="test-report.txt")
         DO jesp=1,nesp_aer-1
	    total_ms=total_mass(jesp)-total_ms_old(jesp)
	    tot_tot=tot_tot+total_aero_mass(jesp)
	    if((total_ms*total_ms).gt.(TINYN)) then
! 	      if(jesp.ne.1.and.jesp.ne.2.and.jesp.ne.11) &
	      print*,"warning! mass not conserved - species:",jesp,&
	      total_ms,total_aero_mass(jesp),concentration_gas(jesp)
	      !concentration_gas(jesp)=concentration_gas(jesp)-total_ms
	    endif
	    if(concentration_gas(jesp).lt.0.d0) then
	      print*,"warning! gas<0 species:",jesp,concentration_gas(jesp)
	      concentration_gas(jesp)=0.d0
	    endif
	    IF(aerosol_species_interact(jesp).GT.0) THEN
		concentration_gas_loc(aerosol_species_interact(jesp))=concentration_gas(jesp)
	    ENDIF
         ENDDO
!          CLOSE(10)
         DO j=1,nbin_aer
            DO jesp=1,nesp_aer
               DLconc_aer(j,jesp) =concentration_mass(j,jesp)
	      if(IsNaN(DLconc_aer(j,jesp)*0.d0)) then
		print*,"Error of infinity/NaN end",DLconc_aer(j,jesp)
		stop
	      endif               
            ENDDO
         ENDDO

	 total_nb=0.d0
!!     Number concentration
         IF(with_number.EQ.1) THEN
            DO j=1,nbin_aer
               DLnum_conc_aer(j) = concentration_number(j)
               total_nb=total_nb+concentration_number(j)
            ENDDO
         ENDIF
         
         
         DO jesp=1,nesp_loc
            DLconc(jesp) = concentration_gas_loc(jesp)
	    if(IsNaN(DLconc(jesp)*0.d0)) then
	      print*,"Error of infinity/NaN end",DLconc(jesp)
	      stop
	    endif
         ENDDO

         DO j=1,nbin_aer
	      if(IsNaN(DLnum_conc_aer(j)*0.d0)) then
		print*,"Error of infinity/NaN end",DLnum_conc_aer(j)
		stop
	      endif
         ENDDO
	
      ENDDO

      CALL free_allocated_memory()

      IF (with_coag.EQ.1) THEN
	call DeallocateCoefficientRepartition()
      ENDIF
      
      END SUBROUTINE aerosol

End module zAerosolSCRAM
