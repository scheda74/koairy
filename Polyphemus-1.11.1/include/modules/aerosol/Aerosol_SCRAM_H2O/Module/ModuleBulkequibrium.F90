!!-----------------------------------------------------------------------
!!     Copyright (C) 2003-2014, ENPC - INRIA - EDF R&D
!!     Author(s): Shupeng Zhu
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the air quality modeling system Polyphemus.
!!
!!     Polyphemus is developed in the INRIA - ENPC joint project-team
!!     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
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
!!
!!     -- DESCRIPTION
!!    This module contains methods to solve particle condensation/evaporation
!!    based on the gas/aerosol bulk equilibrium method.
!!-----------------------------------------------------------------------
Module iBulkequibrium
  use aInitialization
  use cThermodynamics
  implicit none
contains

  subroutine bulkequi_org(nesp_eq)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves gas/aerosol bulk equilibrium for organic species
!     based on H2O method.
!     equilibrium will be established between all size bin and bulk gas
!     mass flux will be redistributed into each cell based on rates
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!------------------------------------------------------------------------   
    implicit none
    
    integer::j,jesp,iter,s,k
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol)
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    double precision:: total_ms(N_aerosol)
    
!!     ******zero init
    do s=1,(N_aerosol-1)
	jesp=List_species(s)
	dq(jesp)=0.D0
	qgas(jesp)=0.d0
	qext(jesp)=0.D0
	qextold(jesp)=0.D0
	qaero(jesp)=0.d0
	ce_kernal_coef_tot(jesp)=0.D0!TINYA
	total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.01
  !     ****** allocate equilibrium species list      

    do s=1,nesp_aec
      eq_species(s)=aec_species(s)
    enddo
    do s=1,nesp_pankow
      eq_species(s+nesp_aec)=pankow_species(s)
    enddo
    do s=1,nesp_pom
      eq_species(s+nesp_aec+nesp_pankow)=poa_species(s)
    enddo
  !compute local equi

  !compute only c/e coefficients
    do s=1,(N_aerosol-1)
      jesp=List_species(s)
      if (aerosol_species_interact(jesp).GT.0) then
	rhop_tmp = 1400 ! kg/m3
	emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
		  
	do j = 1,N_size	!FOR OGANIC
! 	  if(concentration_index(j, 1).eq.1) then
! 	    wet_diam=0.01
! 	  else
 	 wet_diam=wet_diameter(j)
! 	  endif
	call COMPUTE_KELVIN_COEFFICIENT(&
		  Temperature,&          ! temperature (Kelvin)
		  emw_tmp,&       ! ext mol weight (g.mol-1)
		  surface_tension(jesp),&   ! surface tension (N.m-1) from INC
		  wet_diam,&         ! wet aero diameter (µm)
		  rhop_tmp,&      ! aerosol density (kg.m-3)
		  Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
	  call COMPUTE_CONDENSATION_TRANSFER_RATE(&
		diffusion_coef(jesp), &! diffusion coef (m2.s-1)
		quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
		accomodation_coefficient(jesp),& ! accomadation coef (adim)
		wet_diameter(j),   & ! wet aero diameter (Âµm)
		ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
	  !ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
		      !+ce_kernal_coef(j,jesp)*concentration_number(j)
! 	  if(abs(Kelvin_effect(j,jesp)-1.d0).lt.1.d-6) then
! 	    print *,"Kelvin exception",j,jesp,"Kelvin_effect=",Kelvin_effect(j,jesp),wet_diam
! 	    stop
! 	  endif
	  if(Kelvin_effect(j,jesp).lt.1.d0) Kelvin_effect(j,jesp)=1.01
	  ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
		      +ce_kernal_coef(j,jesp)*concentration_number(j)&
		      *(1.d0/(Kelvin_effect(j,jesp)-1.d0))

! 	    if (jesp.eq.30) then
! 	      print*,j,wet_diam,Kelvin_effect(j,jesp),(1.d0/(Kelvin_effect(j,jesp)-1.d0))
! 	      print*,concentration_number(j),concentration_number(j)*(1.d0/(Kelvin_effect(j,jesp)-1.d0))
! 	      print*,ce_kernal_coef(j,jesp),ce_kernal_coef(j,jesp)*concentration_number(j)&
! 			*(1.d0/(Kelvin_effect(j,jesp)-1.d0))
! 	    endif		      
	enddo
      endif
    enddo

    do s=1,(N_aerosol-1)
      jesp=List_species(s)
      do j=1,N_size
	qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
      enddo
      qgas(jesp)=concentration_gas(jesp)!initial gas for H2O
      total_ms(jesp)=qaero(jesp)+qgas(jesp)
!       if(IsNaN(qgas(jesp)*0.d0)) print*,"before",jesp,qgas(jesp)
      qextold(jesp)=qaero(jesp)
!       print*,"initial:",jesp, "aero",qaero(jesp),"gas",qgas(jesp)
    end do
    if(sulfate_computation.eq.0) then
	qgas(ESO4)=concentration_gas(ESO4)
    else
	qgas(ESO4)=0.d0
    endif
    qgas(EH2O)=0.0
    organion = 0.D0
    watorg = 0.D0
    proton = 0.D0
    do iter = 1,NITER_AEC_AQ
      call isoropia_drv(N_aerosol,&
	    qaero,qgas,organion, watorg, proton, lwc, Relative_Humidity, Temperature)

      call aec_drv(N_aerosol,&
	  0, qaero, qgas, proton, lwc, organion,watorg, Relative_Humidity, Temperature,&
	  with_oligomerization, thermodynamic_model)
    enddo

    qaero(N_aerosol)=lwc
    do iter = 1,NITER_AEC_DRY
      call poa_drv(N_aerosol,qaero, qgas,soa_part_coef,Temperature,vaporization_enthalpy)

      call pankow_drv(N_aerosol,qaero, qgas,soa_part_coef)

      call aec_drv(N_aerosol,&
	1, qaero, qgas, proton, lwc, organion,watorg, Relative_Humidity, Temperature,&
	with_oligomerization, thermodynamic_model)

    enddo
    qgas(EH2O)=0.0
    
      !for organic
    do s=1,nesp_eq
      jesp=eq_species(s)
      if(aerosol_species_interact(jesp).GT.0) then
	d_ms=qaero(jesp)+qgas(jesp)-total_ms(jesp)
	qgas(jesp)=qgas(jesp)-d_ms
	if(qgas(jesp).lt.0.d0) then
	  qgas(jesp)=0.d0
	  qaero(jesp)=total_ms(jesp)
	endif
	if (IsNaN(qaero(jesp))) qaero(jesp)=0.d0
	qext(jesp)=qaero(jesp)
	!print*,"aero",jesp,qaero(jesp),qgas(jesp)
	concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
	if(qext(jesp).gt.0.d0) then
	  dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
	else
	  dq(jesp)=-qextold(jesp)
	endif
      endif
    enddo

!     ******redistribute on each cell according to Rates
    call bulkequi_redistribution_anck(concentration_number,concentration_mass,&
    nesp_eq,eq_species,N_size,dq,ce_kernal_coef,ce_kernal_coef_tot,Kelvin_effect)

  end subroutine bulkequi_org

  subroutine bulkequi_inorg(nesp_eq)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves gas/aerosol bulk equilibrium for inorganic species
!     based on ISORROPIA method.
!     equilibrium will be established between all size bin < ICUT and bulk gas
!     mass flux will be redistributed into each cell based on rates
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!------------------------------------------------------------------------   
    implicit none
    integer::j,jesp,s
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::aatoteq,qgasa,qgasi
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam    
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    double precision:: total_ms(N_aerosol)
    
!!     ******zero init
    do s=1,(N_aerosol-1)
      jesp=List_species(s)
      dq(jesp)=0.D0
      qgas(jesp)=0.d0
      qext(jesp)=0.D0
      qextold(jesp)=0.D0
      qaero(jesp)=0.d0
      ce_kernal_coef_tot(jesp)=0.D0
      total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.01
  !     ****** if sulfate computed dynamically avoid it
  !     ****** allocate equilibrium species list
    do jesp=1,nesp_eq
      eq_species(jesp)=isorropia_species(nesp_isorropia-nesp_eq+jesp)
    enddo

  !compute only c/e coefficients

    ce_kernal_coef=0.d0

    do s=1, nesp_isorropia
      jesp=isorropia_species(s)
      if (aerosol_species_interact(jesp).GT.0) then
	    ! compute total ce_kernal_coef coef (s-1) EXCEPT SO4	
	if(jesp.ne.isorropia_species(2)) then
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
!          IF (jesp.NE.ECl) THEN
! #endif
	do j = 1,ICUT	!FOR INOGANIC
	  rhop_tmp = 1400 ! kg/m3
	  emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol

! 	  if(concentration_index(j, 1).eq.1) then
! 	    wet_diam=0.01
! 	  else
 	  wet_diam=wet_diameter(j)
! 	  endif
      if (wet_diam .lt. 1.d-3) then
         write(*,*) "bulkequi_inorg: too small wet_diameter",wet_diam
         stop
      endif !! YK
	  call COMPUTE_KELVIN_COEFFICIENT(&
		    Temperature,&          ! temperature (Kelvin)
		    emw_tmp,&       ! ext mol weight (g.mol-1)
		    surface_tension(jesp),&   ! surface tension (N.m-1) from INC
		    wet_diam,&         ! wet aero diameter (µm)
		    rhop_tmp,&      ! aerosol density (kg.m-3)
		    Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
	  call COMPUTE_CONDENSATION_TRANSFER_RATE(&
		diffusion_coef(jesp), &! diffusion coef (m2.s-1)
		quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
		accomodation_coefficient(jesp),& ! accomadation coef (adim)
		wet_diameter(j),   & ! wet aero diameter (Âµm)
		ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
	    if(Kelvin_effect(j,jesp).lt.1.d0) Kelvin_effect(j,jesp)=1.01
	    ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
			+ce_kernal_coef(j,jesp)*concentration_number(j)&
			*(1.d0/(Kelvin_effect(j,jesp)-1.d0))
	enddo
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
!          ENDIF
! #endif
	endif
      endIF
    enddo

  !     ****** if sulfate computed dynamically avoid it
    if (sulfate_computation.eq.0) then
      jesp=isorropia_species(2)!
      do j=1,N_size
	ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&
	    +ce_kernal_coef(j,jesp)*concentration_number(j)
      end do
      qgasa=0.d0
      qgasi=concentration_gas(jesp)
    else
      qgasi=0.d0
    endif

    !compute total mass of each species
    do s=1,nesp_isorropia!inorganic
      jesp=isorropia_species(s)
      do j=1,ICUT
	qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
      enddo
      qgas(jesp)=concentration_gas(jesp)
      total_ms(jesp)=qaero(jesp)+qgas(jesp)
      qextold(jesp)=qaero(jesp)
    end do

    if(sulfate_computation.eq.0) then
      !compute apparent gas concentration of sulfate
      jesp=isorropia_species(2)!
      do j=1,ICUT
	aatoteq=aatoteq+ce_kernal_coef(j,jesp)*concentration_number(j)
      enddo
      if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	qgas(jesp)=concentration_gas(jesp)*aatoteq/ce_kernal_coef_tot(jesp)
      endif
      qgasa=qgas(jesp)
    else
      qgas(ESO4)=0.d0
    endif

! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
!       qaero(ENa) = 0.D0
!       qaero(ECl) = 0.D0
!       qgas(ECl) = 0.D0
! #endif

    organion= 0.D0
    watorg = 0.D0
    proton= 0.D0

  !eqilibirum for inorganic
    call isoropia_drv(N_aerosol,&
	  qaero,qgas,organion, watorg, proton, lwc, Relative_Humidity, Temperature)
    qaero(N_aerosol)=lwc
    qgas(EH2O)=0.0

  !     ******redistribute on each cell according to Rates

    if(sulfate_computation.eq.0) then
      ce_kernal_coef_tot(ESO4)=aatoteq!for later redistribution
    endif
    do s=1, nesp_eq
      jesp=eq_species(s)
      if(aerosol_species_interact(jesp).GT.0) then      
	d_ms=qaero(jesp)+qgas(jesp)-total_ms(jesp)
	qgas(jesp)=qgas(jesp)-d_ms
	if(qgas(jesp).lt.0.d0) then
	  qgas(jesp)=0.d0
	  qaero(jesp)=total_ms(jesp)
	endif	
	if (IsNaN(qaero(jesp))) qaero(jesp)=0.d0!detect NaN case
	qext(jesp)=qaero(jesp)
	!print*,"aero",jesp,qaero(jesp),qgas(jesp)
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
! 	if(jesp.ne.ECl) then
! #endif
 	concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
! 	endif
! #endif	
	if(qext(jesp).ge.0.d0) then
	  dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
	else
	  print*,"ModBulk inorg qext(jesp)<0!",jesp,qext(jesp)
	  stop
	endif
      endif
    enddo
    
    call bulkequi_redistribution_anck(concentration_number,concentration_mass,&
    nesp_eq,eq_species,N_size,dq,ce_kernal_coef,ce_kernal_coef_tot,Kelvin_effect)
  ! give back initial SO4 gas conc
  ! minus that consumed by equi bins
    if(sulfate_computation.eq.0) then
      concentration_gas(ESO4)=qgasi-qgasa
    else
      concentration_gas(ESO4)=0.d0
    endif
    
  end subroutine bulkequi_inorg

   subroutine bulkequi_redistribution_anck(c_number,c_mass,nesp_eq,&
  eq_species,end_bin,dq,AAi,ce_kernal_coef_tot,ck_ef)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute bulk mass variations into each bins (using mass instead of number)
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!     nesp_eq: number of species at equilibirum
!     eq_species: the list of species pointer at equilibirum
!     end_bin: marks the number of bins concerned during the redistribution
!     dq: bulk mass variations
!     AAi: c/e kernel coefficient          ([m3.s-1]).
!     ce_kernal_coef_tot: sum of c/e kernel coefficient          ([m3.s-1]).
!------------------------------------------------------------------------
    implicit none
    integer::j,s,jesp,end_bin,iclip
    integer::nesp_eq
    integer::eq_species(nesp_eq)
    double precision::totaer,temp_mass,totaa
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::dq(N_aerosol)
    double precision::AAi(N_size,N_aerosol)
    double precision::frac(N_size,N_aerosol)
    double precision::c_number(N_size)
    double precision::c_mass(N_size,N_aerosol)
    double precision::ck_ef(N_size,N_aerosol)
    double precision::sum_frac
    integer::k

    frac=0.d0

    do s=1, nesp_eq
      jesp=eq_species(s)
      !print*,"bulk_rdb",jesp,dq(jesp)
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
!       IF (jesp.NE.ECl) THEN
! #endif
      iclip=0
      sum_frac=0.d0
      do j=1,end_bin!judgment
	  if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	    frac(j,jesp)= AAi(j,jesp)*c_number(j)*&
		(1.d0/(ck_ef(j,jesp)-1.d0))/ce_kernal_coef_tot(jesp)
	  endif
	  sum_frac=sum_frac+frac(j,jesp)
	  temp_mass=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	  if(temp_mass.lt.0.d0) then
	    iclip=1!case of over evaporation
	    !print*,"overevaporation!",j,jesp,c_mass(j,jesp),dq(jesp),frac(j,jesp)
	  endif
      enddo
      !print*, "dq(",jesp,")=",dq(jesp)
      !print *,"sum_fra=",sum_frac
!       if(ce_kernal_coef_tot(jesp).gt.0.d0.and.abs(sum_frac-1).gt.1.d-15) then
! 	print *,"wait exception",jesp,"sum_wait=",sum_frac,ce_kernal_coef_tot(jesp)
! 	stop
!       endif
      totaer=0.d0
      totaa=0.d0
      do j=1,end_bin
	totaer=totaer+c_mass(j,jesp)
	totaa=totaa+AAi(j,jesp)
      enddo
      if(iclip.eq.1) then !over evaporate
	if(totaer.gt.0.d0) then
	  do j=1,end_bin
	    frac(j,jesp)=c_mass(j,jesp)/totaer
	    c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	    !print*,"overevapor2!",j,jesp,c_mass(j,jesp),dq(jesp),frac(j,jesp)
	  enddo
	endif
      else!normal case
	do j=1,end_bin
	  if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	    frac(j,jesp)= AAi(j,jesp)*c_number(j)*&
		(1.d0/(ck_ef(j,jesp)-1.d0))/ce_kernal_coef_tot(jesp)
	  else
	    frac(j,jesp)=0.0
	    if (aerosol_species_interact(jesp).GT.0) frac(j,jesp)=AAi(j,jesp)/totaa
	  endif
	  c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
! 	  if(jesp.eq.9) then
! 	    print*,j,frac(j,jesp),c_mass(j,jesp), dq(jesp),frac(j,jesp),"/",ce_kernal_coef_tot(jesp)
! 	  endif
	enddo
      endif
      !check mass leek
      temp_mass=0.d0
!       OPEN(UNIT=10,ACCESS='APPEND',FILE="test-report.txt")
      do j=1,end_bin
! 	write(unit=10,FMT=*),j,jesp,c_mass(j,jesp),dq(jesp),frac(j,jesp),dq(jesp)*frac(j,jesp)
	temp_mass=temp_mass+c_mass(j,jesp)
      enddo
!       CLOSE(10)
      if(dq(jesp).ne.temp_mass-totaer) &
	concentration_gas(jesp)=concentration_gas(jesp)-(temp_mass-totaer-dq(jesp))
! #ifdef WITHOUT_NACL_IN_THERMODYNAMICS
!       ENDIF
! #endif
      !print*,"tot_mass=",temp_mass
    enddo

  end subroutine bulkequi_redistribution_anck  
  
end Module iBulkequibrium
  
  
