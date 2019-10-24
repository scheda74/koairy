!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Shupeng Zhu
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the SSH-aerosol model.
!!
!!     SSH-aerosol is a free software; you can redistribute it and/or modify
!!     it under the terms of the GNU General Public License as published
!!     by the Free Software Foundation; either version 2 of the License,
!!     or (at your option) any later version.
!!
!!     SSH-aerosol is distributed in the hope that it will be useful, but
!!     WITHOUT ANY WARRANTY; without even the implied warranty of
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!!     General Public License for more details.
!!
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module read configuration file and initialize all global variables
!!-----------------------------------------------------------------------
module aInitialization
    implicit none
    INCLUDE 'CONST.INC'
    INCLUDE 'CONST_A.INC'

    !!part 1: parameters of system dimension
    Integer :: N_gas   !complete gas species number
    integer :: N_size   ! total number of size and composition sections
    integer :: N_groups!Number of groups
    integer :: N_fracmax! maximum number of composition sections per size section
    integer :: N_aerosol !Number of aerosol species
    integer :: N_sizebin!number of  size sections
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(N_gas,N_size,N_groups,N_fracmax,N_aerosol,N_sizebin)
#endif
    integer :: N_organics !Number of organics aerosol species
    integer :: N_inorganic!Number of inorganic aerosol species
    integer :: N_inert!number of inert aerosol species
    integer :: N_liquid!Number of liquid internal species
    integer :: N_solid!Number of solid internal species
    integer :: N_inside_aer!Number of internal species
    integer :: N_hydrophilic!Number of hydrophilic organics aerosol species
    parameter (N_organics=23,N_inorganic=5,N_inert=2,N_liquid=12)
    parameter (N_solid=9,N_inside_aer=21)
    parameter(N_hydrophilic=9)

    !!part 2: parameters of system option    
    integer :: tag_thrm!method for wet diameter computation0 h2o 1 isorropia
    integer :: dynamic_solver  !KDSLV Tag type of solver
    integer :: sulfate_computation !ISULFCOND tag of sulfate condensation method
    integer :: redistribution_method !tag of redistribution method
    integer :: with_coag !Tag gCoagulation
    integer :: with_cond !Tag fCondensation
    integer :: with_nucl !Tag nucleation
    Integer :: aqueous_module!ICLD
    Integer :: with_incloud_scav!IINCLD
    Integer :: with_kelvin_effect!IKELV
    Integer :: with_fixed_density!IDENS
    integer :: ICUT!cutting_bin
    integer :: section_pass
    Integer :: nucl_model!ITERN
    Integer :: wet_diam_estimation!ITHRM
    Integer :: with_oligomerization!IOLIGO
!    Integer :: thermodynamic_model!ITHERMO
    Integer :: ISOAPDYN
    Integer :: with_number!INUM
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(tag_thrm,dynamic_solver,sulfate_computation,redistribution_method)
!$omp threadprivate(with_coag,with_cond,with_nucl,aqueous_module,with_incloud_scav)
!$omp threadprivate(with_kelvin_effect,with_fixed_density,ICUT,nucl_model,wet_diam_estimation)
!$omp threadprivate(section_pass,with_oligomerization,thermodynamic_model,with_number)
#endif
    integer NITER_PKW,NITER_AEC_AQ
    integer NITER_POA,NITER_AEC_DRY    
    double precision ::  ALFHP! percentage of H+ allowed to c/e(0.1)
    double precision ::  EPSER
    double precision ::  TINYM,TINYN,MTSBL
    double precision ::  DTMAX
    double precision :: DTAEROMIN ! minimum timestep for aerosols (in s)
    double precision ::  DMIN,DMAX
    parameter(DMIN = 1.D-3)
    parameter(ALFHP = 0.1D0)
    parameter(TINYM = 1.D-20)
    parameter(TINYN = 1.D-15)
    parameter(DTMAX =10.D0)
    parameter(NITER_AEC_AQ = 1)
    parameter(NITER_AEC_DRY = 1)
    parameter(NITER_PKW = 5)
    parameter(NITER_POA = 10)
    parameter(MTSBL = 1.0D0)

    !!part 3: System pointers
    Integer :: E1,E2,G1,G2!Mark the begin and end of dynamic aerosol (except EH2O)
    Integer :: ENa,ESO4,ENH4,ENO3,ECl,EMD,EBC,EH2O!inorganic pointers
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(E1,E2,G1,G2,ENa,ESO4,ENH4,ENO3,ECl,EMD,EBC,EH2O)
#endif
    Integer :: ictmNH3,ictmHNO3,ictmHCl,ictmSO2,ictmH2O2,ictmHCHO,ictmHNO2
    Integer :: ictmO3,ictmOH,ictmHO2,ictmNO3,ictmNO,ictmNO2,ictmPAN,ictmH2SO4!pointers of cloud species.
    Integer :: IH,INa,INH4,ICl,ISO4,IHSO4,INO3,IH2O,INH3,IHCl,IHNO3,IOH
    Integer :: SNaNO3,SNH4NO3,SNACl,SNH4Cl,SLC,SNa2SO4,SNH42S4,SNaHSO4,SNH4HS4
    Parameter (IH=1,INa=2,INH4=3,ICl=4,ISO4=5,IHSO4=6,INO3=7,IH2O=8,INH3=9)
    Parameter (IHCl=10,IHNO3=11,IOH=12,SNaNO3=13,SNH4NO3=14,SNACl=15)
    Parameter (SNH4Cl=16,SNa2SO4=17,SNH42S4=18,SNaHSO4=19,SNH4HS4=20,SLC=21)
    
    Integer :: nesp, nesp_isorropia, nesp_aec, nesp_pankow, nesp_pom!Number of different species group
    Integer, dimension(:), allocatable :: isorropia_species
    Integer, dimension(:), allocatable :: aec_species
    Integer, dimension(:), allocatable :: pankow_species
    Integer, dimension(:), allocatable :: poa_species
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(nesp, nesp_isorropia, nesp_aec, nesp_pankow, nesp_pom)
#endif

    !!part 4: System state parameters    
    integer :: tagrho
    integer :: tag_coag,tag_cond,tag_nucl
    integer :: kind_composition
    double precision :: timestep_splitting,sub_timestep_splitting
    double precision :: final_time,dtmin! Time step and finnal time
    double precision :: initial_time_splitting,current_sub_time,final_sub_time!current time=initial_time_splitting+current_sub_time
    double precision :: Temperature,Relative_Humidity,Pressure,Humidity
    double precision :: fixed_density,fixed_density_l!density of overall partical
    double precision :: Cut_dim!cuting diameter between equi/dynamic
    double precision :: viscosity!Dynamic viscosity ([kg/m/s]).
    double precision :: air_free_mean_path
    double precision :: total_water!total mass of water
    double precision :: total_IH!total mass of H+
    double precision :: total_PH!overall PH value
    double precision :: n_grow_nucl,n_grow_coag,n_emis
    double precision :: m_grow_cond,m_emis
    double precision :: total_number,o_total_mass,total_mass_t
    double precision :: record_time
    Double precision :: p_fact,k_fact!??
    Double precision :: DQLIMIT
    
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(tagrho,tag_coag,tag_cond,tag_nucl,kind_composition,timestep_splitting)
!$omp threadprivate(sub_timestep_splitting,final_time,dtmin,initial_time_splitting,current_sub_time)
!$omp threadprivate(final_sub_time,Temperature,Relative_Humidity,Pressure,Humidity)
!$omp threadprivate(fixed_density,fixed_density_l,Cut_dim,viscosity,air_free_mean_path)
!$omp threadprivate(total_water,total_IH,total_PH,n_grow_nucl,n_grow_coag,n_emis)
!$omp threadprivate(m_grow_cond,m_emis,total_number,o_total_mass,total_mass_t)
!$omp threadprivate(DQLIMIT,record_time,p_fact,k_fact)
#endif

    !!part5: 1 dimension data array    
    integer, dimension(:), allocatable :: Index_groups!index of which group the species belongs to
    integer, dimension(:), allocatable :: List_species!read species defined in cfg files
    Integer, dimension(:), allocatable :: aerosol_species_interact
    Double precision,dimension(:), allocatable :: density_aer_bin !density of each grid bins
    Double precision,dimension(:), allocatable :: density_aer_size !density of each size section
    Double precision,dimension(:), allocatable :: diam_bound! DBF diameter bounds of each size section
    Double precision,dimension(:), allocatable :: mass_bound! MBF
    Double precision,dimension(:), allocatable :: log_bound!XBF
    Double precision,dimension(:), allocatable :: total_bin_mass!total mass of each size section
    Double precision,dimension(:), allocatable :: size_sect!HSF log size of each section
    Double precision,dimension(:), allocatable :: size_diam_av!DSF average diameter of each size section
    Double precision,dimension(:), allocatable :: size_mass_av!MSF average mass of each size section
    Double precision,dimension(:), allocatable :: size_log_av!XSF
    Double precision,dimension(:), allocatable :: cell_diam_av!!DSF average diameter of each grid cell
    Double precision,dimension(:), allocatable :: cell_mass_av!!MSF average mass of each grid cell
    Double precision,dimension(:), allocatable :: cell_log_av!XSF
    Double precision,dimension(:), allocatable :: total_mass!total mass of each species
    Double precision,dimension(:), allocatable :: mass_total_grid!total mass of each grid cell
    Double precision,dimension(:), allocatable :: total_aero_mass!total aerosol mass of each species
    Double precision,dimension(:), allocatable :: bin_mass!mass concentration of each size section
    Double precision,dimension(:), allocatable :: bin_number!number concentration of each size section
    Double precision,dimension(:), allocatable :: concentration_number_tmp!first order approximation of number
    Double precision,dimension(:), allocatable :: concentration_number!number concentration of each grid cell
    Double precision,dimension(:), allocatable :: concentration_gas! gas concentration of each species
    Double precision,dimension(:), allocatable :: wet_diameter!Aerosol wet diameter (µm). of each grid cell
    Double precision,dimension(:), allocatable :: wet_mass!Aerosol wet mass (µg). of each grid cell
    Double precision,dimension(:), allocatable :: wet_volume!Aerosol wet volume (µm^3). of each grid cell
    Double precision , dimension(:), allocatable :: rho_wet_cell
    Double precision , dimension(:), allocatable :: cell_mass

    !!part6: 2+ dimension data array
    integer, dimension(:,:), allocatable :: concentration_index !matrix from grid index to size and composition index
    integer, dimension(:,:), allocatable :: concentration_index_iv !matrix from size and composition to grid index
    double precision , dimension(:,:), allocatable :: kernel_coagulation
    double precision , dimension(:,:), allocatable :: ce_kernal_coef!c/e kernal
    double precision , dimension(:,:), allocatable :: Kelvin_effect_ext!kelvin effect
    double precision , dimension(:,:), allocatable :: frac_grid !excat fraction of each species in each grid
    double precision , dimension(:,:), allocatable :: concentration_mass
    double precision , dimension(:,:), allocatable :: concentration_mass_tmp!first order apporximation
    double precision , dimension(:,:), allocatable :: concentration_inti!internal inorganic aerosol concentration ([ï¿½g.m-3]).
    double precision , dimension(:,:), allocatable :: dqdt
    double precision , dimension(:,:,:), allocatable :: discretization_composition! multi-array storing discretization of composition

    !! part 7: basic physical and chemical parameters
    double precision :: mass_density_solid(SNaNO3:SLC)!molar weight of internal solids species
    double precision :: molecular_weight_inside(N_liquid)!molar weight of inorganic species in aqueous_phase
    double precision :: molecular_weight_solid(SNaNO3:SLC)!molar weight of solids
    double precision ,dimension(:), allocatable :: saturation_pressure
    double precision ,dimension(:), allocatable :: partition_coefficient
    double precision ,dimension(:), allocatable :: vaporization_enthalpy
    double precision ,dimension(:), allocatable :: accomodation_coefficient
    double precision ,dimension(:), allocatable :: surface_tension
    double precision ,dimension(:), allocatable :: saturation_pressure_mass
    double precision ,dimension(:), allocatable :: saturation_pressure_torr
    double precision ,dimension(:), allocatable :: deliquescence_relative_humidity
    double precision ,dimension(:), allocatable :: molecular_weight_aer! (µg/mol)
    double precision ,dimension(:), allocatable :: molecular_diameter
    double precision ,dimension(:), allocatable :: collision_factor_aer
    double precision ,dimension(:), allocatable :: mass_density_aer!(µg/m3) liquid mass density
    double precision ,dimension(:), allocatable :: quadratic_speed! (m.s-1)
    double precision ,dimension(:), allocatable :: diffusion_coef! (m2.s-1)
    double precision ,dimension(:), allocatable :: soa_sat_conc! (µg.m-3)
    double precision ,dimension(:), allocatable :: soa_part_coef!(m3/microg)
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(mass_density_solid,molecular_weight_inside,molecular_weight_solid)
!$omp threadprivate(isorropia_species,aec_species,pankow_species,poa_species)
!$omp threadprivate(Index_groups,List_species,aerosol_species_interact,density_aer_bin)
!$omp threadprivate(density_aer_size,diam_bound,mass_bound,log_bound,total_bin_mass)
!$omp threadprivate(size_sect,size_diam_av,size_mass_av,size_log_av,cell_diam_av)
!$omp threadprivate(cell_mass_av,cell_log_av,total_mass,mass_total_grid,cell_mass)
!$omp threadprivate(total_aero_mass,bin_mass,bin_number,concentration_number_tmp,rho_wet_cell)
!$omp threadprivate(concentration_number,concentration_gas,wet_diameter,wet_mass,wet_volume)
!$omp threadprivate(concentration_index,concentration_index_iv,kernel_coagulation,ce_kernal_coef)
!$omp threadprivate(Kelvin_effect_ext,frac_grid,concentration_mass,concentration_mass_tmp)
!$omp threadprivate(concentration_inti,dqdt,discretization_composition)
!$omp threadprivate(saturation_pressure,partition_coefficient,vaporization_enthalpy,accomodation_coefficient)
!$omp threadprivate(surface_tension,saturation_pressure_mass,saturation_pressure_torr,molecular_weight_aer)
!$omp threadprivate(deliquescence_relative_humidity,molecular_diameter,collision_factor_aer,mass_density_aer)
!$omp threadprivate(quadratic_speed,diffusion_coef,soa_sat_conc,soa_part_coef)
#endif

 contains

    subroutine Init_global_parameters(nesp_loc,nbin_aer,nsize_section_aer,nesp_aer,&
				     ncomposition_aer,ngroup_aer,&
				     aerosol_species_interact_loc,&
				     nesp_isorropia_loc,isorropia_species_loc,&
				     nesp_aec_loc,aec_species_loc,&
				     nesp_pankow_loc,pankow_species_loc,&
				     nesp_pom_loc,poa_species_loc,&
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
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine initialize all global variables (Only used in 3D model)
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------
    integer b,jesp,j,f,k
    integer :: nesp_aer,nbin_aer,nsize_section_aer,nesp_loc
    integer :: aerosol_species_interact_loc(nesp_aer)
    integer :: nesp_isorropia_loc
    integer :: ncomposition_aer,ngroup_aer
    integer :: isorropia_species_loc(nesp_isorropia_loc)
    integer :: nesp_aec_loc,aec_species_loc(nesp_aec_loc)
    integer :: nesp_pankow_loc,pankow_species_loc(nesp_pankow_loc)
    integer :: nesp_pom_loc,poa_species_loc(nesp_pom_loc)
    integer :: md_species_loc,bc_species_loc
    integer :: nesp_cloud_interact
    integer :: cloud_species_interact(nesp_cloud_interact)
    integer :: noptions_aer,options_aer(noptions_aer)
    integer :: aerosol_species_group_relation(nesp_aer)
    double precision :: saturation_pressure_loc(nesp_aer)
    double precision :: partition_coefficient_loc(nesp_aer)
    double precision :: vaporization_enthalpy_loc(nesp_aer)
    double precision :: accomodation_coefficient_loc(nesp_aer)
    double precision :: surface_tension_loc(nesp_aer)
    double precision :: saturation_pressure_mass_loc(nesp_aer)
    double precision :: saturation_pressure_torr_loc(nesp_aer)
    double precision :: deliquescence_relative_humidity_loc(nesp_aer)
    double precision :: molecular_weight_aer_loc(nesp_aer)
    double precision :: molecular_diameter_aer_loc(nesp_aer)
    double precision :: collision_factor_aer_loc(nesp_aer)
    double precision :: mass_density_aer_loc(nesp_aer)

!! basic parameters (system dimensions)
    N_gas=nesp_loc
    nesp=nesp_loc
    N_aerosol=nesp_aer
    N_size=nbin_aer
    N_fracmax=ncomposition_aer
    N_groups=ngroup_aer
    N_sizebin=nsize_section_aer

    allocate(total_mass(N_aerosol))
    allocate(concentration_gas(N_aerosol))
    allocate(concentration_number(N_size))
    allocate(concentration_number_tmp(N_size))
    allocate(concentration_mass(N_size,N_aerosol))
    allocate(concentration_mass_tmp(N_size,N_aerosol))
    allocate(concentration_inti(N_size,N_inside_aer))
    allocate(cell_mass_av(N_size))
    allocate(cell_mass(N_size))
    allocate(cell_diam_av(N_size))
    allocate(cell_log_av(N_size))
    allocate(total_aero_mass(N_aerosol))
    allocate(mass_total_grid(N_size))
    allocate(Index_groups(N_aerosol-1))
    allocate(wet_mass(N_size))
    allocate(wet_diameter(N_size))
    allocate(wet_volume(N_size))
    allocate(size_diam_av(N_sizebin))
    allocate(size_mass_av(N_sizebin))
    allocate(size_log_av(N_sizebin))
    allocate(size_sect(N_sizebin))
    allocate(diam_bound(N_sizebin+ 1))
    allocate(mass_bound(N_sizebin+ 1))
    allocate(log_bound(N_sizebin+ 1))
    allocate(density_aer_bin(N_size))
    allocate(density_aer_size(N_sizebin))
    allocate(bin_mass(N_sizebin))
    allocate(bin_number(N_sizebin))
    allocate(rho_wet_cell(N_size))
    allocate(kernel_coagulation(N_size,N_size))
    allocate(total_bin_mass(N_sizebin))
    allocate(ce_kernal_coef(N_size,N_aerosol))
    allocate(Kelvin_effect_ext(N_size,N_aerosol))
!    allocate(frac_grid(N_size,N_aerosol)) ! YK
    allocate(frac_grid(N_size,N_groups))
    allocate(dqdt(N_size,N_aerosol))
    
    total_mass=0.d0
    concentration_gas=0.d0
    concentration_number=0.d0
    concentration_number_tmp=0.d0
    concentration_mass=0.d0
    concentration_mass_tmp=0.d0
    concentration_inti=0.d0
    total_aero_mass=0.d0
    cell_mass_av=0.d0
    cell_diam_av=0.d0
    cell_log_av=0.d0
    total_number=0.d0
    mass_total_grid=0.d0
    wet_mass=0.d0
    wet_diameter=0.d0
    wet_volume=0.d0
    size_diam_av=0.d0
    size_mass_av=0.d0
    size_log_av=0.d0
    size_sect=0.d0
    diam_bound=0.d0
    mass_bound=0.d0
    log_bound=0.d0
    density_aer_bin=0.d0
    density_aer_size=0.d0
 
    !!read relations between species and chemical groups
    do jesp =1,(N_aerosol-1)
      Index_groups(jesp)=aerosol_species_group_relation(jesp)+1
      !!c++ to frotran
    enddo

    !!set index relation between bins and size&composition
    allocate(concentration_index_iv(N_sizebin,N_fracmax))
    allocate(concentration_index(N_size,2))
    do b = 1,N_sizebin
      do f = 1,N_fracmax
	j=N_fracmax*(b-1)+f
	concentration_index_iv(b,f)=j
	concentration_index(j,1)=b
	concentration_index(j,2)=f
      enddo
     enddo

!!     Set parameters
    with_coag = options_aer(1)
    with_cond = options_aer(2)
    with_nucl = options_aer(3)
    aqueous_module = options_aer(4)
    with_incloud_scav = options_aer(5)
    with_kelvin_effect = options_aer(6)
    if (options_aer(7) .eq. 1) then
       with_fixed_density = 0
    else
       with_fixed_density = 1
    endif !! Fixed a bug :: YK
    ! with_fixed_density = options_aer(7)
    ICUT = options_aer(8)
    sulfate_computation = options_aer(9)
    dynamic_solver = options_aer(10)
    redistribution_method = options_aer(11)
    nucl_model = options_aer(12)
    wet_diam_estimation = options_aer(13)
    with_oligomerization = options_aer(14)
!    thermodynamic_model = options_aer(15)
    ISOAPDYN = options_aer(15)
    with_number = options_aer(16)
    !     Species to be computed
    IF(ICUT.GT.0.D0) THEN
      section_pass=concentration_index(ICUT,1)
    ELSE
      section_pass=1
    ENDIF
    
    E1=1
    E2=nesp_aer-1 ! to avoid water

    nesp_isorropia=nesp_isorropia_loc
    nesp_aec=nesp_aec_loc
    nesp_pankow=nesp_pankow_loc
    nesp_pom=nesp_pom_loc

    !!initialize pointers
    allocate(isorropia_species(nesp_isorropia_loc))
    allocate(aec_species(nesp_aec_loc))
    allocate(pankow_species(nesp_pankow_loc))
    allocate(poa_species(nesp_pom))
    allocate(aerosol_species_interact(nesp_aer))
    allocate(List_species(N_aerosol))

    ENa=isorropia_species_loc(1)
    ESO4=isorropia_species_loc(2)
    ENH4=isorropia_species_loc(3)
    ENO3=isorropia_species_loc(4)
    ECl=isorropia_species_loc(5)
    EMD=md_species_loc
    EBC=bc_species_loc
    EH2O=nesp_aer ! water always at the end
    
    !initialized the aerosol species pointer list without water
    List_species(1)=EMD
    List_species(2)=EBC

    do jesp=1,nesp_isorropia
      isorropia_species(jesp) = isorropia_species_loc(jesp)
      List_species(2+jesp)=isorropia_species_loc(jesp)
    enddo

    do jesp=1,nesp_aec
      aec_species(jesp) = aec_species_loc(jesp)
      List_species(2+nesp_isorropia+jesp)=aec_species_loc(jesp)
    enddo

    do jesp=1,nesp_pankow
      pankow_species(jesp) = pankow_species_loc(jesp)
      List_species(2+nesp_isorropia+nesp_aec+jesp)=pankow_species_loc(jesp)
    enddo

    do jesp=1,nesp_pom
      poa_species(jesp) = poa_species_loc(jesp)
      List_species(2+nesp_isorropia+nesp_aec+nesp_pankow+jesp)=poa_species_loc(jesp)	
    enddo

    List_species(N_aerosol)=EH2O
    
    G1=ESO4
    G2=ECl

    do jesp=1,nesp_aer
	aerosol_species_interact(jesp)=aerosol_species_interact_loc(jesp)
    enddo

    !! Pointers for cloud species.

    ictmNH3=cloud_species_interact(1)
    ictmHNO3=cloud_species_interact(2)
    ictmHCl=cloud_species_interact(3)
    ictmSO2=cloud_species_interact(4)
    ictmH2O2=cloud_species_interact(5)
    ictmHCHO=cloud_species_interact(6)
    ictmHNO2=cloud_species_interact(7)
    ictmO3=cloud_species_interact(8)
    ictmOH=cloud_species_interact(9)
    ictmHO2=cloud_species_interact(10)
    ictmNO3=cloud_species_interact(11)
    ictmNO=cloud_species_interact(12)
    ictmNO2=cloud_species_interact(13)
    ictmPAN=cloud_species_interact(14)
    ictmH2SO4=cloud_species_interact(15)

    !!initialize basic physical and chemical parameters
    allocate(saturation_pressure(nesp_aer))
    allocate(partition_coefficient(nesp_aer))
    allocate(vaporization_enthalpy(nesp_aer))
    allocate(accomodation_coefficient(nesp_aer))
    allocate(surface_tension(nesp_aer))
    allocate(saturation_pressure_mass(nesp_aer))
    allocate(saturation_pressure_torr(nesp_aer))
    allocate(deliquescence_relative_humidity(nesp_aer))
    allocate(molecular_weight_aer(nesp_aer))
    allocate(molecular_diameter(nesp_aer))
    allocate(collision_factor_aer(nesp_aer))
    allocate(mass_density_aer(nesp_aer))
    allocate(quadratic_speed(nesp_aer))
    allocate(diffusion_coef(nesp_aer))
    allocate(soa_sat_conc(nesp_aer))
    allocate(soa_part_coef(nesp_aer))

    do jesp = 1,nesp_aer
      saturation_pressure(jesp)=saturation_pressure_loc(jesp)
      partition_coefficient(jesp)=partition_coefficient_loc(jesp)
      vaporization_enthalpy(jesp)=vaporization_enthalpy_loc(jesp)
      surface_tension(jesp)=surface_tension_loc(jesp)
      accomodation_coefficient(jesp)=accomodation_coefficient_loc(jesp)
      saturation_pressure_mass(jesp)=saturation_pressure_mass_loc(jesp)
      saturation_pressure_torr(jesp)=saturation_pressure_torr_loc(jesp)
      deliquescence_relative_humidity(jesp)=deliquescence_relative_humidity_loc(jesp)
      molecular_weight_aer(jesp)=molecular_weight_aer_loc(jesp)
      molecular_diameter(jesp)=molecular_diameter_aer_loc(jesp)
      collision_factor_aer(jesp)=collision_factor_aer_loc(jesp)
      mass_density_aer(jesp)=mass_density_aer_loc(jesp)!µg/µm3 ?
      soa_part_coef(jesp)=0.d0
      soa_sat_conc(jesp)=0.d0
      diffusion_coef(jesp)=0.d0
      quadratic_speed(jesp)=0.d0
    enddo

!*       molecular_weight_inside(*)   molar weight of inorganic species *
!*       in aqueous_phase           µg.mol-1                *
    molecular_weight_inside(IH)=1.0D06
    molecular_weight_inside(INa)=23.0D06
    molecular_weight_inside(INH4)=18.0D06
    molecular_weight_inside(ICl)=35.5D06
    molecular_weight_inside(ISO4)=96.0D06
    molecular_weight_inside(IHSO4)=97.0D06
    molecular_weight_inside(INO3)=63.0D06
    molecular_weight_inside(IH2O)=18.0D06
    molecular_weight_inside(INH3)=17.0D06
    molecular_weight_inside(IHCl)=36.5D06
    molecular_weight_inside(IHNO3)=63.0D06
    molecular_weight_inside(IOH)=17.0D06
!      molar weight of solids
    molecular_weight_solid(SNaNO3)=85.0D06
    molecular_weight_solid(SNH4NO3)=80.0D06
    molecular_weight_solid(SNACl)=58.5D06
    molecular_weight_solid(SNH4Cl)=53.5D06
    molecular_weight_solid(SNa2SO4)=142.0D06
    molecular_weight_solid(SNH42S4)=132.0D06
    molecular_weight_solid(SNaHSO4)=120.0D06
    molecular_weight_solid(SNH4HS4)=115.0D06
    molecular_weight_solid(SLC)=247.0D06
!      DENSITIES of solids
    mass_density_solid(SNaNO3)=2.260D-06
    mass_density_solid(SNH4NO3)=1.725D-06
    mass_density_solid(SNACl)=2.165D-06
    mass_density_solid(SNH4Cl)=1.530D-06
    mass_density_solid(SNa2SO4)=2.700D-06
    mass_density_solid(SNH42S4)=1.770D-06
    mass_density_solid(SNaHSO4)=2.740D-06
    mass_density_solid(SNH4HS4)=1.780D-06
    mass_density_solid(SLC)=1.770D-06    

  END subroutine Init_global_parameters

  subroutine free_allocated_memory()
    
    integer ierr

    deallocate(total_mass)
    deallocate(cell_mass)
    deallocate(concentration_gas)
    deallocate(concentration_number)
    deallocate(concentration_number_tmp)
    deallocate(concentration_mass)
    deallocate(concentration_mass_tmp)
    deallocate(concentration_inti)
    deallocate(cell_mass_av)
    deallocate(cell_diam_av)
    deallocate(cell_log_av)
    deallocate(total_aero_mass)
    deallocate(mass_total_grid)
    deallocate(Index_groups)
    deallocate(wet_mass)
    deallocate(wet_diameter)
    deallocate(wet_volume)
    deallocate(size_diam_av)
    deallocate(size_mass_av)
    deallocate(size_log_av)
    deallocate(size_sect)
    deallocate(diam_bound)
    deallocate(mass_bound)
    if (allocated(log_bound)) then
       deallocate(log_bound, stat=ierr)
    endif
    if (allocated(density_aer_bin)) then
       deallocate(density_aer_bin, stat=ierr)
    endif
    if (allocated(density_aer_size)) then
       deallocate(density_aer_size, stat=ierr)
       if (ierr .ne. 0) then
          stop "Deallocation error"
       endif
    endif
    deallocate(bin_mass)
    deallocate(bin_number)
    deallocate(rho_wet_cell)
    deallocate(kernel_coagulation)
    deallocate(total_bin_mass)
    deallocate(ce_kernal_coef)
    deallocate(Kelvin_effect_ext)
    deallocate(frac_grid)
    deallocate(dqdt)
    deallocate(concentration_index_iv)
    deallocate(concentration_index)
    deallocate(isorropia_species)
    deallocate(aec_species)
    deallocate(pankow_species)
    deallocate(poa_species)
    deallocate(aerosol_species_interact)
    deallocate(List_species)
    deallocate(saturation_pressure)
    deallocate(partition_coefficient)
    deallocate(vaporization_enthalpy)
    deallocate(accomodation_coefficient)
    deallocate(surface_tension)
    deallocate(saturation_pressure_mass)
    deallocate(saturation_pressure_torr)
    deallocate(deliquescence_relative_humidity)
    deallocate(molecular_weight_aer)
    deallocate(molecular_diameter)
    deallocate(collision_factor_aer)
    deallocate(mass_density_aer)
    deallocate(quadratic_speed)
    deallocate(diffusion_coef)
    deallocate(soa_sat_conc)
    deallocate(soa_part_coef)
    deallocate(discretization_composition)
  END subroutine free_allocated_memory

end module aInitialization
