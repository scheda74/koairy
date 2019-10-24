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
!!    This module contains methods for for mass and number redistribution
!!    for fixed euler scheme in 3D application.
!!-----------------------------------------------------------------------
Module eRedistribution
  use aInitialization
  use dPhysicalbalance
  use cThermodynamics
  implicit none

contains
  subroutine redistribution_size(scheme)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine provides entries for different size redistribution methods
!     the chosen numerical solver.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     scheme: type of chosen redistribution method
!
!------------------------------------------------------------------------
    implicit none
    integer::k,j,f,jesp,s
    integer:: scheme !redistribution scheme 3 = euler_mass 4 = euler_number 5 = hemen 6 = euler_coupled
    integer:: section_pass!bin include 100nm
    double precision:: Qesp(N_sizebin, N_aerosol)!Temperature mass concentration on each fraction
    double precision:: N(N_sizebin)!Temperature mass concentration on each fraction
    double precision:: totQ(N_sizebin)
    double precision:: d(N_sizebin)
    double precision:: tmp_n
    Qesp=0.d0
    tmp_n=0.d0

    do k=1,N_sizebin
      if(diam_bound(k).lt.1.d-2.and.diam_bound(k+1).ge.1.d-2) then
      section_pass=k
      endif
    enddo
    !see the distribution in one fixed fraction as the same case of internal mixing
    do f=1,N_fracmax
      !extrac the mass and number distribution of each fraction out
      do k=1,N_sizebin
	j=concentration_index_iv(k,f)
	N(k)=concentration_number(j)
	totQ(k)=cell_mass(j)
! 	if(IsNaN(N(k)*0.d0)) then
! 	  print*,"Error of infinity/NaN before rdb",N(k)
! 	  stop
! 	endif	
	do s=1,N_aerosol
	  jesp=List_species(s)
! 	  if(IsNaN(Qesp(k,jesp)*0.d0)) then
! 	    print*,"Error of infinity/NaN before rdb",Qesp(k,jesp)
! 	    stop
! 	  endif	  
	  Qesp(k,jesp)=concentration_mass(j,jesp)
	  d(k)=size_diam_av(k)
	enddo
      enddo
      call redistribution(N_sizebin,N_aerosol,EH2O,diam_bound, d, scheme, &
      section_pass, mass_density_aer, DQLIMIT, Qesp, N, totQ,&
      with_fixed_density)
      !update distribution of each bin in current fraction section size_diam_av
      do k=1,N_sizebin
	j=concentration_index_iv(k,f)
	concentration_number(j)=N(k)
	cell_mass(j)=totQ(k)
! 	if(IsNaN(N(k)*0.d0)) then
! 	  print*,"Error of infinity/NaN after rdb",N(k)
! 	  stop
! 	endif
	do s=1,N_aerosol
	  jesp=List_species(s)
	  concentration_mass(j,jesp)=Qesp(k,jesp)
! 	  if(IsNaN(Qesp(k,jesp)*0.d0)) then
! 	    print*,"Error of infinity/NaN after rdb",Qesp(k,jesp)
! 	    stop
! 	  endif
	enddo
      enddo
    enddo
    !update wet and dry diameter
!     call compute_wet_mass_diameter(1,N_size,concentration_mass,concentration_number,&
! 	  concentration_inti,wet_mass,wet_diameter,wet_volume,cell_diam_av)
	  
  end subroutine redistribution_size

  subroutine redistribution_fraction()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute mass and number
!     based on the new fraction composation of aerosol
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------
    implicit none
    integer::k,i,i1,g,dj,j,jesp,s
    double precision:: mass_groups(N_size,N_groups)!mass by groups
    double precision:: mass_var(N_size,(N_aerosol-1))!mass variation map
    double precision:: numb_var(N_size)!number variation map
    double precision::mass_total,f1,f2
    !calculate redistribution map
!! Initialize to zero the variation in mass and number concentrations due to 
!! redistribution
    do j= 1, N_size  
       numb_var(j)=0.d0
       do jesp=1,(N_aerosol-1)
       	  mass_var(j,jesp)=0.d0
       enddo
       do g= 1, N_groups
	  frac_grid(j,g)=0.d0!Fraction of group g in grid j
	  mass_groups(j,g)=0.d0
       enddo
    enddo

    do j= 1, N_size !!for each grid the number of size sections x fraction sections
      mass_total_grid (j)=0.d0
!      Compute fraction of each species in the grid point before redistribution      
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	g=index_groups(jesp)
	mass_groups(j,g)=mass_groups(j,g)+concentration_mass(j,jesp)
	mass_total_grid (j)=mass_total_grid (j) + concentration_mass(j,jesp)
      enddo

      do g= 1, N_groups
	if(mass_total_grid (j).gt.0d0) then
	  frac_grid(j,g) =mass_groups(j,g)/mass_total_grid (j)
	endif
      enddo
      
      i1 =concentration_index(j, 1)!size bin index
!     Loop on fraction combinations in the size bin i1 to find out where the 
!     original fraction section (before condensation/evaporation) has moved to.
      do i=1,N_fracmax
	dj=0
!       Check for each species whether this fraction bin is the correct one
	do g=1,N_groups-1
          f1=discretization_composition( i, g, 1)
	  f2=discretization_composition( i, g, 2)
          if(f1.eq.0.d0) then
	    if(frac_grid(j,g).ge.f1.and.frac_grid(j,g).le.f2) then
	      dj=dj+1
	    endif
	  else
	    if(frac_grid(j,g).gt.f1.and.frac_grid(j,g).le.f2) then
	      dj=dj+1!dj is the number of matched groups
	    endif
	  endif
	enddo
	
!       The fraction of each species is identified. 
!       i is then the correct fraction combination (after redistribution). 
	if(dj.eq.(N_groups-1)) then
	  k=concentration_index_iv(i1,i)
	  numb_var(k)=numb_var(k)+concentration_number(j)
	  concentration_number(j)=0.d0!minimize numerical disarrange due to difference scale
	  do s=1,(N_aerosol-1)
	    jesp=List_species(s)
	    mass_var(k,jesp)=mass_var(k,jesp)+concentration_mass(j,jesp)
	    concentration_mass(j,jesp)=0.d0
	  enddo
	endif
      enddo
    enddo
    
    !Update mass and number concentrations after redistribution
    do j= 1, N_size !!for each grid point
      concentration_number(j)=concentration_number(j)+numb_var(j)
      mass_total_grid (j)=0.d0
      do g= 1, N_groups
	frac_grid(j,g)=0.d0!Fraction of group g in grid j
	mass_groups(j,g)=0.d0
      enddo      
      !renew Temperature frac of each grid
      do s= 1, (N_aerosol-1)
	jesp=List_species(s)
	g=index_groups(jesp)
	concentration_mass(j,jesp)=concentration_mass(j,jesp)+mass_var(j,jesp)
	mass_total_grid (j) = mass_total_grid (j) + concentration_mass(j,jesp)
	mass_groups(j,g)=mass_groups(j,g)+concentration_mass(j,jesp)
      enddo
      do g= 1, N_groups
	if(mass_total_grid (j).gt.0d0) then
	  frac_grid(j,g) =mass_groups(j,g)/mass_total_grid (j)
	  !print*,"j",j,"g",g,"frac_grid(j,g)",frac_grid(j,g)	  
! 	  if(frac_grid(j,g).gt.1.d0) then
! 	    print*,"j",j,"g",g,"frac_grid(j,g)",frac_grid(j,g),&
! 	    mass_groups(j,g),mass_total_grid (j)
! 	  endif
	endif
      enddo
    enddo
    
    
  end subroutine redistribution_fraction
  
end Module eRedistribution
