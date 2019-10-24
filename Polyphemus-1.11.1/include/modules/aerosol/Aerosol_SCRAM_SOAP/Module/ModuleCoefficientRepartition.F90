!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Hilel Dergaoui, Edouard Debry, Karine Sartelet
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
!!    This module contains methods read coefficient repartition data
!!
!!-----------------------------------------------------------------------
Module bCoefficientRepartition
  use aInitialization
  implicit none

  ! Parameter and variable definitions
  type :: ptr_to_real_array
     integer n
     double precision, dimension(:), pointer :: arr
  end type ptr_to_real_array

  type :: ptr_to_integer_array
     integer n
     integer, dimension(:), pointer :: arr
  end type ptr_to_integer_array

  ! Repartition coefficients.
  type(ptr_to_real_array), dimension(:), allocatable :: repartition_coefficient
  type(ptr_to_real_array), dimension(:), allocatable :: repartition_coefficient_nc
  ! Index of repartition coefficients.
  type(ptr_to_integer_array), dimension(:), allocatable :: index1_repartition_coefficient
  type(ptr_to_integer_array), dimension(:), allocatable :: index1_repartition_coefficient_nc
  type(ptr_to_integer_array), dimension(:), allocatable :: index2_repartition_coefficient
  type(ptr_to_integer_array), dimension(:), allocatable :: index2_repartition_coefficient_nc
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
!$omp threadprivate(repartition_coefficient,repartition_coefficient_nc)
!$omp threadprivate(index1_repartition_coefficient,index1_repartition_coefficient_nc)
!$omp threadprivate(index2_repartition_coefficient,index2_repartition_coefficient_nc)
#endif
contains

  subroutine ReadCoefficient(coef_size,Ncoefficient,index_first,index_second,coefficient)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine reads coefficient repartition data
!     based on original data passed for c++ routine
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     coef_size: total data size of original data
!     Ncoefficient: array of number of coefficient data of each bin
!     index_first: original data of first index of coagulation couple
!     index_second: original data of second index of coagulation couple
!     coefficient: original data of coefficients
!
!------------------------------------------------------------------------   
    implicit none 

    integer:: coef_size
    integer:: Ncoef
    integer:: Nc_id
    integer:: Ncfix
    integer::j, l, i, i1, i2
    INTEGER, intent(in) :: Ncoefficient(N_size)!SZ
    INTEGER, intent(in) :: index_first(coef_size)!SZ
    INTEGER, intent(in) :: index_second(coef_size)!SZ
    DOUBLE PRECISION, intent(in) :: coefficient(coef_size)!SZ
    DOUBLE PRECISION :: Tempcoef

    allocate(index1_repartition_coefficient(N_size))
    allocate(index1_repartition_coefficient_nc(N_size))!to read the oringinal version
    allocate(index2_repartition_coefficient(N_size))
    allocate(index2_repartition_coefficient_nc(N_size))
    allocate(repartition_coefficient(N_size))
    allocate(repartition_coefficient_nc(N_size))

    Nc_id=1!index of original data
    
    do i = 1, N_size
      Ncoef=Ncoefficient(i)
      index1_repartition_coefficient_nc(i)%n = Ncoef
      index2_repartition_coefficient_nc(i)%n = Ncoef
      repartition_coefficient_nc(i)%n= Ncoef
      allocate(index1_repartition_coefficient_nc(i)%arr(Ncoef))
      allocate(index2_repartition_coefficient_nc(i)%arr(Ncoef))
      allocate(repartition_coefficient_nc(i)%arr(Ncoef))
      
	!Problem: In NetCDF index of bin start from 0; However, in fortran index of bin start from 1
      do j= 1, Ncoef
	index1_repartition_coefficient_nc(i)%arr(j)=index_first(Nc_id)+1
	index2_repartition_coefficient_nc(i)%arr(j)=index_second(Nc_id)+1
	repartition_coefficient_nc(i)%arr(j)=coefficient(Nc_id)
	Nc_id=Nc_id+1!move to the next Nc data
      enddo

    end do
    
    do j = 1, N_size
      Ncfix=repartition_coefficient_nc(j)%n
      Ncoef=Ncfix
      do l=1,repartition_coefficient_nc(j)%n! ckeck symmetric case
	i1=index1_repartition_coefficient_nc(j)%arr(l)
	i2=index2_repartition_coefficient_nc(j)%arr(l)
	if(i1/=i2) then
	  Ncoef=Ncoef+1
	endif
      enddo
      repartition_coefficient(j)%n=Ncoef
      index1_repartition_coefficient(j)%n=Ncoef
      index2_repartition_coefficient(j)%n=Ncoef
      allocate(repartition_coefficient(j)%arr(Ncoef))
      allocate(index1_repartition_coefficient(j)%arr(Ncoef))
      allocate(index2_repartition_coefficient(j)%arr(Ncoef))
      do l=1,repartition_coefficient_nc(j)%n! Reassign the coefficient and add missing symmetric coefficients
	i1=index1_repartition_coefficient_nc(j)%arr(l)
	index1_repartition_coefficient(j)%arr(l)=i1
	i2=index2_repartition_coefficient_nc(j)%arr(l)
	index2_repartition_coefficient(j)%arr(l)=i2
	Tempcoef=repartition_coefficient_nc(j)%arr(l)
	repartition_coefficient(j)%arr(l)=Tempcoef
	if(i1/=i2) then
	  Ncfix=Ncfix+1
	  index1_repartition_coefficient(j)%arr(Ncfix)=i2!set symmetric coefficient
	  index2_repartition_coefficient(j)%arr(Ncfix)=i1
	  repartition_coefficient(j)%arr(Ncfix)=Tempcoef
	  if(Ncfix>Ncoef) then
	      print*,"Error: CoefficientRepartition Aarry out of bounds...."
	      STOP
	  endif
	endif
      enddo
    end do
    
  end subroutine ReadCoefficient

  subroutine check_repart_coeff()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine checks the quality of coefficient repartition data
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------
    implicit none
    integer::k,i,j,l
    double precision:: sumcr(N_size,N_size)
    
    sumcr=0.d0
    do k=1,N_size
      do l=1,repartition_coefficient(k)%n
	i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
	j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
	sumcr(i,j)=sumcr(i,j)+repartition_coefficient(k)%arr(l)
      enddo
    enddo
    
    do k=1,N_size 
      do l=1,repartition_coefficient(k)%n
	i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
	j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
	if(sumcr(i,j).ne.1.d0) then
	  if(repartition_coefficient(k)%arr(l)*sumcr(i,j).gt.0.d0) then
	    repartition_coefficient(k)%arr(l)=repartition_coefficient(k)%arr(l)/sumcr(i,j)
	  endif
	endif
      enddo
    enddo    
  end subroutine

  subroutine DeallocateCoefficientRepartition()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine Deallocate Coefficient Repartition arrays
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------  
    implicit none

    integer :: j
    !write (*,'(A, F8.3)') 'Deallocation routine'
    do j=1,N_size
       deallocate (repartition_coefficient(j)%arr)
       deallocate (index1_repartition_coefficient(j)%arr)
       deallocate (index2_repartition_coefficient(j)%arr)
       deallocate (repartition_coefficient_nc(j)%arr)
       deallocate (index1_repartition_coefficient_nc(j)%arr)       
       deallocate (index2_repartition_coefficient_nc(j)%arr)
    end do

    deallocate (repartition_coefficient)
    deallocate (index1_repartition_coefficient)
    deallocate (index2_repartition_coefficient)
    deallocate (repartition_coefficient_nc)
    deallocate (index1_repartition_coefficient_nc)
    deallocate (index2_repartition_coefficient_nc)
    
  end subroutine DeallocateCoefficientRepartition
  
  

end module bCoefficientRepartition
