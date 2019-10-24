subroutine init

!  Main program for simulation with CHIMERE

  use chimere_common
  use master_message_subs
  implicit none
  integer :: ierr

  ! initial time-step values
  dtr=3600d0/(nphour_ref*ichemstep)
  dtr2=2d0*dtr/3d0
  nphour=nphour_ref

  !  All readings and initializations

  call init_mpi

  if(rank==0) call inichimere

!  Run integration
  if(rank==0) then
     call init_integrun
  else
     call init_worker
  endif

end subroutine init
