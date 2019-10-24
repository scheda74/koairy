subroutine initstep_integrun

  !  Integration of the model

  use chimere_consts
  use chimere_params
  use chimere_common
  use message_defs
  use master_message_subs

  implicit none

  !include 'mpif.h'


  !****************************************************************************
  integer :: ksens
  integer :: ierr

  !****************************************************************************

  ndoms = nzdoms*nmdoms

  ksens = 2

  call mpi_barrier(mpi_comm_world,ierr)
  ksens = 2
  call readhour(ksens)
  call checkcfl

end subroutine initstep_integrun
