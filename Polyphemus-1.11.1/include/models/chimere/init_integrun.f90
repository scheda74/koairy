subroutine init_integrun

  !  Integration of the model

  use chimere_consts
  use chimere_params
  use chimere_common
  use message_defs
  use master_message_subs

  implicit none

  !include 'mpif.h'


  !****************************************************************************
  integer :: np,ir
  integer :: ksens
  integer :: ierr
  integer :: iprint

  !****************************************************************************

  ndoms = nzdoms*nmdoms

  if(ndoms==0) then
     print *
     print *, ' *** Error : NZDOMS or NMDOMS should not be null ! ***'
     print *, ' *** Check your configuration script. Exiting.     ***'
     call mpi_finalize(ierr)
     stop
  end if

  iprint = 0
  ihourrun = 0
  thour = 0.5d0/nphour
  ksens = 2

  call master_init_mpi
  call master_send_once
  call mpi_barrier(mpi_comm_world,ierr)

  ! receive initial locvalues from workers for printout
  call master_locvalues
  call master_send_conc
  call master_send_excha
  call mpi_barrier(mpi_comm_world,ierr)
  call master_recv_locvalues

  ! printout of the initial concentrations and meteo parameters
  call outprint(iprint,.true.)

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine init_integrun
