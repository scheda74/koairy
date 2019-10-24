subroutine init_worker

  use chimere_consts
  use worker_params
  use worker_common
  use message_defs
  use worker_message_subs
  use twostep_mod

  implicit none

  integer :: ir,np
  integer :: ierr

  ! initial time-step values
  dtr=3600d0/(nphour_ref*ichemstep)
  dtr2=2d0*dtr/3d0
  nphour=nphour_ref

  !call worker_init_mpi
  call worker_recv_once
  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_comm_rank(wrk_comm,rank,ierr)
  ndoms = nzdoms*nmdoms
  ihourrun = 0
  thour = 0.5d0/nphour

  call worker_recv_conc
  call worker_recv_excha
  call mpi_barrier(mpi_comm_world,ierr)

  !call worker_recv_airmloc
  call zenith
  call locvalues
  call physics2

  ! send locvalues to master for printing of initial state
  call worker_send_locvalues

  ! outprint


  call mpi_barrier(mpi_comm_world,ierr)


end subroutine init_worker
