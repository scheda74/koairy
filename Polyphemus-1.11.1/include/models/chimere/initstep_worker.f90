subroutine initstep_worker

  use chimere_consts
  use worker_params
  use worker_common
  use message_defs
  use worker_message_subs
  use twostep_mod

  implicit none

  integer :: ierr

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine initstep_worker
