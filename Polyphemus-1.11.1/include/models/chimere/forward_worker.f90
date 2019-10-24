subroutine forward_worker

  use chimere_consts
  use worker_params
  use worker_common
  use message_defs
  use worker_message_subs
  use twostep_mod

  use chimere_perturbation

  implicit none

  integer :: ir,np
  integer :: ierr

  call worker_recv_hourly
  call worker_recv_conc
  thour = 0.5d0/nphour
  call mpi_barrier(mpi_comm_world,ierr)

  do np=1,nphour

     call mpi_barrier(mpi_comm_world,ierr) !modif
     call worker_recv_frac_hourly
     call worker_recv_conc_bounds

     call zenith
     call locvalues
     call physics2

     call apply_perturbation

     call mpi_barrier(mpi_comm_world,ierr)

     do ir=1,ichemstep
        call twostep
        call worker_update_halo ! between workers. Master is not involved.
     enddo

     !  Timing update
     thour = thour + dun/nphour
     djul = djul + dun/(nphour*nhourpday)
     if(np.eq.nphour) ihourrun = ihourrun + 1

     ! Printout. Master needs locvalues and conc
     if(np.eq.nphour) then
        call worker_send_locvalues
        call worker_send_conc
        call worker_send_excha

     endif

  enddo ! np=1,nphour


end subroutine forward_worker
