subroutine forward_integrun

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
  integer :: ierr
  integer :: iprint

  !****************************************************************************

  iprint = ihourrun

  thour = 0.5d0/nphour
  call master_send_hourly
  call master_send_conc
  call mpi_barrier(mpi_comm_world,ierr)

  do np=1,nphour
     call master_locvalues
     call mpi_barrier(mpi_comm_world,ierr) !modif
     call master_send_frac_hourly
     call master_send_conc_bounds
     call mpi_barrier(mpi_comm_world,ierr)
     !lmbb: old loop for sequential version
     !do ir=1,ichemstep
     ! twostep running
     ! end twostep
     !enddo

     !  Timing update
     thour = thour + dun/nphour
     djul = djul + dun/(nphour*nhourpday)
     if(np.eq.nphour) ihourrun = ihourrun + 1

     ! Printout. I need locvalues and conc
     if(np.eq.nphour) then
        call master_recv_locvalues
        call master_recv_conc
        call master_recv_excha
        iprint = iprint+1
        call outprint(iprint,.false.)
     end if
  enddo ! np=1,nphour

  call renewhour

  !  Save of all concentrations and accumulated variables

  if(mod(ihourrun,nsaveconcs).eq.0) then
     call write_concs
     ! lmbb: flag for pops
     if(pops.eq.1) call write_excha
  endif
  if(mod(ihourrun,nsavedepos).eq.0) then
     call write_depo
  endif

end subroutine forward_integrun
