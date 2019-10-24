subroutine cleanup

  use chimere_common
  use master_message_subs
  implicit none

!  Final Backup of all concentrations and clean-up

  if (rank==0) then
     call endchimere
     call master_mpi_finalize
  else
     call worker_mpi_finalize
  endif

end subroutine cleanup
