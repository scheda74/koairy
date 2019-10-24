subroutine initstep

  use chimere_common

  implicit none

!  Run integration
  if(rank==0) then
     call initstep_integrun
  else
     call initstep_worker
  endif

end subroutine initstep
