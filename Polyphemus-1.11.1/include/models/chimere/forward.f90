subroutine forward

  use chimere_common
  use chimere_aerosol_anthropic_emission

  implicit none

!  Run integration
  if(rank==0) then
     call compute_aerosol_anthropic_emission
     call forward_integrun
  else
     call forward_worker
  endif

end subroutine forward
