module chimere_aerosol_anthropic_emission

  use chimere_consts
  use chimere_params
  use chimere_common

  implicit none

contains
  ! function computing aerosol anthropic emissions as in readhour.f90
  ! necessary for perturbation to be taken into account.

  subroutine compute_aerosol_anthropic_emission
    implicit none

    integer :: i,ne,nl,ispec,ik,ip
    integer :: izo,ime

    if (aero.eq.1.or.reactive.eq.1.or.carb.eq.1.or.trc.ge.1) then
       !  Particulate emissions
       if(bins.eq.1)then
          emipaloc = dzero
          do ne=1,nemisa
             ip = iprofemipa(ne)
             if(ip.ne.0) then
                ik = ispecemip(ip)
                do i=1,ms
                   do ime=1,nmerid
                      do izo=1,nzonal
                         do nl=1,nlevemis
                            emipaloc(ik,i,izo,ime,nl) = &
                                 emipaloc(ik,i,izo,ime,nl) + emisaloc(ne,izo,ime,nl)*emiprof(ip,i)
                         enddo
                      enddo
                   enddo
                enddo
             endif
          enddo
       endif
    endif

  end subroutine compute_aerosol_anthropic_emission

end module chimere_aerosol_anthropic_emission

