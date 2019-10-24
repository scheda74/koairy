module chimere_perturbation

  use worker_params
  use worker_common

  implicit none

  real,dimension(nspec, 1) :: perturbation_depoloc
  real,dimension(nreac, 1, 1) :: perturbation_rate

contains

  subroutine apply_perturbation
    implicit none

    integer :: nd,nr,izo,ime,ivert

    do ime=1,nmerid
       do izo=1,nzonal
          do nd=1,ndepo
             depoloc(nd,izo,ime) = depoloc(nd,izo,ime) * dble(perturbation_depoloc(nd, 1))
          enddo
       enddo
    enddo

    do nr=1,nreac
       do ivert=1,nverti
          do ime=1,nmerid
             do izo=1,nzonal
                rate(nr,izo,ime,ivert) = rate(nr,izo,ime,ivert) * dble(perturbation_rate(nr, 1, 1))
             enddo
          enddo
       enddo
    enddo
  end subroutine apply_perturbation


end module chimere_perturbation

