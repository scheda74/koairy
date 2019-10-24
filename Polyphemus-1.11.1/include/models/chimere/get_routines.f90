subroutine get_parameters(Nspresc_, Ntabmax_, Nlevphotmax_, Nphotmax_, Ntabuzenmax_, &
     Nspec_, Nreac_, Nkc_, Ms_, Nzonal_, Nmerid_, &
     Nverti_, Nemisa_, Nemisb_, Nlevemis_, Nfam_)

  use chimere_params

  implicit none

  integer, intent(out) :: Nspresc_, Ntabmax_, Nlevphotmax_, Nphotmax_, Ntabuzenmax_, &
       Nspec_, Nreac_, Nkc_, Ms_, Nzonal_, Nmerid_, &
       Nverti_, Nemisa_, Nemisb_, Nlevemis_, Nfam_

  Nspresc_ = nspresc
  Ntabmax_ = ntabmax
  Nlevphotmax_ = nlevphotmax
  Nphotmax_ = nphotmax
  Ntabuzenmax_ = ntabuzenmax
  Nspec_ = nspec
  Nreac_ = nreac
  Nkc_ = nkc
  Ms_ = ms
  Nzonal_ = nzonal
  Nmerid_ = nmerid
  Nverti_ = nverti
  Nemisa_ = nemisa
  Nemisb_ = nemisb
  Nlevemis_ = nlevemis
  Nfam_ = nfam

end subroutine get_parameters

subroutine get_species(species_list_tmp)

  use chimere_params
  use chimere_common

  implicit none

  integer :: i, i1, i2
  character(len=3000), intent(inout) :: species_list_tmp

  i1 = 1
  i2 = 1
  do i = 1,nspectot
     i2 = i1 + len_trim(species(i)%name)
     species_list_tmp(i1:i2) = trim(species(i)%name)
     i1 = i2 + 1
  enddo

  species_list_tmp = trim(species_list_tmp)//CHAR(0)

end subroutine get_species
