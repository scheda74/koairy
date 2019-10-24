SUBROUTINE REDISTRIBUTION(ns, naer, EH2O, dbound, fixed_diameter, &
    fixed_density_aer, IDENS, &
    scheme, section_pass, aerosol_density, DQLIMIT, Qesp, N, Q) 

!!$------------------------------------------------------------------------
!!$     
!!$     -- INPUT VARIABLES
!!$     
!!$     
!!$     ns             : number of sections
!!$     naer           : number of species
!!$     dbound         : list of limit bound diameter [\mu m]
!!$     fixed_diameter : list of mean diameter 
!!$     scheme         : redistribution scheme
!!$                      3 = euler_mass
!!$                      4 = euler_number
!!$                      5 = hemen 
!!$                      6 = euler_coupled
!!$     section_pass   : bin include 100nm  
!!$     aerosol_density : list of liquid mass density of each species
!!$ 
!!$     -- VARIABLES
!!$     
!!$     rho            : density per bin
!!$     
!!$     Eps_diam       : tolerance to the diameter which may be slightly above 
!!$                      the edge of the section
!!$     Eps_dbl_prec   : tolerance lower in the case where the diameter so little 
!!$                      increases or decreases as one considers that it is 
!!$                      not happening 
!!$                    : it can not work with a time not too restrictive 
!!$     grand          : list of 0 or 1
!!$                      1 = cutting with the upper box
!!$                      0 = cutting with the lower box
!!$     d_after_cond   : list of mean diameter after condensation/evaporation
!!$     kloc           : list of bin where is diam
!!$     logdiam        : log(d)
!!$     X              : log(diam)
!!$     alpha          : list of fraction of each species in Q
!!$
!!$     -- INPUT/OUTPUT VARIABLES
!!$      
!!$     N            : Number concentration per bin
!!$     Qesp         : Mass concentration per bin and species
!!$     Q            : Total mass concentration per bin
!!$      
!!$     -- OUTPUT VARIABLES
!!$     
!!$     
!!$------------------------------------------------------------------------

  IMPLICIT NONE
  INCLUDE '../INC/parameuler.inc'

  ! ------ Input
  INTEGER, INTENT(in) :: ns, naer, section_pass, scheme
  DOUBLE PRECISION, DIMENSION(ns+1), INTENT(in) :: dbound
  DOUBLE PRECISION, DIMENSION(ns) :: fixed_diameter
  DOUBLE PRECISION, DIMENSION(naer), INTENT(in) :: aerosol_density
  INTEGER, INTENT(in) :: EH2O
  DOUBLE PRECISION, INTENT(in) :: fixed_density_aer
  INTEGER, INTENT(in) :: idens

  ! ------ Input/Output
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N, Q
  DOUBLE PRECISION, DIMENSION(ns, naer), INTENT(inout) :: Qesp

  ! ------ 
  INTEGER k, js, jaer, loc
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  INTEGER, DIMENSION(ns) :: grand
  INTEGER, DIMENSION(ns) :: kloc
  DOUBLE PRECISION :: dbis, dnew
  DOUBLE PRECISION , DIMENSION(ns) :: X, d_after_cond, logfixed_diameter
  DOUBLE PRECISION, DIMENSION(ns,naer) :: alpha 
  
  DOUBLE PRECISION :: frac, testalpha
  
  !! Calcul of d(k) = dsqrt(dbound(k)*dbound(k+1))
  !! Reestimate d(k) from number and mass concentration BEFORE condensation/evaporation
!  DO k = 1,ns
!    fixed_diameter(k) = dsqrt(dbound(k) * dbound(k+1))
!  ENDDO
  ! DO js = 1, ns
  !    DO jaer = 1, naer
  !       if (Qesp(js, jaer) .LE. TINYM) then
  !          Qesp(js, jaer) = 0.d0
  !       endif
  !    ENDDO
  ! ENDDO
  

  !***** Calcul dry total mass per bin
  Q=0.D0
  DO js = 1, ns
     rho(js) = fixed_density_aer 
     DO jaer = 1, naer
        if (jaer .NE. EH2O) then
           Q(js)= Q(js) + Qesp(js, jaer) 
        endif
     ENDDO
  ENDDO


  !***** Calcul fraction of each composant of mass
  DO js = 1, ns
     testalpha = 0.D0
     DO jaer = 1, naer
        IF (Q(js) .NE. 0.D0) THEN
           alpha(js,jaer) = Qesp(js, jaer)/Q(js)
           testalpha = testalpha + alpha(js,jaer) 
        ELSE
           alpha(js, jaer) = 0.D0
        ENDIF
     ENDDO
!     if (Qesp(js,EH2O).GT.0.d0) then
!        print *,js, "Tot ALPHA",testalpha, "Qh2o", Qesp(js,EH2O), "Qtot", Q(js)
!     endif
  ENDDO
  
  !****** Calcul of new mean diameter after c/e
  DO k = 1, ns
     if (isNaN(N(k))) then
        print *, k, "N(k) NaN dans redist"
     endif
     if (idens .eq. 1) CALL compute_density(ns,naer, EH2O,TINYM,Qesp,aerosol_density,k,rho(k))
     IF ( N(k) .GT. TINYN .and. Q(k) .GT. TINYM) THEN
        dbis = ((Q(k) * 6D0)/(PI * rho(k) * N(k)))**(1D0/3D0)
        IF (dbis .LT. dbound(1)) THEN
           Q(k) = 0.D0
           do jaer = 1,naer
              Qesp(k,jaer)  = 0.d0
           enddo
           N(k) = 0.D0
           dbis = fixed_diameter(1)
           kloc(k) = 1 
        ELSEIF (dbis .GT. dbound(ns+1)) THEN
           dbis = dbound(ns+1)
           kloc(k)= ns 
        ELSE
           loc = 1
           DO WHILE(dbis .GT. dbound(loc+1))
              loc=loc+1
           ENDDO
           kloc(k) = loc
        ENDIF
     ELSE
        dbis = fixed_diameter(k)
        kloc(k)=k
     ENDIF

     d_after_cond(k) = dbis
     if (isnan(d_after_cond(k))) then
        print*, k, "d_after_cond nan dans redistribution", &
             "  Q ",Q(k),"  rho ",rho(k), "  N ", N(k)
     endif
     X(k) = DLOG10(dbis)
     logfixed_diameter(k) = DLOG10(fixed_diameter(kloc(k)))
     IF(dbis .GT. fixed_diameter(kloc(k))) THEN
        grand(k) = 1
     ELSE
        grand(k) = 0 
     ENDIF

  ENDDO

  !****** Select the redistribution method

  SELECT CASE (scheme)
  CASE (3)
     CALL EULER_MASS(ns, naer, eh2o,dbound, grand, &
          X, d_after_cond, fixed_diameter, logfixed_diameter, &
          fixed_density_aer, idens, kloc, aerosol_density, Qesp, N)

  CASE (4)
    CALL EULER_NUMBER(ns, naer,eh2o, dbound, grand, alpha, &
          fixed_diameter, d_after_cond, X, logfixed_diameter, &
          fixed_density_aer, idens, kloc, aerosol_density, rho, Qesp, N)
 
  CASE (5)
    CALL HEMEN(ns, naer, eh2o,grand, section_pass, fixed_diameter, &
          dbound, X, d_after_cond, logfixed_diameter, &
          fixed_density_aer, idens, kloc, alpha, aerosol_density, &
          rho, Qesp, N)

  CASE (6)
    CALL EULER_COUPLED(ns, naer,eh2o,dbound, grand, fixed_diameter, d_after_cond, &
           fixed_density_aer, idens, kloc, alpha, aerosol_density, Qesp, N)

  CASE (7)
    CALL EULER_NUMBER_NORM(ns, naer, eh2o,dbound, grand, alpha, &
          fixed_diameter, d_after_cond, X, logfixed_diameter, &
          fixed_density_aer, idens, kloc, aerosol_density,DQLIMIT, rho, Qesp, N)
 
  CASE (8)
    CALL HEMEN_NORM(ns, naer,  EH2O, grand, section_pass, &
          dbound, X, d_after_cond, fixed_diameter, fixed_density_aer, idens, & 
          kloc, alpha, aerosol_density, DQLIMIT, rho, Qesp, N)

  CASE (9)
    CALL EULER_COUPLED_NORM(ns, naer, eh2o,dbound, grand, fixed_diameter, d_after_cond, &
           fixed_density_aer, idens, kloc, alpha, aerosol_density, DQLIMIT, Qesp, N)

  CASE (10)
    CALL MOVING_DIAM(ns, naer, eh2o, kloc, Qesp, N)  
  
  CASE DEFAULT
     PRINT*, "Please choose from the following redistribution methods : ", &
          "number-conserving, interpolation, euler-mass, euler-number, ",&
          "hemen, euler-coupled."
  END SELECT

! Total dry mass
  Q=0.D0
  DO js = 1, ns
     DO jaer = 1, naer
        if (jaer .NE. EH2O) then
           Q(js)= Q(js) + Qesp(js, jaer) 
        endif
     ENDDO
  ENDDO

  ! Set concentrations to zero if the diameter is too big.
  DO k = 1, ns
     if (idens .eq. 1) CALL compute_density(ns,naer, EH2O,TINYM,Qesp,aerosol_density,k,rho(k))
     IF (N(k) .GT. TINYN .and. Q(k) .GT. TINYM) THEN
        dnew = ((Q(k) * 6D0)/(PI * rho(k) * N(k)))**(1D0/3D0)
        IF (dnew .GT. 1.5D0 * dbound(Ns+1)) THEN
           N(k) = 0.d0
           Q(k) = 0.d0
           DO jaer = 1, naer
              Qesp(k, jaer) = 0.d0
           ENDDO
        ENDIF
     ENDIF
  ENDDO

END SUBROUTINE REDISTRIBUTION

