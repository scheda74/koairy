C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Pierre Plion and Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
C     component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

      SUBROUTINE ROS2CONC(neq,nesp_aer,nbin_aer,q1,q,iq,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,bin_density)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutines solves the GDE with the scond-order Rosenbrock
C     algorithm (ROS2).
C
C     The canonic form of ROS2 is:
C     (I-Gamma*Jn*h) * k1 = fn*h
C     (I-Gamma*J1*h) * k2 = f1*h - Gamma*(Jn+J1)*h*k1
C     c(n+1) = c(n) + (k1+k2)/2
C
C     where c is the variable, f the time derivative and J an
C     approximation of the Jacobian matrix.
C
C     To spare jacobian matrix computation
C     the same approximation can be used at both states
C     (only one LU factorisation):
C     (I-Gamma*Jn*h) * k1 = fn*h
C     (I-Gamma*Jn*h) * k2 = f1*h - 2*Gamma*Jn*h*k1
C     c(n+1) = c(n) + (k1+k2)/2
C
C     A cheaper implementation is furthermore possible
C     (avoiding a matrix-vector product)
C     (I-Gamma*Jn*h) * k1 = fn*h
C     (I-Gamma*Jn*h) * k2 = f1*h - 2*k1
C     c(n+1) = c(n) + (3* k1+k2)/2
C
C     The use of a poor approximation is another way to spare
C     computation. Here only the diagonal part is estimated
C     WITH PHYSICAL ASSUMPTIONS attempting to ensure positivity
C     if f < 0 , it is assumed zero-valued in zero (decreasing may stop)
C     f > 0 , it is assumed zero-valued for total available amount
C     of the considered specie (condensation may stop)
C     Both of these assumptions implies negative values for
C     jacobian diagonal part.
C     With these approximations, the use of canonic form is no longer
C     computational expensive.
C
C     ! Notations.
C     in Verwer (1)   , k is variation rate of variable c
C     in this software, k is first derivative, and then variation
C     (step) for the variable q
C
C     Reference:
C     J.G. Verwer, E.J. Spee, J.G. Blom, W.H. Hundsdorfer
C     "A second order Rosenbrock method applied to photochemical
C     dispersion problems" Report MAS-R9717 August 1997
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     NBIN_AER: number of aerosol bins.
C     IQ: aerosol pointers.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q : second-order evaluation of concentrations ([\mu.g.m^-3]).
C
C
C     -- OUTPUT VARIABLES
C
C     Q1 : first-order evaluation of concentrations ([\mu.g.m^-3]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C     2013/11/27: Added bin_density (Stephanie Deschamps, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Pierre Plion and Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'time.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varq.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q1(neq),q(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)

      DOUBLE PRECISION AA(nesp_aer,nbin_aer)

      INTEGER j
      DOUBLE PRECISION k1(neq),Jdn(neq)
      DOUBLE PRECISION k2(neq),Jd1(neq)
      DOUBLE PRECISION temp1,temp2
      INTEGER js,jesp

      DOUBLE PRECISION bin_density(nbin_aer)


C     parameter for partial implicitation in Rosenbrock's algorithm
C     1+sqrt(2)/2     = 1.7071   R. Djouad / J. Verwer 1st determination
C     1-sqrt(2)/2     = 0.2929               J. Verwer 2nd determination
C     (1+sqrt(1/3))/2 = 0.7887   3rd order for autonomous & linear
C     (1-sqrt(1/3))/2 = 0.2113   idem, 2nd determination
C     1/4            = 0.2500   A-stable & monotonic
C     1/2            = 0.5000   predicted step is still second order
C
C     every gamma values ensures 2nd order scheme
C
C     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
C     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
C     gamma =  1+-sqrt(1/2) ->  L-stability

      DOUBLE PRECISION Gamma
      PARAMETER ( Gamma= 1.7071D0)

C     ******1. compute the derivative k1 in (q ; Tin2)
C     ****** and the approximation of jacobian's diagional.

      CALL FGDE(neq,nesp_aer,nbin_aer,q,iq,k1,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA,
     s     bin_density)

C     Every dynamical variable protected against vanishing

      DO j = 1 , neq
C     if (k1(j).ne.k1(j)) k1(j) = 0.d0
         Jdn(j) = 0.d0
         if ( q(j).gt.1.d-26 .and. k1(j).lt.0.d0 )
     &        Jdn(j) = k1(j) / q(j)
      END DO

C     Condensation limited by available amount (temp1) in both gas phase
C     and small bins (assumed in equilibrium)
C     sum of condensation rates (temp2) is distributed between
C     dynamical bins

      DO jesp = E1, E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
         temp1 = q(IG(jesp))
         DO js = 1, ICUT2
            temp1 = temp1 + q(IQ(jesp,js))
         END DO

         temp2 = 0.d0
         DO js = (ICUT2+1), nbin_aer
            j  = IQ(jesp,js)
            IF (k1(j).GT.0.d0) temp2 = temp2 + k1(j)
         END DO

         DO js = (ICUT2+1), nbin_aer
            j  = IQ(jesp,js)
            IF ( (temp1-q(j)).GT.1.d-26 .AND. k1(j).GT.0.D0)
     &           Jdn(j) = -temp2 / (temp1-q(j))
         END DO
         ENDIF
      END DO

C     ******2. first evaluation

      DO j = 1 , neq
         k1(j) = k1(j) * DT2 / ( 1.d0 - Gamma * Jdn(j) * DT2 )
      END DO

      CALL KLIMIT(neq,nesp_aer,nbin_aer,q,iq,k1,AA)
      DO j = 1 , neq
         q1(j) = DMAX1 ( 0.01d0*q(j) , (q(j)+k1(j)) )
      ENDDO

C     ******3. compute the derivative k2 in (q+k1 ; Tin2+DT2)
C     ******&  the approximation of jacobian's diagonal

      TIN2 = TIN2 + DT2
      CALL FGDE(neq,nesp_aer,nbin_aer,q1,iq,k2,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA, 
     s     bin_density)

      DO j = 1 , neq
         IF (k2(j).ne.k2(j)) k2(j) = 0.d0
         Jd1(j) = 0.d0
         IF ( q1(j).gt.1.d-26 .and. k2(j).lt.0.d0 )
     &        Jd1(j) = k2(j) / q1(j)
      END DO

      DO jesp = E1, E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
         temp1 = q1(IG(jesp))
         DO js = 1, ICUT2
            temp1 = temp1 + q1(IQ(jesp,js))
         END DO

         temp2 = 0.d0
         DO js = (ICUT2+1), nbin_aer
            j = IQ(jesp,js)
            if(k2(j).gt.0.d0) temp2 = temp2 + k2(j)
         END DO
         DO js = (ICUT2+1), nbin_aer
            j = IQ(jesp,js)
            IF ( (temp1-q1(j)).gt.1.d-26 .and. k2(j).gt.0.D0)
     &           Jd1(j) = -temp2  / (temp1-q1(j))
         END DO
         ENDIF
      END DO

C     ******4. second evaluation

      DO j = 1 , neq
         k2(j) = (k2(j)*DT2 - Gamma*DT2*(Jdn(j)+Jd1(j))*k1(j))
     &        / (1.d0      - Gamma*        Jd1(j) *DT2  )
      END DO

      CALL KLIMIT(neq,nesp_aer,nbin_aer,q,iq,k2,AA)
      DO j = 1 , neq
         q(j) = DMAX1 ( 0.01d0*q(j) , (q(j)+0.5d0*(k1(j)+k2(j))) )
      ENDDO

      END
