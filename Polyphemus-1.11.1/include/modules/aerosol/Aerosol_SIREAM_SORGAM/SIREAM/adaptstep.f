C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry and Bruno Sportisse
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

      SUBROUTINE ADAPTSTEP(neq,nbin_aer,q1,q,iq)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes an adaptive timestep by estimating a local
C     error on the basis of a first-order evaluation and a second-order
C     evaluations of the aerosol distribution.
C     This can be only used for second-order algorithms.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ: number of aerosol variables.
C     NBIN_AER: number of aerosol bins.
C     Q1:  first-order evaluation of the aerosol conc. ([\mu.g.m^-3])
C     Q:   second-order evaluation of the aerosol conc. ([\mu.g.m^-3])
C     IQ: aerosol pointers.
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/03/23: cleaning (Bruno Sportisse, CEREA)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Edouard Debry and Bruno Sportisse, CEREA.
C
C------------------------------------------------------------------------


      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'time.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q1(neq),q(neq)
      INTEGER iq(NEXT,nbin_aer)

      INTEGER jj
      DOUBLE PRECISION tmp,n2err,R

C     ******zero init
      n2err=0.D0

C     ******local error estimation
      DO jj=1,nbin_aer
         if(q(jj).gt.TINYN) then
            tmp=(q(jj)-q1(jj))/(q(jj)+TINYN)
            n2err=n2err+tmp*tmp
         endif
      END DO

                                ! only dynamic sizes
      IF(ICUT.NE.nbin_aer) THEN
         DO jj=IQ(E1,ICUT+1),IQ(E2,nbin_aer)
            if(q(jj).gt.TINYM) then
               tmp=(q(jj)-q1(jj))/(q(jj)+TINYM)
               n2err=n2err+tmp*tmp
            endif
         END DO
      ENDIF

      n2err=DSQRT(n2err)


C     ******compute new time step
                                ! first we constrain norm2 error
                                ! in order to prevent division by zero
                                ! and to keep new time step between
                                ! DTMIN and DTMAX defined in time.inc
      R = (1.D+02)/(1.D-05)
      tmp=R*R
      n2err=DMIN1(n2err,EPSER*tmp)
      n2err=DMAX1(EPSER/tmp,n2err)

                                ! formula to compute new time step
      DT2=DT2*DSQRT(EPSER/n2err)

      DT2 = DMIN1( (TOUT2-TIN2) , DMAX1( DTAEROMIN, DT2) )

      END

