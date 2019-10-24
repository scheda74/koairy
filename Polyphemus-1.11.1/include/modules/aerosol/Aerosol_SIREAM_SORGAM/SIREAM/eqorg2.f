C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
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

      SUBROUTINE EQORG2(qtot2,qext2)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine solves the gas/aerosol organic equilibrium for
C     organics by solving the system of algebraic equations with an
C     iterative method. It provides equilibrium values for gas phase
C     and aerosols.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     QTOT2 : total concentration for gas phase and aerosol ([µg.m-3]).
C
C
C     -- INPUT/OUTPUT VARIABLES
C
C     QEXT2 : gas/aerosol concentration at equilibrium ([µg.m-3]).
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
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'emw.inc'

      DOUBLE PRECISION qtot2(NEXT),qext2(NEXT)

      INTEGER jesp,icpt
      DOUBLE PRECISION tmp,tmp2,qold2(NEXT)
      DOUBLE PRECISION cfa,cfb,cfc,q,rms,deter

C     Initialization to zero.

      DO jesp=1,NEXT
         qold2(jesp)=0.D0
      END DO

C     Solve the non linear system

      rms=1.D0
      icpt=0

                                ! convergence loop
      DO WHILE ( rms.GT.EPSEQ.AND.
     &     icpt.LE.MAXIT )

                                ! give back new values
         DO jesp=EARO1,EPOA
            qold2(jesp)=qext2(jesp)
         END DO

                                ! compute useful value
         tmp=0.D0
         DO jesp=EARO1,EPOA
            tmp=tmp+qold2(jesp)/EMW(jesp)
         END DO

                                ! compute new values, except for POA
         DO jesp=EARO1,EBiBmP
            tmp2=tmp-qold2(jesp)/EMW(jesp)

                                ! coef of 2nd degree equation
            cfa=1.D0/EMW(jesp)
            cfb=(QSAT(jesp)-qtot2(jesp))
     &           /EMW(jesp)+tmp2
            cfc=-qtot2(jesp)*tmp2

                                ! one can note cfa>0 and cfc<0
                                ! 2nd degree equation positive sol
            deter=cfb*cfb-4.D0*cfa*cfc
            IF (deter.lt.0.D0) THEN
               WRITE(6,*)'SIREAM (eqorg2.f): sqrt(<0)'
               WRITE(6,*)'cfa=',cfa
               WRITE(6,*)'cfb=',cfb
               WRITE(6,*)'cfc=',cfc
               STOP
            ENDIF

            q=-5.D-01*(cfb+DSIGN(1.D0,cfb)*DSQRT(deter))


            qext2(jesp)=DMAX1(q/cfa,cfc/q)
         END DO

                                ! compute rms between new and old values
         rms=0.D0
         DO jesp=EARO1,EBiBmP
            tmp= (qext2(jesp)-qold2(jesp))
     &           /(qold2(jesp)+TINYM)
            rms=rms+tmp*tmp
         END DO

         rms=DSQRT(rms)
         icpt=icpt+1
      END DO

      END

