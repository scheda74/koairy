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

      SUBROUTINE MASSCNSRV(neq,nbin_aer,q,iq)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the gas-phase concentration by mass
C     conservation.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     iq: index of aerosol species in q(*) vector.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q : gas/aerosol concentration ([µg.m-3]).
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
      INCLUDE 'dynaero.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(NEXT,nbin_aer)

      INTEGER js,ji,jj,jesp
      DOUBLE PRECISION qaero

C     ******Compute gas conc by mass conservation.

      DO jesp=G1,G2
         ji=IG(jesp)

         qaero=0.D0
         DO js=1,nbin_aer
            jj=IQ(jesp,js)
            qaero=qaero+q(jj)
         END DO

         q(ji)=QTOT(jesp)-qaero

         IF(q(ji) .lt. 0.d0) then
            DO js = 1,nbin_aer
               q(IQ(jesp,js)) = q(IQ(jesp,js))-(dabs(q(ji)))*
     &              q(IQ(jesp,js))/qaero
            ENDDO

            qaero=0.D0
            DO js=1,nbin_aer
               jj=IQ(jesp,js)
               qaero=qaero+q(jj)
            END DO

            q(ji)=QTOT(jesp)-qaero
         ENDIF

      END DO

      END

