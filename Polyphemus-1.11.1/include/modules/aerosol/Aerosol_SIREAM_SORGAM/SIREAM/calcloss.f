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

      SUBROUTINE CALCLOSS(nbin_aer,qnum,kercg,loss)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the coagulation loss coefficient for the GDE.
C     If dc/dt=P(c)-L(c)c, this computes L(c).
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nbin_aer: number of aerosol bins.
C     QNUM: aerosol number density ([#part.m^-3]).
C     kercg: aerosol coagulation kernel ([m^3.s^-1]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     LOSS: coagulation loss coefficient ([s^-1]).
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

      INTEGER nbin_aer
      DOUBLE PRECISION qnum(nbin_aer),loss(nbin_aer)
      DOUBLE PRECISION kercg(nbin_aer,nbin_aer)

      INTEGER js,jk
      DOUBLE PRECISION lsum

C     ******zero init
      DO js=1,nbin_aer
         loss(js)=0.D0
      END DO

C     ******compute loss coefficient
      DO js=1,nbin_aer
         lsum=0.D0

         DO jk=1,nbin_aer
            lsum=lsum+KERCG(js,jk)*qnum(jk)
         END DO

         loss(js)=lsum
      END DO

      END

