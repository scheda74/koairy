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

      SUBROUTINE CALCGAIN(nbin_aer,qnum,qesp,fcg,gain2,couples_coag,
     s      first_index_coag,second_index_coag)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the coagulation gain for the GDE.
C     If dc/dt=P(c)-L(c)c, this computes P(c).
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nbin_aer: number of aerosol bins.
C     QNUM: aerosol number density ([#aero.m^-3]).
C     QESP: aerosol mass density ([\mu g.m^-3]).
C     fcg: coagulation coefficients ([#aero.m^3.s^-1]).
C     GAIN2: aerosol mass flux ([\mu g.m^-3.s^-1]).
C     couples_coag: coagulation couples for each bin.
C     first_index_coag: first index of coagulation couples.
C     second_index_coag: second index of coagulation couples.
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     GAIN2: coagulation gain ([µg.m-3.s-1]).
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
C     2004: Edouard Debry, 2004.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nbin_aer
      DOUBLE PRECISION qnum(nbin_aer),qesp(nbin_aer)
      DOUBLE PRECISION gain2(nbin_aer),fcg(nbin_aer,nbin_aer,nbin_aer)
      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)

      INTEGER js,j1,j2,jl
      DOUBLE PRECISION gsum

C     ******zero init
      DO js=1,nbin_aer
         gain2(js)=0.D0
      END DO

C     ******compute coagulation gain2
      DO js=1,nbin_aer
         gsum=0.D0

         DO jl=1,couples_coag(js)

            j1=first_index_coag(js,jl)
            j2=second_index_coag(js,jl)

            IF (j1.EQ.j2) THEN
               gsum=gsum+fcg(j1,j1,js)
     &              *qnum(j1)*qesp(j1)
            ELSE
               gsum=gsum+fcg(j1,j2,js)
     &              *( qnum(j1)*qesp(j2)
     &              +qnum(j2)*qesp(j1) )
            ENDIF
         END DO

         gain2(js)=gsum
      END DO

      END

