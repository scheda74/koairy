C-----------------------------------------------------------------------
C     Copyright (C) 2001-2007, ENPC - INRIA - EDF R&D
C
C     This file is part of the air quality modeling system Polyphemus.
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



      SUBROUTINE global_rosdiff_y (Ny,Dcy,Dmy,DLmatyl,DLmaty,DLmatyu
     $     ,DLky,rhoy,DLconc,ts,tf,LITERDIFFY)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along Y for one timestep.
C     The numerical scheme is a second order Rosenbrock.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLKY: array of K_y (diffusion coefficient) along Y.
C     RHOY: array of rho_y (air density) along Y.
C     TS: initial time.
C     TF: final time.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along Y.
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
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Denis Quélo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER Ny, Ji
      DOUBLE PRECISION Dcy(Ny), Dmy(Ny)
      DOUBLE PRECISION DLmaty(Ny)
      DOUBLE PRECISION DLmatyu(Ny)
      DOUBLE PRECISION DLmatyl(Ny)
      DOUBLE PRECISION DLconcbis(Ny), DLconc(Ny)
      DOUBLE PRECISION DLky(Ny+1), rhoy(Ny)
      DOUBLE PRECISION DLstep,ts,tf
      DOUBLE PRECISION Igamma
      DOUBLE PRECISION DLk1(Ny), DLk2(Ny)
      DOUBLE PRECISION DLb1(Ny), DLb2(Ny)

      DOUBLE PRECISION DLconc_new(Ny)
      DOUBLE PRECISION Threshold

      LOGICAL LITERDIFFY

      DOUBLE PRECISION DLdy(Ny),DLdyu(Ny),DLdyl(Ny)

      Threshold = 0.D0

      DLstep=tf-ts
      Igamma = 1.D0 + 1.D0 / DSQRT(2.D0)

C     -- ROSENBROCK
C     The numerical scheme reads:
C     DLCONC = DLCONC + 3/2 (TF-TS) K1 + 1/2 (TF-TS) K2
C     where K1 and K2 are computed as follows.

C------------------------------------------------------------------------
C     1 - Compute K1

C     div(K_y \partial_y c).
      CALL global_fexdiff_y (Ny,Dcy,Dmy,DLconc,DLky,rhoy,DLb1)

      IF (LITERDIFFY) THEN

C     Jacobian matrix.
         CALL global_jacddiffdc_y(Ny,Dcy,Dmy,DLky,rhoy,DLdyl,DLdy,DLdyu)

         DO Ji=1,Ny
            DLmaty(Ji) = 1.D0-Igamma*DLstep*DLdy(Ji)
            DLmatyu(Ji) = -Igamma*DLstep*DLdyu(Ji)
            DLmatyl(Ji) = -Igamma*DLstep*DLdyl(Ji)
         ENDDO

      ENDIF

      CALL global_solvtridiag_y (Ny,DLmatyl,DLmaty,DLmatyu,DLk1,DLb1)

C------------------------------------------------------------------------
C     2 - Compute K2

      DO Ji=1,Ny
         DLconcbis(Ji) = DLconc(Ji) + DLstep * DLk1(Ji)
         IF (DLconcbis(Ji) .LE. Threshold) THEN
            DLconcbis(Ji) = Threshold
            DLk1(Ji) = (DLconcbis(Ji) - DLconc(Ji)) / DLstep
         ENDIF
      ENDDO

      CALL global_fexdiff_y (Ny,Dcy,Dmy,DLconcbis,DLky,rhoy,DLb2)

      DO Ji=1,Ny
         DLb2(Ji) =  DLb2(Ji) - 2.D0*DLk1(Ji)
      ENDDO

      CALL global_solvtridiag_y (Ny,DLmatyl,DLmaty,DLmatyu,DLk2,DLb2)

C------------------------------------------------------------------------
C     3 - Compute DLCONC

      DO Ji=1,Ny
         DLconc_new(Ji) = DLconc(Ji)+1.5D0*DLstep*DLk1(Ji)
     &        +0.5D0*DLstep*DLk2(Ji)
         IF (DLconc_new(Ji) .LE. Threshold) THEN
            DLconc(Ji) = Threshold
         ELSE
            DLconc(Ji) = DLconc_new(Ji)
         ENDIF
      ENDDO

      END
