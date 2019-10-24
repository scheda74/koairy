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



      SUBROUTINE rosdiff_x (Nx,Dcx,Dmx,DLmatxl,DLmatx,DLmatxu,DLkx,
     $     rhox,DLconc,ts,tf,LITERDIFFX)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along X for one timestep.
C     The numerical scheme is a second order Rosenbrock.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLKX: array of K_x (diffusion coefficient) along X.
C     RHOX: array of rho_x (air density) along X.
C     TS: initial time.
C     TF: final time.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along X.
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

      INTEGER Nx, Ji
      DOUBLE PRECISION Dcx(Nx),Dmx(Nx)
      DOUBLE PRECISION DLmatx(Nx)
      DOUBLE PRECISION DLmatxu(Nx)
      DOUBLE PRECISION DLmatxl(Nx)
      DOUBLE PRECISION DLconcbis(Nx), DLconc(Nx)
      DOUBLE PRECISION DLkx(Nx+1), rhox(Nx)
      DOUBLE PRECISION DLstep,tf,ts
      DOUBLE PRECISION Igamma
      DOUBLE PRECISION DLk1(Nx), DLk2(Nx)
      DOUBLE PRECISION DLb1(Nx), DLb2(Nx)

      DOUBLE PRECISION DLconc_new(Nx)
      DOUBLE PRECISION Threshold

      LOGICAL LITERDIFFX

      DOUBLE PRECISION DLdx(Nx),DLdxu(Nx),DLdxl(Nx)

      Threshold = 0.D0

      DLstep=tf-ts
      Igamma = 1.D0+1.D0/DSQRT(2.D0)

C     -- ROSENBROCK
C     The numerical scheme reads:
C     DLCONC = DLCONC + 3/2 (TF-TS) K1 + 1/2 (TF-TS) K2
C     where K1 and K2 are computed as follows.

C------------------------------------------------------------------------
C     1 - Compute K1

C     div(K_x \partial_x c).
      CALL fexdiff_x (Nx,Dcx,Dmx,DLconc,DLkx,rhox,DLb1)

      IF (LITERDIFFX) THEN

C     Jacobian matrix.
         CALL jacddiffdc_x(Nx,Dcx,Dmx,DLkx,rhox,DLdxl,DLdx,DLdxu)

         DO Ji=1,Nx
            DLmatx(Ji) = 1.D0-Igamma*DLstep*DLdx(Ji)
            DLmatxu(Ji) = -Igamma*DLstep*DLdxu(Ji)
            DLmatxl(Ji) = -Igamma*DLstep*DLdxl(Ji)
         ENDDO

      ENDIF

      CALL solvtridiag_x (Nx,DLmatxl,DLmatx,DLmatxu,DLk1,DLb1)

C------------------------------------------------------------------------
C     2 - Compute K2

      DO Ji=1,Nx
         DLconcbis(Ji) = DLconc(Ji) + DLstep * DLk1(Ji)
         IF (DLconcbis(Ji) .LE. Threshold) THEN
            DLconcbis(Ji) = Threshold
            DLk1(Ji) = (DLconcbis(Ji) - DLconc(Ji)) / DLstep
         ENDIF
      ENDDO

      CALL fexdiff_x (Nx,Dcx,Dmx,DLconcbis,DLkx,rhox,DLb2)

      DO Ji=1,Nx
         DLb2(Ji) =  DLb2(Ji)-2.D0*DLk1(Ji)
      ENDDO

      CALL solvtridiag_x (Nx,DLmatxl,DLmatx,DLmatxu,DLk2,DLb2)

C------------------------------------------------------------------------
C     3 - Compute DLCONC

      DO Ji=1,Nx
         DLconc_new(Ji) = DLconc(Ji)+1.5D0*DLstep*DLk1(Ji)
     &        +0.5D0*DLstep*DLk2(Ji)
         IF (DLconc_new(Ji) .LE. Threshold) THEN
            DLconc(Ji) = Threshold
         ELSE
            DLconc(Ji) = DLconc_new(Ji)
         ENDIF
      ENDDO

      END
