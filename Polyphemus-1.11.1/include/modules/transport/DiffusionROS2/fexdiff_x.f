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



      SUBROUTINE FEXDIFF_X(Nx,Dcx,Dmx,DLconc,DLkx,rhox,DLr)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes div(rho_x K_x \partial_x (c/rho_x)) with a three
C     point scheme.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along X.
C     DLKX: array of K_x (diffusion coefficient) along X.
C     RHOX: array of rho_x (air density) along X.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DLR: div(rho_x K_x \partial_x (c/rho_x) ).
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
C     Denis Quélo, Jaouad Boutahar, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER Nx, Ji
      DOUBLE PRECISION Dcx(Nx), Dmx(Nx)
      DOUBLE PRECISION DLconc(Nx), DLr(Nx), DLkx(Nx+1), rhox(Nx)
      DOUBLE PRECISION Flux(Nx+1)

      DO Ji=1,Nx+1
         IF (Ji.EQ.1) THEN      ! Boundary condition (left).
            Flux(Ji) = 0.D0
         ELSEIF (Ji.EQ.Nx+1) THEN ! Boundary condition (right).
            Flux(Ji) = 0.D0
         ELSE
            Flux(Ji) = DLkx(Ji) / Dcx(Ji)
     s           * (DLconc(Ji)/rhox(Ji) - DLconc(Ji-1)/rhox(Ji-1))
         ENDIF
      ENDDO

      DO Ji=1,Nx
         DLr(Ji) = (Flux(Ji+1)-Flux(Ji)) / Dmx(Ji)
      ENDDO

      END
