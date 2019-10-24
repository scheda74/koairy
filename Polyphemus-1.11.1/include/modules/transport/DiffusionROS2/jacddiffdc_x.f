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



      SUBROUTINE jacddiffdc_x(Nx,Dcx,Dmx,DLkx,rhox,DLdxl,DLdx,DLdxu)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the jacobian matrix of the diffusion
C     scheme.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLKX: array of K_x (diffusion coefficient) along X.
C     RHOX: array of rho_x (air density) along X.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     The jacobian matrix is a tridiagonal matrix. This routine
C     computes the three diagonals DLDXL (low), DLDX (main diagonal)
C     and DLDXU (up).
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

      INTEGER Nx
      DOUBLE PRECISION Dcx(Nx), Dmx(Nx)
      DOUBLE PRECISION DLkx(Nx+1),rhox(Nx), tmp
      DOUBLE PRECISION DfluxmDc(Nx+1), DfluxpDc(Nx+1)
      INTEGER Ji

      DOUBLE PRECISION DLdx(Nx),DLdxu(Nx),DLdxl(Nx)

      DO Ji = 1, Nx+1
         IF (Ji.EQ.1) THEN      ! Boundary condition (left).
            DfluxmDc(Ji) = 0.D0
            DfluxpDc(Ji) = 0.D0
         ELSEIF (Ji.EQ.Nx+1) THEN ! Boundary condition (right).
            DfluxmDc(Ji) = 0.D0
            DfluxpDc(Ji) = 0.D0
         ELSE
            tmp = DLkx(Ji) / Dcx(Ji)
            DfluxmDc(Ji) = tmp / rhox(Ji-1)
            DfluxpDc(Ji) = tmp / rhox(Ji)
         ENDIF
      ENDDO

      DO Ji = 1, Nx
         DLdxl(Ji) = DfluxmDc(Ji) / Dmx(Ji)
         DLdxu(Ji) = DfluxpDc(Ji+1) / Dmx(Ji)
         DLdx(Ji) = - ( DfluxpDc(Ji) + DfluxmDc(Ji+1) ) / Dmx(Ji)
      ENDDO

      END
