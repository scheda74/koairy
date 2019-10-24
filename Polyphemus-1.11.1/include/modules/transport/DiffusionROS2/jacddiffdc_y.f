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



      SUBROUTINE jacddiffdc_y(Ny,Dcy,Dmy,DLky,rhoy,DLdyl,DLdy,DLdyu)

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
C     DLKY: array of K_y (diffusion coefficient) along Y.
C     RHOY: array of rho_y (air density) along Y.
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
C     computes the three diagonals DLDYL (low), DLDY (main diagonal)
C     and DLDYU (up).
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

      INTEGER Ny
      DOUBLE PRECISION Dcy(Ny), Dmy(Ny)
      DOUBLE PRECISION DLky(Ny+1),rhoy(Ny), tmp
      DOUBLE PRECISION DfluxmDc(Ny+1), DfluxpDc(Ny+1)
      INTEGER Jj

      DOUBLE PRECISION DLdy(Ny),DLdyu(Ny),DLdyl(Ny)

      DO Jj = 1, Ny+1
         IF (Jj.EQ.1) THEN      ! Boundary condition (left).
            DfluxmDc(Jj) = 0.D0
            DfluxpDc(Jj) = 0.D0
         ELSEIF (Jj.EQ.Ny+1) THEN ! Boundary condition (right).
            DfluxmDc(Jj) = 0.D0
            DfluxpDc(Jj) = 0.D0
         ELSE
            tmp = DLky(Jj) / Dcy(Jj)
            DfluxmDc(Jj) = tmp / rhoy(Jj-1)
            DfluxpDc(Jj) = tmp / rhoy(Jj)
         ENDIF
      ENDDO

      DO Jj = 1, Ny
         DLdyl(Jj) = DfluxmDc(Jj) / Dmy(Jj)
         DLdyu(Jj) = DfluxpDc(Jj+1) / Dmy(Jj)
         DLdy(Jj) = - ( DfluxpDc(Jj) + DfluxmDc(Jj+1) ) / Dmy(Jj)
      ENDDO

      END
