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



      SUBROUTINE FEXDIFF_Y(Ny,Dcy,Dmy,DLconc,DLky,rhoy,DLr)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes div(rho_y K_y \partial_y (c/rho_y)) with a three
C     point scheme.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along Y.
C     DLKY: array of K_y (diffusion coefficient) along Y.
C     RHOY: array of rho_y (air density) along Y.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DLR: div(rho_y K_y \partial_y (c/rho_y) ).
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

      INTEGER Ny, Jj
      DOUBLE PRECISION Dcy(Ny), Dmy(Ny)
      DOUBLE PRECISION DLconc(Ny), DLr(Ny), DLky(Ny+1), rhoy(Ny)
      DOUBLE PRECISION Flux(Ny+1)

      DO Jj=1,Ny+1
         IF (Jj.EQ.1) THEN      ! Boundary condition (left).
            Flux(Jj) = 0.D0
         ELSEIF (Jj.EQ.Ny+1) THEN ! Boundary condition (right).
            Flux(Jj) = 0.D0
         ELSE
            Flux(Jj) = DLky(Jj) / Dcy(Jj)
     s           * (DLconc(Jj)/rhoy(Jj) - DLconc(Jj-1)/rhoy(Jj-1))
         ENDIF
      ENDDO

      DO Jj=1,Ny
         DLr(Jj) = (Flux(Jj+1)-Flux(Jj)) / Dmy(Jj)
      ENDDO

      END
