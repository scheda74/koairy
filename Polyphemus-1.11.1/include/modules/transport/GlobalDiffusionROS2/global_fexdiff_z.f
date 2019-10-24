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



      SUBROUTINE global_FEXDIFF_Z (Nz,DLconc,DLr,DLkz,
     $     Zdm,Zdc,Dcdep,Dcemis,rhoz)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes div(K_z \partial_z c) with a three point
C     scheme.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along Z.
C     DLKz: array of K_z (diffusion coefficient) along Z.
C     ZDM: array of distances between cell interfaces along Z.
C     ZDC: array of distances between nodes along Z.
C     DCDEP: deposition velocity.
C     DCEMIS: surfacic emission.
C     RHOZ: array of rho_z (air density) along Z.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DLR: div(rho_z K_z \partial_z (c/rho_z)).
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

      INTEGER Nz, Jk
      DOUBLE PRECISION DLconc(Nz), DLr(Nz)
      DOUBLE PRECISION Zdm(Nz)
      DOUBLE PRECISION Zdc(Nz)
      DOUBLE PRECISION DLkz(Nz+1)
      DOUBLE PRECISION Dcdep,Dcemis
      DOUBLE PRECISION rhoz(Nz)
      DOUBLE PRECISION Flux(Nz+1)

      DO Jk=1,Nz+1
         IF (Jk.EQ.1) THEN      ! Boundary condition (left).
            Flux(Jk) = Dcdep*DLconc(Jk)-Dcemis
         ELSEIF (Jk.EQ.Nz+1) THEN ! Boundary condition (right).
            Flux(Jk) = 0.D0
         ELSE
            Flux(Jk) = DLkz(Jk) / Zdc(Jk)
     s           * (DLconc(Jk)/rhoz(Jk) - DLconc(Jk-1)/rhoz(Jk-1))
         ENDIF
      ENDDO

      DO Jk=1,Nz
         DLr(Jk) = (Flux(Jk+1)-Flux(Jk)) / Zdm(Jk)
      ENDDO

      END
