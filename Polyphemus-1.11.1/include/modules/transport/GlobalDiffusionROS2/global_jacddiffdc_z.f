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



      SUBROUTINE global_jacddiffdc_z (Nz,DLdz,DLdzu,DLdzl,DLkz,
     $     Zdm,Zdc,Dcdep,rhoz)

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
C     DLKz: array of K_z (diffusion coefficient) along Z.
C     ZDM: array of distances between cell interfaces along Z.
C     ZDC: array of distances between nodes along Z.
C     DCDEP: deposition velocity.
C     RHOZ: array of rho_z (air density) along Z.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DLDZ: main diagonal elements of the jacobian matrix.
C     DLDZU: super-diagonal elements of the jacobian matrix.
C     DLDZL: sub-diagonal elements of the jacobian matrix.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     The jacobian matrix is a tridiagonal matrix. This routine
C     computes the three diagonals DLDZL (low), DLDZ (main diagonal)
C     and DLDZU (up).
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

      INTEGER Nz, Jk

      DOUBLE PRECISION DLdz(Nz)
      DOUBLE PRECISION DLdzu(Nz)
      DOUBLE PRECISION DLdzl(Nz)
      DOUBLE PRECISION Zdm(Nz)
      DOUBLE PRECISION Zdc(Nz)
      DOUBLE PRECISION DLkz(Nz+1)
      DOUBLE PRECISION Dcdep
      DOUBLE PRECISION rhoz(Nz)
      DOUBLE PRECISION DfluxmDc(Nz+1), DfluxpDc(Nz+1), tmp

      DO Jk = 1, Nz+1
         IF (Jk.EQ.1) THEN      ! Boundary condition (left).
            DfluxmDc(Jk) = 0.D0
            DfluxpDc(Jk) = Dcdep
         ELSEIF (Jk.EQ.Nz+1) THEN ! Boundary condition (right).
            DfluxmDc(Jk) = 0.D0
            DfluxpDc(Jk) = 0.D0
         ELSE
            tmp = DLkz(Jk) / Zdc(Jk)
            DfluxmDc(Jk) = tmp / rhoz(Jk-1)
            DfluxpDc(Jk) = tmp / rhoz(Jk)
         ENDIF
      ENDDO

      DO Jk = 1, Nz
         DLdzl(Jk) = DfluxmDc(Jk) / Zdm(Jk)
         DLdzu(Jk) = DfluxpDc(Jk+1) / Zdm(Jk)
         DLdz(Jk) = - ( DfluxpDc(Jk) + DfluxmDc(Jk+1) ) / Zdm(Jk)
      ENDDO

      END
