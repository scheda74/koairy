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



      SUBROUTINE solvtridiag_x (Nx,DLmxl,DLmx,DLmxu,DLx,DLb)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Solves M * X = b where b is known and M is a known tridiagonal
C     matrix.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLMXL: sub-diagonal elements of the matrix M.
C     DLMX: main diagonal elements of the matrix M.
C     DLMXU: super-diagonal elements of the matrix M.
C     DLB: right-hand side of the equation (i.e. vector b).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DLX: output vector X = M^(-1) * b.
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
C     Denis Quélo, CEREA, October 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER Nx, Ji
      DOUBLE PRECISION DLx(Nx), DLb(Nx)
      DOUBLE PRECISION DLmx(Nx)
      DOUBLE PRECISION DLmxu(Nx)
      DOUBLE PRECISION DLmxl(Nx)
      DOUBLE PRECISION bet,gam(Nx)


      bet=DLmx(1)
      DLx(1)=DLb(1)/bet

C------------------------------------------------------------------------
C     1 - Decomposition and forward substitution

      DO Ji=1,Nx
         IF(Ji.NE.1)THEN
            gam(Ji)=DLmxu(Ji-1)/bet
            bet=DLmx(Ji)-DLmxl(Ji)*gam(Ji)
            IF (bet.EQ.0.D0) THEN
               WRITE (*,*) 'solvdiff failed'
            ENDIF
            DLx(Ji)=(DLb(Ji)-DLmxl(Ji)*DLx(Ji-1))/bet
         ENDIF
      ENDDO

C------------------------------------------------------------------------
C     2 - Back substitution

      DO Ji=Nx,1,-1
         IF(Ji.NE.Nx)THEN
            DLx(Ji)=DLx(Ji)-gam(Ji+1)*DLx(Ji+1)
         ENDIF
      ENDDO

      END
