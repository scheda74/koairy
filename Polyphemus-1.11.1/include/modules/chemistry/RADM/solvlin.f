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



      SUBROUTINE solvlin_radm (nesp, Kindlu,DLa,DLalu,DLx,DLb)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves DLA * DLX = DLB where DLA is an input matrix,
C     and DLB is an input vector.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     KINDLU: 0 if DLALU is not a LU factorization of DLA. If KINDLU is
C     # not zero, DLALU is assumed to be a LU factorization of DLA.
C     DLA: matrix (NESP x NESP).
C     DLB: right-hand-side vector (NEPS) of the equation to be solved.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLALU: if KINDLU is not zero, DLALU is an LU factorisation of DLA.
C     # Otherwise, on exit, DLALU is an LU factorization of DLA.
C     IPVT: pivot indices; for 1 <= i <= NESP, row i of the
C     # matrix was interchanged with row IPVT(i).
C
C     -- OUTPUT VARIABLES
C
C     DLX: solution of DLA * DLX = DLB.
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

      integer nesp

      INTEGER Kindlu
      INTEGER Ji, Jj
      DOUBLE PRECISION DLa(NESP,NESP)
      DOUBLE PRECISION DLalu(NESP,NESP)
      DOUBLE PRECISION DLx(NESP), DLb(NESP)

      DO Ji=1,NESP
         DLx(Ji)=DLb(Ji)
      ENDDO

C------------------------------------------------------------------------
C     1 - Solve DLa * Dlx = Dlb

      IF (Kindlu .EQ. 0) THEN   ! DLalu is not
                                ! an LU factorization of DLa.
         DO Jj=1,NESP
            DO Ji=1,NESP
               DLalu(Ji,Jj)=DLa(Ji,Jj)
            ENDDO
         ENDDO

         CALL LU_decompose_radm(nesp,DLalu)
         CALL LU_solve_radm(nesp,DLalu,DLx)

      ELSE                      ! DLalu is an LU factorization of DLa.

         CALL LU_solve_radm(nesp,DLalu,DLx)

      ENDIF

      END
