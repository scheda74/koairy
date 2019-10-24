C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
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


      SUBROUTINE SOLVLINCL (ns, Kindlu,a,alu,x,b,
     $     accl,aluccl,xccl,bccl)

C     Routine adjointe de SOLVLIN

      IMPLICIT NONE

      INTEGER ns, Kindlu
      DOUBLE PRECISION B(ns), BCCL(ns)
      DOUBLE PRECISION A(ns,NS),ACCL(NS,NS)
      DOUBLE PRECISION ALU(NS,NS), ALUCCL(NS,NS)
      DOUBLE PRECISION X(NS), XCCL(NS)

      integer i,j

C----------------------------------------------------------------
C     1. A^T z = xccl. sortie: z=xccl

      CALL LU_solve_tr_racm(NS,alu,xccl)

C----------------------------------------------------------------
C     2. bccl = bccl + z

      do i=1,NS
         bCCL(i)=bCCL(i)+xccl(i)
      enddo

C----------------------------------------------------------------
C     3. Accl = Accl + x z^T

      do j=1,NS
         do i=1,NS
            ALUCCL(i,j)=ALUCCL(i,j)-xccl(i)*x(j)
            IF (Kindlu.EQ.0) THEN
               ACCL(i,j)=ACCL(i,j)+ALUCCL(i,j)
               ALUCCL(i,j)=0.D0
            ENDIF
         enddo
      enddo

C----------------------------------------------------------------
C     4. xccl=0

      do i = 1,NS
         xccl(i) = 0.d0
      enddo

      end
