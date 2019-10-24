C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
C     component of the air quality modeling system Polyphemus.
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

      SUBROUTINE SPLINT(n,xa,ya,y2a,x,y)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine performs a cubic spline evaluation.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     N  : size of the arrays.
C     XA : description of the grid.
C     YA : values of function along the grid.
C     Y2A: second-order time derivative along the grid.
C     X  : point at which the function has to be computed.
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     Y : function evaluation at X.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)

      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h

      klo=1
      khi=n

      DO WHILE (khi-klo.GT.1)
         k=(khi+klo)/2
         IF (xa(k).GT.x) THEN
            khi=k
         ELSE
            klo=k
         ENDIF
      END DO

      h=xa(khi)-xa(klo)
      IF (h.EQ.0.D0) THEN
         PRINT *,'Warning from splint.f: < bad xa input in splint >'
      ENDIF

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y= a*ya(klo)+b*ya(khi)
     &     +( a*(a*a-1.D0)*y2a(klo)
     &     +b*(b*b-1.D0)*y2a(khi) )
     &     *(h*h)/6.D0

      END

