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

      SUBROUTINE SPLINE(n,x,y,y2)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the second-order derivative of Y(*) at
C     points X(*) for cubic splines.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     
C     N  : number of points.
C     X  : points of collocation.
C     Y  : values at points.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     Y2 : second-order derivative of Y(*). 
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
      DOUBLE PRECISION x(n),y(n),y2(n)

      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(n)
      
      y2(1)=0.D0
      u(1)=0.D0

      DO i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*Y2(i-1)+2.D0
         
         y2(i)=(sig- 1.D0)/p
         u(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
     &        -(y(i)-y(i-1))/(x(i)-x(i-1))
         u(i)=u(i)*6.D0/(x(i+1)-x(i-1))
         u(i)=(u(i)-sig*u(i-1))/p
      END DO
      
      qn=0.D0
      un=0.D0
      
      
      y2(N)=(un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
      
      DO k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      END DO

      END
