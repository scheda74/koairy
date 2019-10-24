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

      SUBROUTINE SIZEBND(nbin_aer,XSF,XBF,XSD,XBD,MBD,DBD,HSD)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes size bounds.     
C     All the variables (I/O) are in common blocks.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     XSF: logarithm of fixed mean aerosol mass in bins. ([]) 
C     XBF: logarithm of fixed aerosol mass at bounds. ([])
C     
C     -- INPUT/OUTPUT VARIABLES
C 
C     XSD: logarithm of moving mean aerosol mass in bins. ([])
C     XBD: logarithm of moving aerosol mass at bounds. ([])
C     HSD: logarithm width of each bin ([])
C     MBD: moving aerosol mass at bounds. ([\mu.g])
C     DBD: moving aerosol diameter at bounds. ([\mu.m])
C     
C     
C     -- OUTPUT VARIABLES
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

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'

      INTEGER nbin_aer,jb,js
      DOUBLE PRECISION dx(nbin_aer),d2x(nbin_aer)
      DOUBLE PRECISION dxb(nbin_aer+1)

      DOUBLE PRECISION XSF(nbin_aer),XBF(nbin_aer+1)
      DOUBLE PRECISION XSD(nbin_aer),XBD(nbin_aer+1)
      DOUBLE PRECISION MBD(nbin_aer+1),DBD(nbin_aer+1)
      DOUBLE PRECISION HSD(nbin_aer)

C     ******zero init
      DO js=1,nbin_aer
         dx(js)  =0.D0
         d2x(js)=0.D0
      END DO

      DO jb=1,nbin_aer+1
         dxb(jb)=0.D0
      END DO

C     ******bound mass calculation 
      DO js=1,nbin_aer
         dx(js)=XSD(js)-XSF(js)
      END DO

      CALL SPLINE(nbin_aer,XSF,dx,d2x)

      DO jb=1,nbin_aer+1
         CALL SPLINT(nbin_aer,XSF,dx,d2x,
     &        XBF(jb),dxb(jb))
         
         XBD(jb)=XBF(jb)+dxb(jb)
      END DO

C     ******nucleation bound
      IF (INUCL.EQ.1) XBD(1)=XBF(1)
                                ! the first size is kept equal to that of
                                ! nucleation : an improvement would be 
                                ! to use the computed size given by
                                ! nucleation routines, but coagulation
                                ! coefficients should be re-computed

C     ******compute mass and diameter
      DO jb=1,nbin_aer+1
         MBD(jb)=DEXP(XBD(jb))
         DBD(jb)=(MBD(jb)/RHOA/cst_pi6)**cst_FRAC3
      END DO

C     ******sectional width
      DO js=1,nbin_aer
         HSD(js)=XBD(js+1)-XBD(js)
      END DO

      END
