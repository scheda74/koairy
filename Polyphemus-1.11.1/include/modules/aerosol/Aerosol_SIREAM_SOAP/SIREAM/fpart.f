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

      SUBROUTINE FPART(Nbin_aer, bound_mass, bound_mass_log,
     &                 bin_width, j1, j2, xx, fx)
C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes the partition functions between two given 
C     bins on the basis of the Fernandez-Diaz algorithm.
C     The best reference is the PhD work of Edouard Debry (Chapter 7,
C     section 7.3.2.).
C     This is not exactly the partition function (the true value is
C     obtained by dividing by the widths of the coagulating bins).
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     (J1,J2): index of the coagulating bins.
C     XX     : point at which the function has to be computed.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     FX: evaluation of the function.
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

      INTEGER Nbin_aer
      DOUBLE PRECISION bound_mass(Nbin_aer+1)
      DOUBLE PRECISION bound_mass_log(Nbin_aer+1)
      DOUBLE PRECISION bin_width(Nbin_aer)

      INTEGER j1,j2
      DOUBLE PRECISION xx,fx

      INTEGER jl
      DOUBLE PRECISION m(2),tmp,mx,x(2)

      m(1)=bound_mass(j1+1)+bound_mass(j2)
      m(2)=bound_mass(j1)+bound_mass(j2+1)
      
      IF (m(1).GT.m(2)) THEN 
         tmp=m(1)
         m(1)=m(2)
         m(2)=tmp
         jl=j2
      ELSE
         jl=j1
      ENDIF

      mx=DEXP(xx)
      
C     There are three cases according to the position of the point
C     XX in the resulting bin after coagulation.
C     For the two first cases, j1 et j2 play a symmetric role.

      IF (mx.LE.m(1)) THEN
         x(1)=bound_mass(j1)/mx
         x(2)=bound_mass(j2)/mx

         fx= 2.D0*xx-bound_mass_log(j1)-bound_mass_log(j2)
     &        +DLOG( (1.D0-x(1))*(1.D0-x(2)) )

      ELSEIF (mx.GT.m(2)) THEN
         x(1)=bound_mass(j1+1)/mx
         x(2)=bound_mass(j2+1)/mx
         
         fx= bound_mass_log(j1+1)+bound_mass_log(j2+1)-2.D0*xx
     &        -DLOG( (1.D0-x(1))*(1.D0-x(2)) )
      ELSE
         x(1)=bound_mass(jl)/mx
         x(2)=bound_mass(jl+1)/mx
         
         fx= bin_width(jl) + DLOG( (1.D0-x(1))/(1.D0-x(2)) )
      ENDIF

      END
