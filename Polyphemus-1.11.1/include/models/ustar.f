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


      SUBROUTINE ustar (Richardson, FirstLevelWindModule, Zref, z0, z0t,
     $     FrictionVelocity)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This function describes the effect of wind stability.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Richardson: Richardson number.
C     FirstLevelWindModule: First level wind module. ([m/s])
C     Zref: Height of the first vertical node. ([m])
C     z0: Dynamical roughness height. ([m])
C     z0t: Thermal roughness height. ([m])
C
C     -- OUTPUT VARIABLES
C
C     FrictionVelocity: Friction velocity. ([m/s])
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- Inputs.

      DOUBLE PRECISION Richardson
      DOUBLE PRECISION FirstLevelWindModule
      DOUBLE PRECISION Zref, z0, z0t


C     -- Local variables.

      DOUBLE PRECISION Zsz0t
      DOUBLE PRECISION alu, alt
      DOUBLE PRECISION AbsRichardson
      DOUBLE PRECISION Stability
C     Constants of J.F. LOUIS's formulas.
      DOUBLE PRECISION b, d, c
      PARAMETER (b=5.d0,c=b,d=b)
C     Van Karman constant (diffusion).
      DOUBLE PRECISION Ka

C     -- Outputs.

      DOUBLE PRECISION FrictionVelocity

C------------------------------------------------------------------------
C     0. Setup

      Ka = 0.4D0
      Zsz0t = (Zref+z0t)/z0t

      alu = Ka / DLOG( (Zref+z0)/z0 )
      alt = Ka / DLOG( Zsz0t )
      AbsRichardson=dabs(Richardson)

C------------------------------------------------------------------------
C     1. Dynamical stability function

      IF(Richardson.LE.0.D0)THEN
         Stability = 1.d0 + 2.d0*b*AbsRichardson /
     $        (1.d0 + 3.d0*b*c*alu*alt*sqrt(AbsRichardson)*
     $        sqrt(1.d0-1.d0/Zsz0t)*((Zsz0t**0.33d0-1.d0)**1.5D0))
      ELSE
         Stability = 1.d0 / (1.d0 + 2.d0*b*AbsRichardson /
     $        (1.d0+d*AbsRichardson)**.5d0)
      ENDIF

C------------------------------------------------------------------------
C     2. Friction velocity

      FrictionVelocity = alu*FirstLevelWindModule*sqrt(Stability)
      FrictionVelocity = DMAX1(FrictionVelocity,0.001d0)

      END
