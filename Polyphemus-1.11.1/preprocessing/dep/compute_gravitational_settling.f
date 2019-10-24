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



      SUBROUTINE COMPUTE_GRAVITATIONAL_SETTLING(Temperature, Pressure,
     $     Density, Diameter, AirMeanFreePath, SettlingVelocity)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes gravitational settling velocities.
C     (This is the Stokes velocity.)
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Temperature: Temperature. ([K])
C     Pressure: Pressure. ([Pa])
C     Density: Particle density. ([kg/m3])
C     Diameter: Particle wet diameter. ([m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C     SettlingVelocity:  Gravitational settling velocity. ([m/s])
C     AirMeanFreePath : Air mean free path. ([m])
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- Input.
      DOUBLE PRECISION Temperature
      DOUBLE PRECISION Pressure
      DOUBLE PRECISION Density
      DOUBLE PRECISION Diameter

C     -- Local variables.
      DOUBLE PRECISION RHO,VIM
      DOUBLE PRECISION CC
      DOUBLE PRECISION VSTOKES
      DOUBLE PRECISION VG0
      DOUBLE PRECISION dVG,dFVG,FVG1,FVG2
      INTEGER I
C     Molar mass of air. ([kg.mol-1])
      DOUBLE PRECISION MMair
C     Perfect gas constant. ([J.mol-1.K-1])
      DOUBLE PRECISION RGAS

C     -- Output.
      DOUBLE PRECISION AirMeanFreePath
      DOUBLE PRECISION SettlingVelocity


C     I-Compute:
      MMair = 2.897D-02
      RGAS = 8.314D0

C     1-Air density (GP).
      RHO = Pressure/(RGAS / MMair * Temperature)

C     2-Air mean free path ([m]) and dynamic viscosity.
      call COMPUTE_AIR_FREE_MEAN_PATH(Temperature, Pressure,
     $     AirMeanFreePath, VIM)
      AirMeanFreePath = AirMeanFreePath*1.D-6

C     3-Cuningham correction.
      call compute_CC(AirMeanFreePath, Diameter, CC)

C     4-Settling velocity in Stokes's region with Cuningham correction.
      call compute_VSTOKES(Diameter,Density,CC,VIM,VSTOKES)

C     II-Compute settling velocity for all Reynolds Number.

      SettlingVelocity=VSTOKES
      VG0=0.D0

      DO I=1,100
         IF ( ABS(SettlingVelocity-VG0).GT.1.D-10 ) THEN

            VG0=SettlingVelocity
            dVG=1.D-8*SettlingVelocity
            call compute_FVG(Diameter, Density, VIM, rho, CC,
     $           SettlingVelocity+dVG, FVG1)
            call compute_FVG(Diameter, Density, VIM, rho, CC,
     $           SettlingVelocity, FVG2)

            dFVG = (FVG1-FVG2)/dVG

            SettlingVelocity=SettlingVelocity - FVG2/dFVG

         ELSE

            EXIT

         ENDIF
      ENDDO                     !Newton

      END


C------------------------------------------------------------------------
      subroutine compute_FVG(D,rhop,mu,rho,CC,U,FVG)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the function f(.) such that the gravitational
C     settling velocity for a particle (density, diameter) meets f(v)=0.
C     This equation is then solved by a Newton algorithm.
C     Pandis/Seinfeld 1998, page 467, (8:44).
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     D: particle diameter ([m]).
C     RHOP: particle density ([kg/m3]).
C     CC: Cunningham correction factor
C     MU: air viscosity ([kg/m/s]).
C     RHO: air density ([kg/m^3]).
C     U : candidate velocity ([m/s]).
C
C     -- OUTPUT VARIABLES
C
C     FVSET : evaluation of f(U) ([m/s]).
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION D,rhop,mu,rho,CC,U
      DOUBLE PRECISION CD,Re,nuair,FVG
C     Gravity acceleration. ([m/s^2])
      DOUBLE PRECISION g

      g = 9.81d0

C     Kinematic viscosity
      nuair=mu/rho

C     Particle Reynolds
      re=D*U/nuair

      call compute_CD(Re,CD)
      FVG=U-dsqrt(4.D0*g*D*CC*rhop/
     &     (3.D0*CD*rho))

      RETURN
      END
