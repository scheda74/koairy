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


      SUBROUTINE compute_dry_deposition(Date, Nsection, Nland, Diameter,
     $     Density,LUC, Zref,  FirstLevelWindModule
     $     ,Temperature, SurfaceTemperature, Pressure, SurfacePressure
     $     ,SnowHeight, RoughnessHeight, Gamma, Alpha, SmallRadius
     $     ,LargeRadius, DepositionVelocity)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes particle dry deposition velocity.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Date: Time past from the begining of the
C     first year simulated. ([s])
C     Nsection: Number of size section.
C     Nland: Number of land use coverage.
C     Diameter: Representative diameter of aerosol particles. ([m])
C     Density: Density of aerosol particles. ([Kg/m^3])
C     LUC: Land Use Coverage percentages in the cell. ([0,1])
C     Zref: Height of the first vertical node. ([m])
C     FirstLevelWindModule: First level wind module in the cell. ([m/s])
C     Temperature: Temperature in the cell. ([K])
C     SurfaceTemperature: Surface temperature in the cell. ([K])
C     Pressure: Pressure in the cell. ([Pa])
C     SurfacePressure: Surface pressure in the cell. ([Pa])
C     SnowHeight: Snow height in the cell. ([m])
C     RoughnessHeight: Roughness height in the cell. ([m])
C     Gamma: Zhang model coefficient.
C     Alpha: Zhang model coefficient.
C     SmallRadius: Characteristic radius of small receptors. ([m])
C     LargeRadius: Characteristic radius of large receptors. ([m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DepositionVelocity: Dry deposition velocity. ([m/s])
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- Input.
      DOUBLE PRECISION Date
      INTEGER Nsection,Nland
      DOUBLE PRECISION Diameter(Nsection)
      DOUBLE PRECISION Density(Nsection)
      DOUBLE PRECISION LUC(Nland)
      DOUBLE PRECISION Zref
      DOUBLE PRECISION FirstLevelWindModule
      DOUBLE PRECISION Temperature
      DOUBLE PRECISION SurfaceTemperature
      DOUBLE PRECISION Pressure
      DOUBLE PRECISION SurfacePressure
      DOUBLE PRECISION SnowHeight
      DOUBLE PRECISION RoughnessHeight(Nland,5)
      DOUBLE PRECISION Gamma(Nland)
      DOUBLE PRECISION Alpha(Nland)
      DOUBLE PRECISION SmallRadius(Nland,5)
      DOUBLE PRECISION LargeRadius(Nland,5)

C     -- Local variables.
C     Time past from the begining of the year. ([s])
      DOUBLE PRECISION Time
C     Friction velocity. ([m/s])
      DOUBLE PRECISION FrictionVelocity(Nland)
C     Aerodynamic resistance. ([s/m])
      DOUBLE PRECISION Ra(Nland)
C     Gravitationnal settling velocity. ([m/s])
      DOUBLE PRECISION SettlingVelocity(Nsection)
C     Air mean free path. ([m])
      DOUBLE PRECISION AirMeanFreePath
C     Surface resistance. ([s/m])
      DOUBLE PRECISION Rs(Nland)
C     Section index.
      INTEGER Is
C     LUC index.
      INTEGER Iland
C     Season index.
      INTEGER Iseason

C     -- Output.
      DOUBLE PRECISION DepositionVelocity(Nsection)

C----------------------------------------------------------------------
C     1 - Setup.

      DO Is = 1,Nsection
         DepositionVelocity(Is) = 0.d0
      ENDDO

C----------------------------------------------------------------------
C     2 - Computing season index.

      Time=MOD(Date,31536000.d0)

      IF((Time.LE.5097600.d0).OR.(Time.GT.26179200.d0))THEN
C     Iseason=1 winter; november, december, january,february.
         Iseason=1
      ELSEIF(Time.LE.10368000.d0)THEN
C     Iseason=2 spring; march, april.
         Iseason=2
      ELSEIF(Time.LE.20908800.d0)THEN
C     Iseason=3 summer; may june july august.
         Iseason=3
      ELSE
C     Iseason=4 fall; september october.
         Iseason=4
      ENDIF

C------------------------------------------------------------------------
C     3 - Computing aerodynamic resistance.

      CALL aero_rst(Iseason, Nland, LUC, Zref,
     $     Temperature, SurfaceTemperature, Pressure, SurfacePressure,
     $     SnowHeight, FirstLevelWindModule, RoughnessHeight,
     $     FrictionVelocity, Ra)

      DO Is=1,Nsection

C------------------------------------------------------------------------
C     5 - Computing gravitational settling velocity.

         CALL COMPUTE_GRAVITATIONAL_SETTLING(Temperature, Pressure,
     &        Density(Is), Diameter(Is), AirMeanFreePath,
     $        SettlingVelocity(IS))

C------------------------------------------------------------------------
C     6 - Computing surface resistance.

         CALL surface_rst(Nland, LUC, Iseason, Temperature, Pressure,
     $        SettlingVelocity(Is), FrictionVelocity, AirMeanFreePath,
     $        Diameter(Is), Gamma, Alpha, SmallRadius, LargeRadius,
     $        Rs)

C------------------------------------------------------------------------
C     7 - Computing Dry Deposition velocity.

         DO Iland=1,Nland
            IF (LUC(Iland).GT.0.D0) THEN
               DepositionVelocity(Is) = DepositionVelocity(Is) +
     $              LUC(Iland) / (Ra(Iland) + Rs(Iland) +
     $              Ra(Iland)*Rs(Iland)*SettlingVelocity(Is))
            ENDIF
         ENDDO                  !End of cycle on LUC.
         DepositionVelocity(Is) = DepositionVelocity(Is) +
     &        SettlingVelocity(Is)
      ENDDO                     !End of cycle on sections.

      END
