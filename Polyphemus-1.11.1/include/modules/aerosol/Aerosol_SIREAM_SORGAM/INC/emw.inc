C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM),
C     which is a component of the air quality modeling system
C     Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       LIST OF MOLAR WEIGHT EXPRESSED IN µg.mol-1
C       FOR EXTERNAL SPECIES
C
C       EMW(C)   molar weight of external species
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DOUBLE PRECISION EMW(NEXT)
CCCCCC
	DATA EMW / 28.0D06,	! Mineral Dust (Silicium here)
     &             12.0D06,	! Black Carbon
     &             23.0D06,	! Na
     &             98.0D06,	! H2SO4
     &             17.0D06,     ! NH3
     &             63.0D06,	! HNO3
     &             36.5D06,    	! HCl
     &             150.0D06,   	! ARO1
     &             150.0D06,    ! ARO2
     &             140.0D06,	! ALK1
     &             140.0D06,	! OLE1
     &             184.0D06,	! API1
     &             184.0D06,	! API2
     &             200.0D06,	! LIM1
     &             200.0D06,	! LIM2
     &             167.0D06,    ! AnClP
     &             118.0D06,    ! BiISO1
     &             136.0D06,    ! BiISO2
     &             236.0D06,    ! BiBmP
     &             168.5D06,    ! POA
     &             18.0D06 /	! H2O
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
