C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey and Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM),
C     which is a component of the air quality modeling system
C     Polyphemus.
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
C
C     Modif : parameters for RACM2 mechanism (YS: 28/08/2009)
C-----------------------------------------------------------------------

      integer ictmNH3,ictmHNO3,ictmHCl,ictmSO2
      integer ictmH2O2,ictmHCHO,ictmHNO2,ictmO3
      integer ictmOH,ictmHO2,ictmNO3,ictmNO,ictmNO2
      integer ictmPAN,ictmH2SO4,ictmN2O5
	integer ictmCVARO1,ictmCVARO2,ictmCVALK1
	integer ictmCVOLE1,ictmCVAPI1,ictmCVAPI2
	integer ictmCVLIM1,ictmCVLIM2,ictmAnClP
        integer ictmBiISO1,ictmBiISO2,ictmBiBmP
      parameter (ictmNH3   = 1)
      parameter (ictmHNO3  = 111)
      parameter (ictmHCl   = 2)
      parameter (ictmSO2   = 4)
      parameter (ictmH2O2  = 53)
      parameter (ictmHCHO  = 128)
      parameter (ictmHNO2  = 28)
      parameter (ictmO3    = 97)
      parameter (ictmOH    = 131)
      parameter (ictmHO2   = 130)
      parameter (ictmNO3   = 129)
      parameter (ictmNO    = 125)
      parameter (ictmNO2   = 123)
      parameter (ictmPAN   = 46)
      parameter (ictmH2SO4 = 5)
      parameter (ictmN2O5  = 20)
	parameter (ictmCVARO1= 90)
	parameter (ictmCVARO2= 91)
	parameter (ictmCVALK1= 17)
	parameter (ictmCVOLE1= 52)
	parameter (ictmCVAPI1= 26)
	parameter (ictmCVAPI2= 27)
	parameter (ictmCVLIM1= 37)
	parameter (ictmCVLIM2= 38)
	parameter (ictmAnClP= 103)
	parameter (ictmBiISO1= 23)
	parameter (ictmBiISO2= 24)
	parameter (ictmBiBmP= 19)
C     parameter (ictm =)

