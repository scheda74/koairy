C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, NPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C
C     CEREA (http://www.enpc.fr/cerea/) is a joint laboratory of ENPC
C     (http://www.enpc.fr/) and EDF R&D (http://www.edf.fr/).
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
C-----------------------------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       C
C	PARAMETERS NEEDED FOR ORGANICS 		 C
C       C
C       C
C       PSAT(C)    soa sat pres (Pascals)          C
C       DHVAP(C)   vaporization enthalpy (J.mol-1) C
C       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC
	DOUBLE PRECISION DHVAP(EARO1:EBiBmP)
	DOUBLE PRECISION PSAT(EARO1:EBiBmP)

	DATA PSAT / 5.7D-05,	! ARO1
     &              1.6D-03,	! ARO2
     &              5.0D-06,	! ALK1
     &              5.0D-06,	! OLE1
     &              4.0D-06,	! API1
     &              1.7D-04,	! API2
     &              2.5D-05,	! LIM1
     &              1.2D-04, 	! LIM2
     &              2.67D-07,   ! AnClP
     &              3.77D-03,   ! BiISO1
     &              1.74D-05,   ! BiISO2
     &              4.00D-05 /  ! BiBmP

	DATA DHVAP / 88.0D03,   ! ARO1
     &               88.0D03,	! ARO2
     &               88.0D03,	! ALK1
     &               88.0D03,	! OLE1
     &               88.0D03,	! API1
     &               88.0D03,	! API2
     &               88.0D03,	! LIM1
     &               88.0D03,	! LIM2
     &               88.0D03,   ! AnClP
     &               42.0D03,   ! BiISO1
     &               42.0D03,   ! BiISO2
     &               175.0D03 / ! BiBmP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
