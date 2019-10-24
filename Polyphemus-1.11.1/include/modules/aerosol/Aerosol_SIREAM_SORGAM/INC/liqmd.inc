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

**************************************************
*       *
*       LIST OF PURE LIQUID MASS DENSITIES         *
*       EXPRESSED IN µg.µm-3                       *
*       *
**************************************************
	DOUBLE PRECISION LMD(NEXT)
******
	DATA LMD / 2.33D-06,	!Mineral Dust (Silicium here)
     &             2.25D-06,	!Black Carbon
     &             0.97D-06,	!Na
     &             1.84D-06,	!H2SO4
     &             0.91D-06,	!NH3
     &             1.50D-06,	!HNO3
     &             1.15D-06,	!HCl
     &             1.30D-06,	!ARO1
     &             1.30D-06,	!ARO2
     &             1.30D-06,	!ALK1
     &             1.30D-06,	!OLE1
     &             1.30D-06,	!API1
     &             1.30D-06,	!API2
     &             1.30D-06,	!LIM1
     &             1.30D-06,	!LIM2
     &             1.30D-06,	!AnClP
     &             1.30D-06,	!BiISO1
     &             1.30D-06,	!BiISO2
     &             1.30D-06,	!BiBmP
     &             1.30D-06,	!POA
     &             1.00D-06 /	!H2O
**************************************************
