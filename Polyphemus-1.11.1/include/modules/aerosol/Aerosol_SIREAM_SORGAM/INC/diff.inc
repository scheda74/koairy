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
*     *
*     PARAMETERS NEEDED FOR                      *
*     DIFFUSION COEFFICIENTS                     *
*     *
*     SIGM(*)    molecular diameter (A°)         *
*     PARM(*)    collision factor                *
*     *
*     for inorganic : Hirschfelder et al. (1954) *
*     for organic   : Tao, Y. and McMurry, P.H.  *
*     Envir. Sci. Tech. (1989)   *
*     vol 23, pp 1519-1523       *
*     *
*     for na, poa and h2o parameters set to      *
*     arbitrary high values because not used     *
*     *
**************************************************
      DOUBLE PRECISION SIGM(NEXT),PARM(NEXT)
******
      DATA SIGM/ 1.00D10,       !Dust
     &     1.00D10,             !BC
     &     1.00D10,             !Na
     &     5.500D0,             !H2SO4
     &     2.900D0,             !NH3
     &     3.300D0,             !HNO3
     &     3.339D0,             !HCl
     &     8.390D0,             !ARO1
     &     8.390D0,             !ARO2
     &     8.390D0,             !ALK1
     &     8.390D0,             !OLE1
     &     8.390D0,             !API1
     &     8.390D0,             !API2
     &     8.390D0,             !LIM1
     &     8.390D0,             !LIM2
     &     8.390D0,             !AnClP
     &     8.390D0,             !BiISO1
     &     8.390D0,             !BiISO2
     &     8.390D0,             !BiBmP
     &     1.00D10,             !POA
     &     1.00D10 /            !H2O

      DATA PARM/ 1.000D10,      !Dust
     &     1.000D10,            !BC
     &     1.000D10,            !Na
     &     7.730D01,            !H2SO4
     &     5.583D02,            !NH3
     &     4.759D02,            !HNO3
     &     3.447D02,            !HCl
     &     6.870D02,            !ARO1
     &     6.870D02,            !ARO2
     &     6.870D02,            !ALK1
     &     6.870D02,            !OLE1
     &     6.870D02,            !API1
     &     6.870D02,            !API2
     &     6.870D02,            !LIM1
     &     6.870D02,            !LIM2
     &     6.870D02,            !AnClP
     &     6.870D02,            !BiISO1
     &     6.870D02,            !BiISO2
     &     6.870D02,            !BiBmP
     &     1.000D10,            !POA
     &     1.000D10 /           !H2O
**************************************************
