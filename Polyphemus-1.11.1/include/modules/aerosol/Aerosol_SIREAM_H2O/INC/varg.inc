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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       VQMG(*)   gas quadratic mean speed (m.s-1)
C       DIFFG(*)  gas diffus coefficient (m2.s-1)
C       AIRFMP    air free mean path (�m)
C       QSAT(*)   soa sat conc (�g.m-3)
C       STICK(*)  gas sticking coef (adim)
C       SIGMA(*)  aero fixed surf tension (N.m-1)
C       PSATREF(C)    ref soa sat pres (Pascals)  298K
C       KPARTREF(C)  ref soa part coef (m3/microg)  298K
C       DHVAP(C)   vaporization enthalpy (J.mol-1)
C       QSATREF(*) ref soa sat conc (�g.m-3)
C       TSATREF(*) ref soa sat conc (torr)
C       SIGM(*)    molecular diameter (A�)
C	  PARM(*)    collision factor
C
C       for inorganic : Hirschfelder et al. (1954)
C       for organic   : Tao, Y. and McMurry, P.H.
C       Envir. Sci. Tech. (1989)
C       vol 23, pp 1519-1523
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DOUBLE PRECISION DIFFG(NEXT),VQMG(NEXT)
	DOUBLE PRECISION AIRFMP, KPARTREF(NEXT)
	DOUBLE PRECISION QSAT(NEXT),KPART(NEXT)
	DOUBLE PRECISION STICK(NEXT), SIGMA(NEXT)
	DOUBLE PRECISION DHVAP(NEXT), PSATREF(NEXT)
	DOUBLE PRECISION QSATREF(NEXT),DRH(NEXT)
	DOUBLE PRECISION TSATREF(NEXT)
	DOUBLE PRECISION SIGM(NEXT),PARM(NEXT)

	COMMON /gasvar/DIFFG,VQMG,QSAT,AIRFMP,DRH,
     &               KPART,STICK,SIGMA,QSATREF,
     &               DHVAP,PSATREF,KPARTREF,TSATREF,
     &               SIGM,PARM
!$OMP   THREADPRIVATE(/gasvar/)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
