C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, PC - INRIA - EDF R&D
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
C     NS2       current number of sctns
C     NB2       current number of bnds
C
C     ICOAG     integer=1 if coagulation
C     ICOND     integer=1 if condensation
C     INUCL     integer=1 if nucleation
C
C
C     KDSLV     kind of dynamic solver
C     1=etr,2=ros2 3=ebi
C
C     IKELV     integer=1 if kelvin effect
C     IDENS     integer=1 if varying_density
C
C     ICUT      cutting size :
C     ICUT last eq size,
C     ICUT+1 first dynamic size
C
C     QTOT(C)   total aero species mass conc (µg.m-3)
C
C     RHOA      fixed aerosol_density (µg.µm-3)
C
C     ITHRM     integer ruling the number of
C     calls to thermodynamic models
C     ( i.e. routine step.f )
C     =0 means computed every time
C     >1 means computed once per call
C     to siream model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER IKELV,IDENS,KDSLV,
     &     ICOAG,ICOND,INUCL,
     &     ICUT,IREDIST,
     &     ISOURC,ITHRM,ITERN,
     &     ISULFCOND,MTHRM

      DOUBLE PRECISION RHOA,QTOT(NEXT)

      COMMON /dyn_aero/QTOT,RHOA,ICUT,
     &     IKELV,IDENS,
     &     KDSLV,ICOAG,IREDIST,
     &     ICOND,INUCL,
     &     ISOURC,ITHRM,
     &     ITERN,ISULFCOND,MTHRM
!$OMP THREADPRIVATE(/dyn_aero/)
	INCLUDE 'paraero.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
