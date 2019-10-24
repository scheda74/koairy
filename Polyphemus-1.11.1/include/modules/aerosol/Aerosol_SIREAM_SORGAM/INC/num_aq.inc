C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University. It is a
C     component of the air quality modeling system Polyphemus.
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

C     Numerical setup for aqueous-phase model.

C     Tiny parameters for outputs.
      double precision TINYAQ
      parameter (TINYAQ   = 1.d-20)

C     Tiny parameters for outputs of AQOPERATOR.f.
      double precision TINYAQ2
      parameter (TINYAQ2  = 1.d-12)

C     Subcycling timestep.
      integer NITSUBAQ
      parameter (NITSUBAQ = 10)

C     For sulfate balance conservation.

C     Tolerance error
      double precision RTOLSULF
      parameter (RTOLSULF = 1.d-2)

C     Number of possible restarts if sulfate conservation not met.
      integer NITVSRM
      parameter (NITVSRM = 10)

C     For PH computation.

C     Default value if no convergence.
      double precision PHDEF
      parameter (PHDEF = 4.5d0)

C     Relative tolerance for convergence.
      double precision RTOLPH
      parameter (RTOLPH  = 1.d-5)

C     maximal number of iterations for bisection method.
      integer NIT_PH
      parameter (NIT_PH = 1000)
