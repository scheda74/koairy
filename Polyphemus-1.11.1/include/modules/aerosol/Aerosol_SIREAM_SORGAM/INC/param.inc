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

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Parameters for RACM mechanism .
C
C     NS : Number of bins.
C     NB : Number of bin bounds.
C     NEORG : Number of organic species.
C     NEINORG : Number of inorganic species.
C     NILIQ, NISLD: isorropia int species.
C     NEXT : Number of external species (aerosol composition).
C     NINTIS : Number of internal species.
C
C------------------------------------------------------------------------
      INTEGER NEORG,NEINORG,NINERT,NEXT,
     &     NILIQ,NISLD,NINTIS

      PARAMETER (
     &     NEORG=13,
     &     NEINORG=5,
     &     NINERT=2,
     &     NILIQ=12,
     &     NISLD=9,
     &     NEXT=NEINORG+NEORG+NINERT+1,
     &     NINTIS=NILIQ+NISLD)
