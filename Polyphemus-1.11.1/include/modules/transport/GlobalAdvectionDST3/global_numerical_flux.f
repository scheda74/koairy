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



      SUBROUTINE global_numerical_flux(nu, c_0, c_1, c_2, flux)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the numerical flux between cells #1 and #2.
C     It uses a third-order upwind-biased scheme with a Sweby-type flux
C     limiter. Refer to "An efficient horizontal advection scheme for
C     the modelling of global transport of constituents", W. Hundsdorfer
C     and E.J. Spee.
C
C     `````````````````flux````````
C     `````````````````---->```````
C     ----o----|----o----|----o----
C     ```c_0```````c_1```````c_2```
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nu: courant number at the interface between cells #1 and #2.
C     C_0, C_1, C_2: concentrations in cells #0, #1, #2 and #3.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     FLUX: the numerical flux between cell #1 and cell #2.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Vivien Mallet, CEREA, May 2003.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION c_0, c_1, c_2
      DOUBLE PRECISION flux

      DOUBLE PRECISION nu, diff_c1, diff_c2, sgndiff_c2
      DOUBLE PRECISION phi, m1, m2

C     Upwind.
      flux = nu * c_1

      diff_c1 = c_1 - c_0
      diff_c2 = c_2 - c_1

      sgndiff_c2 = sign(1.D0, diff_c2)
      diff_c1 = diff_c1 * sgndiff_c2

      IF (diff_c1.GT.0.D0) THEN

         diff_c2 = abs(diff_c2)

         m1 = (1.D0 - nu) / 6.D0
     s        * ( (2.D0 - nu) * diff_c2 + (1.D0 + nu) * diff_c1 )
         m2 = (1.D0 - nu) * diff_c1

         phi = min (diff_c2, m1)
         phi = min (phi * nu, m2)

         flux = flux + phi * sgndiff_c2

      ENDIF

      END
