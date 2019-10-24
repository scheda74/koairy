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



      SUBROUTINE SPLIT_ADVECTION(Ncell, Dmx, U, CLESP, Zclx, dt, C)

C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Vivien Mallet, CEREA, April 2007.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      integer ncell, clesp

      double precision dmx(ncell)

      DOUBLE PRECISION C(NCELL)

      DOUBLE PRECISION dt

      DOUBLE PRECISION U(NCELL+1)

      DOUBLE PRECISION Zclx(2)

      DOUBLE PRECISION fluxxm(NCELL+1)
      DOUBLE PRECISION fluxxp(NCELL+1)

      INTEGER Ji

      DOUBLE PRECISION nu_m, nu_p, wind
      DOUBLE PRECISION C_0, C_1, C_2, C_3
      DOUBLE PRECISION Fx_0, Fx_1

C------------------------------------------------------------------------
C     Boundary conditions.

      if (clesp.eq.1) then

         if (u(1).le.0.d0) zclx(1) = c(1)
         if (u(ncell+1).ge.0.d0) zclx(2) = c(ncell)

      else

         zclx(1) = c(1)
         zclx(2) = c(ncell)

      endif

C------------------------------------------------------------------------
C     Fluxes computations.

      DO Ji = 1, Ncell+1

         wind = U(Ji)
         IF (Ji.NE.(Ncell+1)) nu_m = wind * dt / Dmx(Ji)
         IF (Ji.NE.1) nu_p = wind * dt / Dmx(Ji-1)

         IF (Ji.EQ.1) THEN
            IF (wind.GE.0.D0) THEN
               fluxxm(Ji) = nu_m * Zclx(1)
            ELSE
               fluxxm(Ji) = nu_m * C(Ji)
            ENDIF

         ELSEIF (Ji.EQ.(Ncell+1)) THEN
            IF (wind.GE.0.D0) THEN
               fluxxp(Ji) = nu_p * C(Ji-1)
            ELSE
               fluxxp(Ji) = nu_p * Zclx(2)
            ENDIF

         ELSE

            C_1 = C(Ji-1)
            C_2 = C(Ji)

            IF (wind.GE.0.D0) THEN
               IF (Ji.EQ.2) THEN
                  C_0 = Zclx(1)
               ELSE
                  C_0 = C(Ji-2)
               ENDIF
               CALL split_numerical_flux(nu_m, C_0, C_1, C_2, Fx_0)

               IF (nu_p .NE. nu_m) THEN
                  CALL split_numerical_flux(nu_p, C_0, C_1, C_2, Fx_1)
               ELSE
                  Fx_1 = Fx_0
               ENDIF

               fluxxm(Ji) = Fx_0
               fluxxp(Ji) = Fx_1

            ELSE
               IF (Ji.EQ.Ncell) THEN
                  C_3 = Zclx(2)
               ELSE
                  C_3 = C(Ji+1)
               ENDIF

               CALL split_numerical_flux(-nu_m, C_3, C_2, C_1, Fx_0)

               IF (nu_p .NE. nu_m) THEN
                  CALL split_numerical_flux(-nu_p, C_3, C_2, C_1, Fx_1)
               ELSE
                  Fx_1 = Fx_0
               ENDIF

               fluxxm(Ji) = -Fx_0
               fluxxp(Ji) = -Fx_1

            ENDIF

         ENDIF

      ENDDO

C------------------------------------------------------------------------
C     Updates.

      DO Ji = 1, Ncell

         C(Ji) = C(Ji) - Fluxxp(Ji+1) + Fluxxm(Ji)

         IF (C(Ji) .LT. 0.D0) THEN
            C(Ji) = 0.D0
         ENDIF

      ENDDO

      END
