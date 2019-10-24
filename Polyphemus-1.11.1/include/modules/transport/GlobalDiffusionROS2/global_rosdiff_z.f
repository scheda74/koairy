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



      SUBROUTINE global_rosdiff_z (Nz,DLconc,
     $     ts,tf,Dkz,Dkzf,Zdm,Zdc,
     s     Dcdep,Dcemis,Dcdepf,Dcemisf,rhoz)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along Z for one timestep.
C     The numerical scheme is a second order Rosenbrock.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     DLKZ: array of K_z (diffusion coefficient) along Z.
C     RHOZ: array of rho_z (air density) along Z.
C     TS: initial time.
C     TF: final time.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of chemical concentrations along Z.
C
C     -- OUTPUT VARIABLES
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
C     Denis Quélo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER  Nz, Ji
      DOUBLE PRECISION DLconcbis(Nz),DLconc(Nz)
      DOUBLE PRECISION DLstep,ts,tf
      DOUBLE PRECISION Zdm(Nz),Zdc(Nz)
      DOUBLE PRECISION Igamma
      DOUBLE PRECISION DLk1(Nz),DLk2(Nz)
      DOUBLE PRECISION DLb1(Nz),DLb2(Nz)
      DOUBLE PRECISION DLmatz(Nz)
      DOUBLE PRECISION DLmatzu(Nz)
      DOUBLE PRECISION DLmatzl(Nz)
      DOUBLE PRECISION DLdz(Nz)
      DOUBLE PRECISION DLdzu(Nz)
      DOUBLE PRECISION DLdzl(Nz)
      DOUBLE PRECISION Dkz(Nz+1),Dkzf(Nz+1)
      DOUBLE PRECISION Dcdep,Dcdepf
      DOUBLE PRECISION Dcemis,Dcemisf
      DOUBLE PRECISION rhoz(Nz)

      DOUBLE PRECISION DLconc_new(Nz)
      DOUBLE PRECISION Threshold

      Threshold = 0.D0

      DLstep=tf-ts
      Igamma = 1.D0+1.D0/DSQRT(2.D0)

C     -- ROSENBROCK
C     The numerical scheme reads:
C     DLCONC = DLCONC + 3/2 (TF-TS) K1 + 1/2 (TF-TS) K2
C     where K1 and K2 are computed as follows.

C------------------------------------------------------------------------
C     1 - Compute K1

C     div(K_z \partial_z c).
      CALL global_fexdiff_z (Nz,DLconc,DLb1,Dkz,Zdm,Zdc,
     $     Dcdep,Dcemis,rhoz)

C     Jacobian matrix.
      CALL global_jacddiffdc_z (Nz,DLdz,DLdzu,DLdzl,Dkz,Zdm,Zdc,Dcdep
     $     ,rhoz)

      DO Ji=1,Nz
         DLmatz(Ji) = 1.D0-Igamma*DLstep*DLdz(Ji)
         DLmatzu(Ji) = -Igamma*DLstep*DLdzu(Ji)
         DLmatzl(Ji) = -Igamma*DLstep*DLdzl(Ji)
      ENDDO

      CALL global_solvtridiag_z (Nz,DLmatzl,DLmatz,DLmatzu,DLk1,DLb1)

C------------------------------------------------------------------------
C     2 - Compute K2

      DO Ji=1,Nz
         DLconcbis(Ji) = DLconc(Ji) + DLstep * DLk1(Ji)
         IF (DLconcbis(Ji) .LE. Threshold) THEN
            DLconcbis(Ji) = Threshold
            DLk1(Ji) = (DLconcbis(Ji) - DLconc(Ji)) / DLstep
         ENDIF
      ENDDO

      CALL global_fexdiff_z (Nz,DLconcbis,DLb2,Dkzf,Zdm,Zdc,
     $     Dcdepf,Dcemisf,rhoz)

      DO Ji=1,Nz
         DLb2(Ji) =  DLb2(Ji) - 2.D0*DLk1(Ji)
      ENDDO

      CALL global_solvtridiag_z (Nz,DLmatzl,DLmatz,DLmatzu,DLk2,DLb2)

C------------------------------------------------------------------------
C     3 - Compute DLCONC

      DO Ji=1,Nz
         DLconc_new(Ji) = DLconc(Ji)+1.5D0*DLstep*DLk1(Ji)
     &        +0.5D0*DLstep*DLk2(Ji)
         IF (DLconc_new(Ji) .LE. Threshold) THEN
            DLconc(Ji) = Threshold
         ELSE
            DLconc(Ji) = DLconc_new(Ji)
         ENDIF
      ENDDO

      END

