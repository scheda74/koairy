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



      SUBROUTINE global_DIFF_Y(Nx,Ny,Nz,ts,tf,
     $     DKy,Dcy,Dmy,rho,ZC)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along Y for one timestep.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     RHO: array of rho (air density).
C     TS: initial time.
C     TF: final time.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     ZC: 3D field of concentrations.
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
C     Jaouad Boutahar, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER Nx,Ny,Nz
      DOUBLE PRECISION DLmaty(Ny)
      DOUBLE PRECISION DLmatyu(Ny)
      DOUBLE PRECISION DLmatyl(Ny)
      DOUBLE PRECISION Dcy(Ny)
      DOUBLE PRECISION Dmy(Ny)
      DOUBLE PRECISION DKy(NX,NY+1,NZ)
      DOUBLE PRECISION ZC(NX,NY,NZ),ZCy(NY),Zky(NY+1)
      DOUBLE PRECISION rho(NX,NY,NZ),rhoy(Ny)
      DOUBLE PRECISION ts,tf
      INTEGER Ji,Jj,Jk

      LOGICAL LITERDIFFY


      IF (NY.NE.1) THEN

         LITERDIFFY=.TRUE.

         DO Jk=1,NZ
            DO Ji=1,NX

               DO Jj=1,NY
                  ZCy(Jj)=ZC(Ji,Jj,Jk)
                  rhoy(Jj)=rho(Ji,Jj,Jk)
               ENDDO

               DO Jj=1,NY+1
                  Zky(Jj)=Dky(Ji,Jj,Jk)
               ENDDO

               CALL global_rosdiff_y (Ny,Dcy,Dmy,DLmatyl,DLmaty,DLmatyu,
     $              Zky,rhoy,ZCy,ts,tf,LITERDIFFY)

               LITERDIFFY=.FALSE.

               DO Jj=1,Ny
                  ZC(Ji,Jj,Jk)=ZCy(Jj)
               ENDDO

            ENDDO
         ENDDO

      ENDIF

      END
