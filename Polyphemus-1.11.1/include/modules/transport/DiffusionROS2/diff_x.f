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



      SUBROUTINE DIFF_X(Nx,Ny,Nz,ts,tf,
     $     DKx,Dcx,Dmx,rho,ZC)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along X for one timestep.
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
      DOUBLE PRECISION DLmatx(Nx)
      DOUBLE PRECISION DLmatxu(Nx)
      DOUBLE PRECISION DLmatxl(Nx)
      DOUBLE PRECISION Dcx(Nx)
      DOUBLE PRECISION Dmx(Nx)
      DOUBLE PRECISION DKx(NX+1,NY,NZ)
      DOUBLE PRECISION ZC(NX,NY,NZ),ZCx(NX),Zkx(Nx+1)
      DOUBLE PRECISION rho(NX,NY,NZ),rhox(Nx)
      DOUBLE PRECISION ts,tf
      INTEGER Ji,Jj,Jk

      LOGICAL LITERDIFFX

      IF (Nx.NE.1) Then

         LITERDIFFX=.TRUE.

         DO Jk=1,NZ
            DO Jj=1,Ny

               DO Ji=1,Nx
                  ZCx(Ji)=ZC(Ji,Jj,Jk)
                  rhox(Ji)=rho(Ji,Jj,Jk)
               ENDDO

               DO Ji=1,Nx+1
                  Zkx(Ji)=Dkx(Ji,Jj,Jk)
               ENDDO

               CALL rosdiff_x (Nx,Dcx,Dmx,DLmatxl,DLmatx,DLmatxu,
     $              Zkx,rhox,ZCx,ts,tf,LITERDIFFX)

               LITERDIFFX=.FALSE.

               DO Ji=1,Nx
                  ZC(Ji,Jj,Jk)=ZCx(Ji)
               ENDDO

            ENDDO
         ENDDO

      ENDIF

      END
