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



      SUBROUTINE ADVECTION(Nx, Ny, Nz, Dmx, Dmy, Dmz, U, V, W, CLESP,
     $     Zclx, Zcly, Zclz, dt, C)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves advection along X, Y and Z for one timestep.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     U: 3D wind field along X ([m/s]).
C     V: 3D wind field along Y ([m/s]).
C     W: 3D wind field along Z ([m/s]).
C     ZCLX: 2D field for boundary conditions along X.
C     ZCLY: 2D field for boundary conditions along Y.
C     ZCLZ: 2D field for boundary conditions along Z.
C     DT: timestep ([s]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     C: 3D field of concentrations for a given chemical species.
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
C     Vivien Mallet, CEREA, May 2003.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- INCLUDE FILES
C     PARADOM: parameters for all forced fields.
C     COMDOM: global variables for mesh and forcing fields.

      integer nx, ny, nz, clesp

      double precision dmx(nx), dmy(ny), dmz(nz)

      DOUBLE PRECISION C(NX,NY,NZ)

      DOUBLE PRECISION dt

      DOUBLE PRECISION U(NX+1,NY,NZ), V(NX,NY+1,NZ), W(NX,NY,NZ+1)

      DOUBLE PRECISION Zclx(2,Ny,Nz), Zcly(Nx,2,Nz), Zclz(Nx,Ny)

      DOUBLE PRECISION fluxxm(NX+1,NY,NZ)
      DOUBLE PRECISION fluxxp(NX+1,NY,NZ)
      DOUBLE PRECISION fluxym(NX,NY+1,NZ)
      DOUBLE PRECISION fluxyp(NX,NY+1,NZ)
      DOUBLE PRECISION fluxzm(NX,NY,NZ+1)
      DOUBLE PRECISION fluxzp(NX,NY,NZ+1)

      INTEGER Ji,Jj,Jk,j,k

      DOUBLE PRECISION nu_m, nu_p, wind
      DOUBLE PRECISION C_0, C_1, C_2, C_3
      DOUBLE PRECISION Fx_0, Fx_1, Fy_0, Fy_1, Fz_0, Fz_1

C------------------------------------------------------------------------
C     Boundary conditions.

      if (clesp.eq.1) then

         do k=1,nz
            do j=1,ny
               if (u(1,j,k).le.0.d0) zclx(1,j,k) = c(1,j,k)
               if (u(nx+1,j,k).ge.0.d0) zclx(2,j,k) = c(nx,j,k)
            enddo
         enddo

         do k=1,nz
            do j=1,nx
               if (v(j,1,k).le.0.d0) zcly(j,1,k) = c(j,1,k)
               if (v(j,ny+1,k).ge.0.d0) zcly(j,2,k) = c(j,ny,k)
            enddo
         enddo

         do k=1,ny
            do j=1,nx
               if (w(j,k,nz+1).ge.0.d0) zclz(j,k) = c(j,k,nz)
            enddo
         enddo

      else

         do k=1,nz
            do j=1,ny
               zclx(1,j,k) = c(1,j,k)
               zclx(2,j,k) = c(nx,j,k)
            enddo
         enddo

         do k=1,nz
            do j=1,nx
               zcly(j,1,k) = c(j,1,k)
               zcly(j,2,k) = c(j,ny,k)
            enddo
         enddo

         do k=1,ny
            do j=1,nx
               zclz(j,k) = c(j,k,nz)
            enddo
         enddo

      endif


C------------------------------------------------------------------------
C     Fluxes computations along x.

      DO Jk = 1, Nz
         DO Jj = 1, Ny
            DO Ji = 1, Nx+1

               wind = U(Ji, Jj, Jk)
               IF (Ji.NE.(Nx+1)) nu_m = wind * dt / Dmx(Ji)
               IF (Ji.NE.1)      nu_p = wind * dt / Dmx(Ji-1)

               IF (Ji.EQ.1) THEN
                  IF (wind.GE.0.D0) THEN
                     fluxxm(Ji,Jj,Jk) = nu_m * Zclx(1,Jj,Jk)
                  ELSE
                     fluxxm(Ji,Jj,Jk) = nu_m * C(Ji,Jj,Jk)
                  ENDIF

               ELSEIF (Ji.EQ.(Nx+1)) THEN
                  IF (wind.GE.0.D0) THEN
                     fluxxp(Ji,Jj,Jk) = nu_p * C(Ji-1,Jj,Jk)
                  ELSE
                     fluxxp(Ji,Jj,Jk) = nu_p * Zclx(2,Jj,Jk)
                  ENDIF

               ELSE

                  C_1 = C(Ji-1,Jj,Jk)
                  C_2 = C(Ji,Jj,Jk)

                  IF (wind.GE.0.D0) THEN
                     IF (Ji.EQ.2) THEN
                        C_0 = Zclx(1,Jj,Jk)
                     ELSE
                        C_0 = C(Ji-2,Jj,Jk)
                     ENDIF
                     CALL numerical_flux( nu_m, C_0, C_1, C_2, Fx_0 )

                     IF (nu_p .NE. nu_m) THEN
                        CALL numerical_flux( nu_p, C_0, C_1, C_2, Fx_1 )
                     ELSE
                        Fx_1 = Fx_0
                     ENDIF

                     fluxxm(Ji,Jj,Jk) = Fx_0
                     fluxxp(Ji,Jj,Jk) = Fx_1

                  ELSE
                     IF (Ji.EQ.Nx) THEN
                        C_3 = Zclx(2,Jj,Jk)
                     ELSE
                        C_3 = C(Ji+1,Jj,Jk)
                     ENDIF

                     CALL numerical_flux( -nu_m, C_3, C_2, C_1, Fx_0 )

                     IF (nu_p .NE. nu_m) THEN
                        CALL numerical_flux( -nu_p, C_3, C_2, C_1, Fx_1)
                     ELSE
                        Fx_1 = Fx_0
                     ENDIF

                     fluxxm(Ji,Jj,Jk) = -Fx_0
                     fluxxp(Ji,Jj,Jk) = -Fx_1

                  ENDIF

               ENDIF

            ENDDO
         ENDDO
      ENDDO

C------------------------------------------------------------------------
C     Fluxes computations along y.

      DO Jk = 1, Nz
         DO Jj = 1, Ny+1
            DO Ji = 1, Nx

               wind = V(Ji, Jj, Jk)
               IF (Jj.NE.(Ny+1)) nu_m = wind * dt / Dmy(Jj)
               IF (Jj.NE.1)      nu_p = wind * dt / Dmy(Jj-1)

               IF (Jj.EQ.1) THEN
                  IF (wind.GE.0.D0) THEN
                     fluxym(Ji,Jj,Jk) = nu_m * Zcly(Ji,1,Jk)
                  ELSE
                     fluxym(Ji,Jj,Jk) = nu_m * C(Ji,Jj,Jk)
                  ENDIF

               ELSEIF (Jj.EQ.(Ny+1)) THEN
                  IF (wind.GE.0.D0) THEN
                     fluxyp(Ji,Jj,Jk) = nu_p * C(Ji,Jj-1,Jk)
                  ELSE
                     fluxyp(Ji,Jj,Jk) = nu_p * Zcly(Ji,2,Jk)
                  ENDIF

               ELSE

                  C_1 = C(Ji,Jj-1,Jk)
                  C_2 = C(Ji,Jj,Jk)

                  IF (wind.GE.0.D0) THEN
                     IF (Jj.EQ.2) THEN
                        C_0 = Zcly(Ji,1,Jk)
                     ELSE
                        C_0 = C(Ji,Jj-2,Jk)
                     ENDIF

                     CALL numerical_flux( nu_m, C_0, C_1, C_2, Fy_0 )
                     CALL numerical_flux( nu_p, C_0, C_1, C_2, Fy_1 )

                     fluxym(Ji,Jj,Jk) = Fy_0
                     fluxyp(Ji,Jj,Jk) = Fy_1

                  ELSE
                     IF (Jj.EQ.Ny) THEN
                        C_3 = Zcly(Ji,2,Jk)
                     ELSE
                        C_3 = C(Ji,Jj+1,Jk)
                     ENDIF

                     CALL numerical_flux( -nu_m, C_3, C_2, C_1, Fy_0 )
                     CALL numerical_flux( -nu_p, C_3, C_2, C_1, Fy_1 )

                     fluxym(Ji,Jj,Jk) = -Fy_0
                     fluxyp(Ji,Jj,Jk) = -Fy_1

                  ENDIF

               ENDIF

            ENDDO
         ENDDO
      ENDDO

C------------------------------------------------------------------------
C     Fluxes computations along z.

      DO Jk = 1, Nz+1
         DO Jj = 1, Ny
            DO Ji = 1, Nx

               wind = W(Ji, Jj, Jk)
               IF (Jk.NE.(Nz+1)) nu_m = wind * dt / Dmz(Jk)
               IF (Jk.NE.1)      nu_p = wind * dt / Dmz(Jk-1)

               IF (Jk.EQ.1) THEN
                  fluxzm(Ji,Jj,Jk) = 0.D0

               ELSEIF (Jk.EQ.(Nz+1)) THEN
                  IF (wind.GE.0.D0) THEN
                     fluxzp(Ji,Jj,Jk) = nu_p * C(Ji,Jj,Jk-1)
                  ELSE
                     fluxzp(Ji,Jj,Jk) = nu_p * Zclz(Ji,Jj)
                  ENDIF

               ELSE

                  C_1 = C(Ji,Jj,Jk-1)
                  C_2 = C(Ji,Jj,Jk)

                  IF (wind.GE.0.D0) THEN
                     IF (Jk.EQ.2) THEN
                        C_0 = C_1
                     ELSE
                        C_0 = C(Ji,Jj,Jk-2)
                     ENDIF

                     CALL numerical_flux( nu_m, C_0, C_1, C_2, Fz_0 )
                     CALL numerical_flux( nu_p, C_0, C_1, C_2, Fz_1 )

                     fluxzm(Ji,Jj,Jk) = Fz_0
                     fluxzp(Ji,Jj,Jk) = Fz_1

                  ELSE
                     IF (Jk.EQ.Nz) THEN
                        C_3 = Zclz(Ji,Jj)
                     ELSE
                        C_3 = C(Ji,Jj,Jk+1)
                     ENDIF

                     CALL numerical_flux( -nu_m, C_3, C_2, C_1, Fz_0 )
                     CALL numerical_flux( -nu_p, C_3, C_2, C_1, Fz_1 )

                     fluxzm(Ji,Jj,Jk) = -Fz_0
                     fluxzp(Ji,Jj,Jk) = -Fz_1

                  ENDIF

               ENDIF

            ENDDO
         ENDDO
      ENDDO

C------------------------------------------------------------------------
C     Updates.

      DO Jk = 1, Nz
         DO Jj = 1, Ny
            DO Ji = 1, Nx

               C(Ji, Jj, Jk) = C(Ji, Jj, Jk)
     $              - Fluxxp(Ji+1,Jj,Jk) + Fluxxm(Ji,Jj,Jk)
     $              - Fluxyp(Ji,Jj+1,Jk) + Fluxym(Ji,Jj,Jk)
     $              - Fluxzp(Ji,Jj,Jk+1) + Fluxzm(Ji,Jj,Jk)

               IF (C(Ji, Jj, Jk) .LT. 0.D0) THEN
                  C(Ji, Jj, Jk) = 0.D0
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      END
