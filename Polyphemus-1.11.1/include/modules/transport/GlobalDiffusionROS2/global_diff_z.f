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



      SUBROUTINE global_DIFF_Z (Nx,Ny,Nz,ts,tf,Rkz,Cdep,Cemis,
     s     Rkzf,Cdepf,Cemisf,Dcz,Dmz,rho,C)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves diffusion along Z for one timestep.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     TS: initial time.
C     TF: final time.
C     RKZ: 3D Kz field at initial time.
C     CDEP: 2D deposition velocities.
C     CEMIS: .
C     RKZF: 3D Kz field at final time.
C     CDEPF: .
C     CEMISF: .
C     RHO: array of rho (air density).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     C: 3D field of concentrations.
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
      DOUBLE PRECISION Dcz(Nz)
      DOUBLE PRECISION Dmz(Nz)
      DOUBLE PRECISION C(NX,NY,NZ),ZCz(NZ)
      DOUBLE PRECISION Rkz(NX,NY,NZ+1),Rkzf(NX,NY,NZ+1)
      DOUBLE PRECISION Cdep(Nx,Ny),Cdepf(Nx,Ny)
      DOUBLE PRECISION Cemis(Nx,Ny),Cemisf(Nx,Ny)
      DOUBLE PRECISION rho(NX,NY,NZ),rhoz(Nz)

      DOUBLE PRECISION Dkz(Nz+1),Dkzf(Nz+1)
      DOUBLE PRECISION Dcdep,Dcdepf
      DOUBLE PRECISION Dcemis,Dcemisf
      DOUBLE PRECISION ts,tf

      INTEGER Ji,Jj,Jk


      DO Jj=1,NY
         DO Ji=1,NX

            DO Jk=1,NZ
               ZCz(Jk)=C(Ji,Jj,Jk)
               rhoz(Jk)=rho(Ji,Jj,Jk)
            ENDDO

            DO Jk=1,NZ+1
               Dkz(Jk)=Rkz(Ji,Jj,Jk)
               Dkzf(Jk)=Rkzf(Ji,Jj,Jk)
            ENDDO

            Dcdep=Cdep(Ji,Jj)
            Dcdepf=Cdepf(Ji,Jj)
            Dcemis=Cemis(Ji,Jj)
            Dcemisf=Cemisf(Ji,Jj)

            CALL global_rosdiff_z (Nz,ZCz,ts,tf,Dkz,Dkzf,dmz,dcz,
     s           Dcdep,Dcemis,Dcdepf,Dcemisf,rhoz)

            DO Jk=1,NZ
               C(Ji,Jj,Jk)=ZCz(Jk)
            ENDDO

         ENDDO
      ENDDO

      END
