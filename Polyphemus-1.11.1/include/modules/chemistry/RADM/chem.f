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



      SUBROUTINE chem_radm (i1, i2, Nx,Ny,Nz,nesp,nr,nrphot,nreactphot
     $     ,nemis
     $     ,nemisspecies,nzemis,Wmol,ts, DLattenuation,DLhumid,DLtemp
     $     ,DLpress,DLCsourc, DLCphotolysis_rates, delta_t
     $     ,DLattenuationf
     $     ,DLhumidf,DLtempf, DLpressf,DLCsourcf,DLCphotolysis_ratesf
     $     ,ncycle,dlon,dlat,DLconc
     $     ,option_photolysis)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes one timestep for gas-phase chemistry RADM.
C     Chemical kinetics is solved in each grid cell.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     i1: Starting index of the subdomain along x (parallel case).
C     i2: Ending index of the subdomain along x (parallel case).
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLATTENUATION: 3D cloud attenuation field at initial time.
C     DLHUMID: 3D specific humidity field at initial time ([%]).
C     DLTEMP: 3D temperature field at initial time ([K]).
C     DLPRESS: 3D pressure field at initial time ([Pa]).
C     DLCFORC: array of 2D forced chemical concentration fields at
C     # initial time ([\mu.g/m^3]).
C     DLCSOURC: array of chemical volumic emissions at initial time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATES: photochemical kinetic rates
C     # at initial time ([s^{-1}]).
C     DELTA_T: time step ([s]).
C
C     The same variables are defined at final time of the timestep.
C     'f' is then put at the end of the name.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of 3D chemical concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
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
C     2002/02/26: new treatment of sources (Jaouad Boutahar, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Denis Quélo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION ts,delta_t
      DOUBLE PRECISION tschem,tfchem

      integer i1, i2,nx,ny,nz,nesp,nr,nrphot,nemis,nzemis

      DOUBLE PRECISION DLconc(NX,NY,NZ,nesp),ZC(NESP)
      DOUBLE PRECISION DLtemp(NX,NY,NZ),DLtempf(NX,NY,NZ)
      DOUBLE PRECISION DLattenuation(NX,NY,NZ)
      DOUBLE PRECISION DLattenuationf(NX,NY,NZ)
      DOUBLE PRECISION DLhumid(NX,NY,NZ),DLhumidf(NX,NY,NZ)
      DOUBLE PRECISION DLCsourc(Nx,Ny,Nzemis,Nemis)
      DOUBLE PRECISION DLCsourcf(Nx,Ny,Nzemis,Nemis)
      DOUBLE PRECISION ZCsourc(NESP)
      DOUBLE PRECISION ZCsourcf(NESP)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)
      DOUBLE PRECISION DLpress(NX,NY,NZ),DLpressf(NX,NY,NZ)
      DOUBLE PRECISION DLCphotolysis_rates(Nx,Ny,Nz,NRphot)
      DOUBLE PRECISION DLCphotolysis_ratesf(Nx,Ny,Nz,NRphot)

      double precision dlon(nx),dlat(ny)

      integer ncycle
      double precision convers_factor(nesp)
      double precision convers_factor_jac(nesp,nesp)
      double precision Wmol(nesp)

      DOUBLE PRECISION Zangzen,Zangzenf
      DOUBLE PRECISION Zatt,Zattf

      DOUBLE PRECISION muzero,DLmuzero
      EXTERNAL muzero

      double precision Navog
      double precision pi

      integer nreactphot(nrphot)
      integer nemisspecies(nemis)

      INTEGER Jt,Ji,Jj,Jk,Jsp,i

      INTEGER option_photolysis

C     Constants.

      pi = 3.14159265358979323846D0

      Navog = 6.02213D+23

C     Conversion \mu.g/m3 to molecules/cm3.

      DO Jsp = 1, Nesp
         convers_factor(Jsp) = Navog * 1D-12 / Wmol(Jsp)
      ENDDO

      DO Jj = 1, Nesp
         DO Ji = 1, Nesp
            convers_factor_jac(Ji,Jj) = Wmol(Ji) / Wmol(Jj)
         ENDDO
      ENDDO

      DO Jsp=1,Nesp
         ZCsourc(jsp)=0.D0
         ZCsourcf(jsp)=0.d0
      ENDDO

C     Loop on grid cells.

      DO Jk=1,NZ
         DO Jj=1,NY
            DO Ji=i1,i2

C     Spatial extraction for volumic sources.

               if (jk.le.nzemis) then
                  DO Jsp=1,Nemis
                     ZCsourc(nemisspecies(jsp)+1)=DLCsourc(Ji,Jj,Jk,Jsp)
                     ZCsourcf(nemisspecies(jsp)+1)=
     $                    DLCsourcf(Ji,Jj,Jk,Jsp)
                  ENDDO
               else
                  DO Jsp=1,Nemis
                     ZCsourc(nemisspecies(jsp)+1)=0.D0
                     ZCsourcf(nemisspecies(jsp)+1)=0.D0
                  ENDDO
               endif

C     Cloud attenuation.

               Zatt = DLattenuation(Ji, Jj, Jk)
               Zattf = DLattenuationf(Ji, Jj, Jk)

C     Projection.

               DO Jsp=1,Nesp
                  ZC(Jsp) = DLconc(Ji,Jj,Jk,Jsp)
               ENDDO

C     Integration of chemistry (eventually with subcycling).

               DO Jt=1,Ncycle
                  tschem=ts+(Jt-1)*delta_t/Ncycle
                  tfchem=tschem+delta_t/Ncycle

C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f
                  DLmuzero=muzero(tschem,Dlon(ji),Dlat(jj))
                  Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
                  DLmuzero=muzero(tfchem,Dlon(ji),Dlat(jj))
                  Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)
                  CALL Kinetic_radm(nr, DLRKi,DLtemp(Ji,Jj,Jk),
     s                 DLhumid(Ji,Jj,Jk),DLpress(Ji,Jj,Jk),
     s                 Zangzen,Zatt,option_photolysis)
                  CALL Kinetic_radm(nr, DLRKf,DLtempf(Ji,Jj,Jk),
     s                 DLhumidf(Ji,Jj,Jk),DLpressf(Ji,Jj,Jk),
     s                 Zangzenf,Zatt,option_photolysis)

C     If option_photolysis is 2,
C     photolytic reactions may be read.
            IF (option_photolysis.eq.2) then
                  DO i=1,Nrphot
                     DLRKi(Nreactphot(i)+1) = Zatt *
     $                    DLCphotolysis_rates(Ji, Jj, Jk, i)
                     DLRKf(Nreactphot(i)+1) = Zattf *
     $                    DLCphotolysis_ratesf(Ji, Jj, Jk, i)
                  ENDDO
            ENDIF

                  CALL roschem_radm(nesp, nr, ZC,ZCsourc,ZCsourcf,
     $                 convers_factor, convers_factor_jac,tschem
     $                 ,tfchem,DLRki,DLRkf)
               ENDDO

C     Storage in the 3D array of chemical concentrations.

               DO i=1,NESP
                  DLconc(Ji,Jj,Jk,i) = ZC(i)
               ENDDO

            ENDDO               ! loop x.
         ENDDO                  ! loop y.
      ENDDO                     ! loop z.

      END
