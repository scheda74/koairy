C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Debry Edouard
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
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

      subroutine aerosol(i1, i2, j1, j2, nx, ny, nz, nesp, LWCmin, ts,
     $     DLhumid, DLtemp, DLpress, delta_t, DLconc, noptions_aer,
     $     options_aer, nesp_aer, nbin_aer, ncycle_aer, DLLWC, DLrain,
     $     bin_bound_aer, fixed_density_aer, density_aer, couples_coag,
     $     first_index_coag, second_index_coag, coefficient_coag,
     $     vertical_interface, DLconc_aer, Wet_Deposition,
     $     Wet_Deposition_aer, pH)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes one timestep for aerosol chemistry SIREAM.
C     Chemical kinetics is solved in each grid cell.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NX: number of cells in X direction (longitude).
C     NY: number of cells in Y direction (latitude).
C     NZ: number of vertical levels.
C     NESP: number of gas species.
C     LWCmin: air liquid water content threshold.
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLHUMID: 3D specific humidity field at initial time ([%]).
C     DLTEMP: 3D temperature field at initial time ([K]).
C     DLPRESS: 3D pressure field at initial time ([Pa]).
C     DELTA_T: time step ([s]).
C     NOPTIONS_AER: number of aerosol module options.
C     OPTIONS_AER: 1D list of aerosol module options.
C     NESP_AER: number of aerosol species.
C     NBIN_AER: number of aerosol bins.
C     NCYCLE_AER: number of cycle in aerosol computation.
C     DLLWC: 3D air liquid water content ([fraction]).
C     DLRAIN: 2D field of rain rate at initial time.
C     BIN_BOUND_AER: aerosol diameters at bin bounds.
C     FIXED_DENSITY_AER: fixed aerosol density ([kg/m^3]).
C     DENSITY_AER: size variable aerosol density ([kg/m^3]).
C     COUPLES_COAG: coagulation couples for each bin.
C     FIRST_INDEX_COAG: first bin index of coagulation couples.
C     SECOND_INDEX_COAG: second bin index of coagulation couples.
C     COEFFICIENT_COAG: coagulation partition coefficient.
C     VERTICAL_INTERFACE: 1D field of vertical interface height ([m]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of 3D gas concentrations ([\mu.g/m^3]).
C     DLCONC_AER: array of 3D aerosol concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
C
C     -- OUTPUT VARIABLES
C
C     Wet_Deposition: 2D wet fluxes of gaseous species due to in-cloud
C     scavenging ([\mu.g/m^2/s]).
C     Wet_Deposition_aer: 2D wet fluxes of particulate species due to
C     in-cloud scavenging ([\mu.g/m^2/s]).
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
C     2006/09/29: updated header (Edouard Debry).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Edouard Debry, CEREA, September 2006.
C
C------------------------------------------------------------------------

      implicit none

      include 'CONST.INC'
      include 'CONST_A.INC'
      include 'param.inc'
      include 'dynaero.inc'

      double precision ts, delta_t

      integer nx, ny, nz, nesp
      integer i1, i2, j1, j2

      double precision DLconc(nx, ny, nz, nesp), DLtemp(nx, ny, nz)
      double precision DLhumid(nx, ny, nz), DLpress(nx, ny, nz)

      integer Jt, Ji, Jj, Jk, Jsp, i, Jb

      integer nesp_aer, nbin_aer, ncycle_aer
      integer noptions_aer, options_aer(noptions_aer)
      double precision DLLWC(nx, ny, nz), DLrain(nx, ny)
      double precision bin_bound_aer(nbin_aer + 1)
      double precision density_aer(nbin_aer), fixed_density_aer
      integer couples_coag(nbin_aer)
      integer first_index_coag(nbin_aer, 4 * nbin_aer)
      integer second_index_coag(nbin_aer, 4 * nbin_aer)
      double precision coefficient_coag(nbin_aer, nbin_aer, nbin_aer)
      double precision vertical_interface(nz + 1)
      double precision layer_height(nz)
      double precision DLconc_aer(nx, ny, nz, nbin_aer, nesp_aer)

      double precision Wet_Deposition(nx, ny, nesp)
      double precision Wet_Deposition_aer(nx, ny, nbin_aer, nesp_aer)

      integer ICLD, IINCLD
      double precision rho_aero(nbin_aer)

      double precision lwca(nx, ny, nz), pH(nx, ny, nz)
      double precision tschem_aer, tfchem_aer, dtchem_aer
      double precision ZA(nesp + nbin_aer * nesp_aer)
      double precision rain_rate
      double precision qscav(nesp + nbin_aer * nesp_aer)
      double precision LWCmin

      double precision XSF(nbin_aer), MSF(nbin_aer)
      double precision DSF(nbin_aer), XBF(nbin_aer + 1)
      double precision MBF(nbin_aer + 1), DBF(nbin_aer + 1)
      double precision HSF(nbin_aer)

      integer iq(NEXT, nbin_aer)

C     Set parameters.
      ICOAG = options_aer(1)
      ICOND = options_aer(2)
      INUCL = options_aer(3)
      ICLD = options_aer(4)
      IINCLD = options_aer(5)
      IKELV = options_aer(6)
      IDENS = options_aer(7)
      ICUT = options_aer(8)
      ISULFCOND = options_aer(9)
      KDSLV = options_aer(10)
      IREDIST = options_aer(11)
      ITERN = options_aer(12)
      ITHRM = options_aer(13)
      MTHRM = options_aer(14)

C     Set pointers.
      call initpoint(nbin_aer, iq)

C     Aerosol density converted from kg / m^3 to microg / microm^3.
      RHOA = fixed_density_aer * 1.d-9
      do Jb = 1, nbin_aer
         rho_aero(Jb) = density_aer(Jb) * 1.d-9
      enddo

C     Aerosol discretization converted in microm.
      do Jb = 1, nbin_aer + 1
         DBF(Jb) = bin_bound_aer(Jb) * 1.d6
         MBF(Jb) = RHOA * cst_pi6 * DBF(Jb)**3
         XBF(Jb) = dlog(MBF(Jb))
      enddo

      do Jb = 1, nbin_aer
         XSF(Jb) = (XBF(Jb) + XBF(Jb + 1)) * 5.d-1
         MSF(Jb) = dsqrt(MBF(Jb) * MBF(Jb + 1))
         DSF(Jb) = dsqrt(DBF(Jb) * DBF(Jb + 1))
      enddo

C     Width of each fixed bin.
      do Jb = 1, nbin_aer
         HSF(Jb) = XBF(Jb + 1) - XBF(Jb)
      enddo

C     Cloud liquid water content (g/m^3).
      do Jk = 1, nz
         do Jj = j1, j2
            do Ji = i1, i2
               lwca(Ji, Jj, Jk) = DLLWC(Ji, Jj, Jk) * 1.d3
     $              * DLpress(Ji, Jj, Jk) / 101325.d0 * 28.97d0 / Pr
     $              / DLtemp(Ji, Jj, Jk)
            enddo
         enddo
      enddo

C     Initialize layer_height.
      do Jk = 1, nz
         layer_height(Jk) = vertical_interface(Jk + 1) -
     $        vertical_interface(Jk)
      enddo

C     Initialize pH.
      do Jk = 1, nz
         do Jj = j1, j2
            do Ji = i1, i2
               pH(Ji, Jj, Jk) = 0.d0
            enddo
         enddo
      enddo

C     Initialize in-cloud wet fluxes in case the aqueous module needs it.
      if (ICLD.ge.1) then
         do Jsp = 1, nesp
            do Jj = j1, j2
               do Ji = i1, i2
                  Wet_Deposition(Ji, Jj, Jsp) = 0.d0
               enddo
            enddo
         enddo

         do Jsp = 1, nesp_aer
            do Jb = 1, nbin_aer
               do Jj = j1, j2
                  do Ji = i1, i2
                     Wet_Deposition_aer(Ji, Jj, Jb, Jsp) = 0.d0
                  enddo
               enddo
            enddo
         enddo
      endif

C     Loop on grid cells.
      do Jt = 1, ncycle_aer
         tschem_aer = ts + (Jt - 1) * delta_t / ncycle_aer
         tfchem_aer = tschem_aer + delta_t / ncycle_aer
         dtchem_aer = tfchem_aer - tschem_aer
         do Jk = 1, nz
            do Jj = j1, j2
               do Ji = i1, i2
                  do Jsp = 1, nesp
                     ZA(Jsp) = DLconc(Ji, Jj, Jk, Jsp)
                     if (ZA(Jsp).lt.0.d0) ZA(Jsp) = 0.d0
                  enddo

                  do Jb = 1, nbin_aer
                     do Jsp = 1, nesp_aer
                        i = nesp + (Jsp - 1) * nbin_aer + Jb
                        ZA(i) = DLconc_aer(Ji, Jj, Jk, Jb, Jsp)
                        ZA(i) = dmax1(ZA(i), TINYM)
                     enddo
                  enddo

                  if (lwca(Ji, Jj, Jk).ge.LWCmin) then

                     if (IINCLD.eq.1) then
                        rain_rate = DLrain(Ji, Jj)
                     else
                        rain_rate = 0.d0
                     endif

                     if (ICLD.eq.1) then ! Use VSRM chemical mechanism.

                        call vsrmchem(nesp, nesp_aer, nbin_aer,
     $                       rho_aero, RHOA, DBF, DSF, XBF, ZA,
     $                       DLhumid(Ji, Jj, Jk), DLpress(Ji, Jj, Jk),
     $                       DLtemp(Ji, Jj, Jk), lwca(Ji, Jj, Jk),
     $                       tschem_aer, tfchem_aer, rain_rate,
     $                       pH(Ji, Jj, Jk), qscav)

                        do Jsp = 1, nesp
                           Wet_Deposition(Ji, Jj, Jsp) =
     $                          Wet_Deposition(Ji, Jj, Jsp)
     $                          + qscav(Jsp) * layer_height(Jk)
     $                          / dtchem_aer
                        enddo

                        do Jsp = 1, nesp_aer
                           do Jb = 1, nbin_aer
                              Wet_Deposition_aer(Ji, Jj, Jb, Jsp) =
     $                             Wet_Deposition_aer(Ji, Jj, Jb, Jsp)
     $                             + qscav(nesp + nbin_aer*(Jsp-1) + Jb)
     $                             * layer_height(Jk) / dtchem_aer
                           enddo
                        enddo

                     else if(ICLD.eq.2) then ! Use "simple aqueous" chemical mechanism.

                        call simple_aqueous_module(nesp, nesp_aer,
     $                       nbin_aer, rho_aero, RHOA, DBF, DSF, XBF,
     $                       ZA, DLhumid(Ji, Jj, Jk),
     $                       DLpress(Ji, Jj, Jk), DLtemp(Ji, Jj, Jk),
     $                       lwca(Ji, Jj, Jk), tschem_aer, tfchem_aer,
     $                       rain_rate, pH(Ji, Jj, Jk), qscav)

                        do Jsp = 1, nesp
                           Wet_Deposition(Ji, Jj, Jsp) =
     $                          Wet_Deposition(Ji, Jj, Jsp)
     $                          + qscav(Jsp) * layer_height(Jk)
     $                          / dtchem_aer
                        enddo

                        do Jsp = 1, nesp_aer
                           do Jb = 1, nbin_aer
                              Wet_Deposition_aer(Ji, Jj, Jb, Jsp) =
     $                             Wet_Deposition_aer(Ji, Jj, Jb, Jsp)
     $                             + qscav(nesp + nbin_aer*(Jsp-1) + Jb)
     $                             * layer_height(Jk) / dtchem_aer
                           enddo
                        enddo

                     else

                        call aerodyn(nesp, nesp_aer, nbin_aer,
     $                       DLtemp(Ji, Jj, Jk), DLpress(Ji, Jj, Jk),
     $                       DLhumid(Ji, Jj, Jk), tschem_aer,
     $                       tfchem_aer, ZA, couples_coag,
     $                       first_index_coag, second_index_coag,
     $                       coefficient_coag, XSF, MSF, DSF, XBF, MBF,
     $                       DBF, HSF, iq)
                     endif

                  else

                     call aerodyn(nesp, nesp_aer, nbin_aer,
     $                    DLtemp(Ji, Jj, Jk), DLpress(Ji, Jj, Jk),
     $                    DLhumid(Ji, Jj, Jk), tschem_aer,
     $                    tfchem_aer, ZA, couples_coag,
     $                    first_index_coag, second_index_coag,
     $                    coefficient_coag, XSF, MSF, DSF, XBF, MBF,
     $                    DBF, HSF, iq)

                  endif

                  do Jb = 1, nbin_aer
                     do Jsp = 1, nesp_aer
                        i = nesp + (Jsp - 1) * nbin_aer + Jb
                        DLconc_aer(Ji, Jj, Jk, Jb, Jsp) = ZA(i)
                     enddo
                  enddo

                  do Jsp = 1, nesp
                     DLconc(Ji, Jj, Jk, Jsp) = ZA(Jsp)
                  enddo

               enddo            ! loop x.
            enddo               ! loop y.
         enddo                  ! loop z.

         if (ICLD.ge.1.and.IINCLD.eq.1) then
            call settling(i1, i2, j1, j2, nx, ny, nz, nesp_aer,
     $           nbin_aer, DSF, vertical_interface, DLconc_aer,
     $           lwca, tschem_aer, tfchem_aer)
         endif
      enddo

      end
