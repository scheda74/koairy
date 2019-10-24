C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007 CEREA
C     Author(s): Marilyne Tombette
C
C     CEREA (http://www.enpc.fr/cerea/) is a joint laboratory of ENPC
C     (http://www.enpc.fr/) and EDF R&D (http://www.edf.fr/).
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University,
C     which is a component of the air quality modeling system Polyphemus.
C
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C
C     For more information, please see the Polyphemus web site:
C            http://www.enpc.fr/cerea/polyphemus/
C-----------------------------------------------------------------------

      subroutine redistaq_simple(nbin_aer, DSF_AERO, fixed_rho_aero,
     $     dnew, aerosol, number)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine redistributes the aerosol mass and number
C     concentration over the fixed aerosol sections.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     nbin_aer: number of section of the aerosol distribution. ([\#])
C     DSF_AERO: center diameters of the fixed sections. ([\mu m])
C     fixed_rho_aero: fixed aerosol density. ([\mu g/m^3])
C     dnew: center diameters of the sections. ([\mu m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     aerosol: aerosol concentration ([\mu g/m^3]).
C     number: number concentration ([\#/m^3]).
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
C     2007/10 Marilyne Tombette, CEREA,
C     2007/12 Cleaning and correct, Yelva Roustan, CEREA.
C
C------------------------------------------------------------------------

        IMPLICIT NONE

        include 'aerpar_simple.inc'
        include 'droppar_simple.inc'
        include 'paraero.inc'
        include 'CONST.INC'
        include 'CONST_A.INC'

        integer nbin_aer
        double precision aerosol(nbin_aer, Naers)
        double precision number(nbin_aer)
        double precision aerosol_redist(nbin_aer, Naers)
        double precision number_redist(nbin_aer)

        integer jsec, j
        integer ispe, i
        double precision dnew(nbin_aer)
        double precision aeroorig(Naers), aeronew(Naers)
        double precision DSF_AERO(nbin_aer)
        double precision fixed_rho_aero
        double precision qtotal
        double precision weight1, weight2

c     Initialize concentration (mass and number)
        do ispe = 1,Naers
           do jsec = 1,nbin_aer
              aerosol_redist(jsec, ispe) = 0.d0
           enddo
        enddo

        do jsec = 1,nbin_aer
           number_redist(jsec) = 0.d0
        end do

        do jsec = 1,nbin_aer

           call locate(nbin_aer, DSF_AERO, dnew(jsec), j)
                                ! redistribute over fixed sections
           if (j.eq.0) then
              qtotal = 0.d0
              do ispe = 1,Naers
                 aerosol_redist(1, ispe) = aerosol_redist(1, ispe)
     $                + aerosol(jsec, ispe)
                                ! compute total mass in the section
                 qtotal = qtotal + aerosol_redist(1, ispe)
              enddo
                                ! The number concentration is then modified.
                                ! This should not be the case if the mass
                                ! and the number are followed.
                                ! (the diameter could be in this case?)
              number_redist(1) = qtotal / (cst_pi6 * fixed_rho_aero)
     $             / DSF_AERO(1)**3.d0
           elseif (j.eq.nbin_aer) then
              qtotal = 0.d0
              do ispe = 1,Naers
                 aerosol_redist(nbin_aer, ispe) =
     $              aerosol_redist(nbin_aer, ispe) + aerosol(jsec, ispe)
                                ! compute total mass in the section
                 qtotal = qtotal + aerosol_redist(nbin_aer, ispe)
              enddo
                                ! The number concentration is then modified.
                                ! This should not be the case if the mass
                                ! and the number are followed.
                                ! (the diameter could be in this case?)
              number_redist(nbin_aer) = qtotal /
     $             (cst_pi6 * fixed_rho_aero) / DSF_AERO(nbin_aer)**3.d0
           else
              weight1 = ( 1 - (DSF_AERO(j+1) / dnew(jsec))**3.d0 )
     $             / (1 - (DSF_AERO(j+1) / DSF_AERO(j))**3.d0)
              weight2 = (1 - (DSF_AERO(j) / dnew(jsec))**3.d0)
     $             / (1 - (DSF_AERO(j) / DSF_AERO(j+1))**3.d0)

              do ispe = 1,Naers
                 aerosol_redist(j, ispe) = aerosol_redist(j, ispe)
     $                + aerosol(jsec, ispe) * weight1
                 aerosol_redist(j+1, ispe) = aerosol_redist(j+1, ispe)
     $                + aerosol(jsec, ispe) * weight2
              enddo

              number_redist(j) = number_redist(j)
     $             + number(jsec)
     $             * (dnew(jsec)**3.d0 - DSF_AERO(j+1)**3.d0)
     $             / (DSF_AERO(j)**3.d0 - DSF_AERO(j+1)**3.d0)
              number_redist(j+1) = number_redist(j+1)
     $             + number(jsec)
     $             * (dnew(jsec)**3.d0 - DSF_AERO(j)**3.d0)
     $             / (DSF_AERO(j+1)**3.d0 - DSF_AERO(j)**3.d0)

           endif
        enddo

C     Turn back to conc vector.

        do ispe = 1,Naers
           do jsec = 1,nbin_aer
	      aerosol(jsec, ispe) = aerosol_redist(jsec, ispe)
	   enddo
	enddo

        do jsec = 1,nbin_aer
           number(jsec) = number_redist(jsec)
        enddo

	end


