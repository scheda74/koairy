C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Karine Sartelet and Edouard Debry
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

      SUBROUTINE FNUCL(neq,nbin_aer,q,iq,dqdt)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the source terms for the system of
C     Ordinary Differential Equations defined by nucleation.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     Q   : aerosol concentration ([µg.m-3]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C
C     DQDT : time derivative ([µg.m-3.s-1]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Karine Sartelet and Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'emw.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'dynaero.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q(neq),dqdt(neq)
      INTEGER iq(NEXT,nbin_aer)

      INTEGER jj
      DOUBLE PRECISION nanucl,qanucl,na
      DOUBLE PRECISION jnucl,ntot,dpnucl
      DOUBLE PRECISION xstar,mSO4
      DOUBLE PRECISION ntotnh3,mr

C     Compute gas mass conservation

      CALL MASSCNSRV(neq,nbin_aer,q,iq)

      IF(ITERN.EQ.1) THEN            ! sulfuric-acid-ammonia-water nucl'n
C     mr should be in ppt
         mr = 10.d12 * q(IG(ENH3))/(PRES*MMair/RGAS*TEMP)
         jj=IG(ESO4)

         na= q(jj)*1.D-06       ! convert to µg.cm-3
     &        /EMW(ESO4)        ! to mol.m-3
     &        *Navog            ! to #molec.m-3

         CALL COMPUTE_TERNARY_NUCLEATION(RH,     ! relative humidity
     &        TEMP,             ! temperature (Kelvin)
     &        na,               ! gas h2so4 conc (#molec.cm-3)
     &        mr,
     &        jnucl,            ! nucleation rate (#part.cm-3.s-1)
     &        ntot,             ! number of molec of h2so4 in nucleus
     &        ntotnh3,          ! number of molec of nh3 in nucleus
     &        dpnucl )          ! nucleation diameter (nm)

                                ! nucleation rate (#part.m-3.s-1)
         jnucl=jnucl*1.D06
         jj=IQ(ESO4,1)
         dqdt(1) =dqdt(1) +jnucl ! #part.m-3.s-1
         dqdt(jj)=dqdt(jj)+jnucl*ntot/Navog*EMW(ESO4) ! µg.m-3.s-1
         jj=IQ(ENH3,1)
         dqdt(jj)=dqdt(jj)+jnucl*ntotnh3/Navog*EMW(ENH3)

      ELSE                      !sulfuric-acid-water nucl'n

C     Compute H2SO4 threshold concentration

         CALL NA_THRESHOLD_VEAHKAMAKI(RH,TEMP,nanucl) !#molec.cm-3

         qanucl= nanucl*1.D06   ! convert to #molec.m-3
     &        /Navog            ! to mol.m-3
     &        *EMW(ESO4)        ! to µg.mol-1

C     Compute nucleation kernel if qSO4 exceed qanucl

         jj=IG(ESO4)
         IF (q(jj).GE.qanucl) THEN

            na= q(jj)*1.D-06    ! convert to µg.cm-3
     &           /EMW(ESO4)     ! to mol.m-3
     &           *Navog         ! to #molec.m-3

            CALL COMPUTE_BINARY_NUCLEATION_KERNEL( RH, ! relative humidity
     &           TEMP,          ! temperature (Kelvin)
     &           na,            ! gas h2so4 conc (#molec.cm-3)
     &           jnucl,         ! nucleation rate (#part.cm-3.s-1)
     &           ntot,          ! num of molec in nucleus
     &           xstar,         ! mol fraction of h2so4
     &           dpnucl )       ! nucleation diameter (nm)

                                ! nucleation rate (#part.m-3.s-1)
            jnucl=jnucl*1.D06

                                ! h2so4 mass in nucleus (µg)
            mSO4= ntot          ! #molec
     &           /Navog         ! Avogadro number (adim)
     &           *xstar         ! mol fraction of h2so4
     &           *EMW(ESO4)     ! mol weight µg.mol-1

            jj=IQ(ESO4,1)
            dqdt(1) =dqdt(1) +jnucl ! #part.m-3.s-1
            dqdt(jj)=dqdt(jj)+jnucl*mSO4 ! µg.m-3.s-1

         ENDIF
      ENDIF

      END

