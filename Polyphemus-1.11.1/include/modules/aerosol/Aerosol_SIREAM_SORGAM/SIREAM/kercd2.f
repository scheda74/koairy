C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
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

      SUBROUTINE KERCD2(nbin_aer,is1,is2,neq,q,iq,QT,
     s            QINT,QGEQ,VSW,DSW,AA,CKV,KERCD)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine computes the condensation/evaporation time
C     derivatives for the GDE along a moving grid.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     IS1 : first bin to compute.
C     IS2 : last bin to compute.
C     NEQ : number of equations.
C     Q   : gas/aerosol concentration ([µg.m-3]).
C     VSW  : wet volume aerosol concentration ([\mu.m^3.m^-3]).
C     AA   : condensation coefficient             ([m^3.s^-1]).
C     DSW  : wet aerosol diameter             ([\mu.m]).
C     QINT : internal inorganic concentration ([\mu.g.m^-3]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     QGEQ : equilibrium gas concentration    ([\mu.g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C     CKV  : kelvin coefficient                  ([]).
C     KERCD: condensation kernel ([\mu.g.s^-1])
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
C     2004: Edouard Debry, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'paraero.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'

      INTEGER nbin_aer
      INTEGER is1,is2,neq,iq(NEXT,nbin_aer)
      DOUBLE PRECISION q(neq)

      INTEGER js,jesp,jj
      DOUBLE PRECISION qext(NEXT)
      DOUBLE PRECISION qTMPA(NEXT),QINTMPA(NINTIS),QGEQTMPA(NEXT)
      DOUBLE PRECISION AATMPA(NEXT),CKVTMPA(NEXT),KERCDTMPA(NEXT)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),QGEQ(NEXT,nbin_aer)

      DOUBLE PRECISION KERCD(NEXT,nbin_aer),CKV(NEXT,nbin_aer)
      DOUBLE PRECISION AA(NEXT,nbin_aer)

C     Compute c/e kernels in each section.

      DO js=is1,is2
                                ! zero init
         DO jesp=1,NEXT
            qext(jesp)=0.D0
         END DO

                                ! test if enough aerosols and mass
         IF ( q(js).GT.TINYN.AND.
     &        QT(js).GT.TINYM ) THEN

            DO jesp=E1,E2
               jj=IQ(jesp,js)
               qext(jesp)=q(jj)
            END DO

                                ! do not forget water
            jj=IQ(EH2O,js)
            qext(EH2O)=q(jj)

                                ! we prevent evaporation when conc
                                ! are too near from zero
            DO jesp=G1,G2
               IF (qext(jesp).LE.TINYM) THEN
                  QGEQ(jesp,js)=0.D0
               ENDIF
            END DO

            Do jj=1,NEXT
               qTMPA(jj) = q(IG1+jj-1)
               QGEQTMPA(jj) = QGEQ(jj,js)
            Enddo
            Do jj=1,NINTIS
               QINTMPA(jj) = QINT(jj,js)
            Enddo

            CALL KERCOND( q(js), ! num aero conc (#aero.m-3)
     &           qext,          ! extern aero conc (µg)
     &           qTMPA,         ! bulk gas conc (µg.m-3)
     &           QINTMPA,       ! int inorg conc (µg.m-3)
     &           VSW(js),       ! vol aero conc (µm3)
     &           DSW(js),       ! wet aero diam (µm)
     &           QGEQTMPA,      ! equi gas conc (µg)
     &           AATMPA,        ! c/e kernel coef (m3.s-1)
     &           CKVTMPA,       ! kelvin effect coef (adim)
     &           KERCDTMPA )    ! c/e kernel (µg.s-1)

            DO jj=1,NEXT
               QGEQ(jj,js) = QGEQTMPA(jj)
               AA(jj,js) = AATMPA(jj)
               CKV(jj,js) = CKVTMPA(jj)
               KERCD(jj,js) = KERCDTMPA(jj)
            ENDDO

         ELSE
                                ! if nothing in aerosols
                                ! kelvin effect set to one
                                ! mass transfer variables to zero
            DO jesp=G1,G2
               CKV(jesp,js)=1.D0
               AA(jesp,js)=0.D0
               QGEQ(jesp,js)=0.D0
               KERCD(jesp,js)=0.D0
            END DO
         ENDIF
      END DO

      END

