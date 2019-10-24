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

      SUBROUTINE REDIST(neq,nesp_aer,nbin_aer,q,iq,QT,
     s      XBF,MBF,DBF,HSF,XBD,MBD,DBD,HSD)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine redistributes the number and mass aerosol
C     concentrations onto the fixed sectional grid.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     Q   : gas/aerosol concentration ([\mu.g.m^-3]).
C     IQ: aerosol pointers.
C     XBF: neperian logarithm fixed aerosol bound dry mass ([]).
C     MBF: fixed aerosol bound dry mass ([\mu g]).
C     DBF: fixed aerosol bound dry diameter [\mu m]).
C     HSF: fixed width of aerosol bins ([adim]).
C     XBD: neperian logarithm moving aerosol bound dry mass ([]).
C     MBD: moving aerosol bound dry mass ([\mu g]).
C     DBD: moving aerosol bound dry diameter [\mu m]).
C     HSD: moving width of aerosol bins ([adim]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     QT: total aerosol concentration per bin ([\mu g.m^-3]).
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
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C     2005/08/11: "New" redistribution (Marilyne Tombette, CEREA).
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
      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer,nesp_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER ji,jj,jesp,j1lo,j1hi,js1,js2,icpt
      DOUBLE PRECISION qnew(neq),qtnew(nbin_aer)
      DOUBLE PRECISION x2lo,x2hi,xmin,xmax
      DOUBLE PRECISION xx,tmp1,tmp2,frac

      DOUBLE PRECISION XBF(nbin_aer+1)
      DOUBLE PRECISION MBF(nbin_aer+1),DBF(nbin_aer+1)
      DOUBLE PRECISION HSF(nbin_aer)
      DOUBLE PRECISION XBD(nbin_aer+1)
      DOUBLE PRECISION MBD(nbin_aer+1),DBD(nbin_aer+1)
      DOUBLE PRECISION HSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)

C     ******zero init
      DO jj=1,neq
         qnew(jj)=0.D0
      END DO

      DO js1=1,nbin_aer
         qtnew(js1)=0.D0
      END DO

C     ******compute redistribution
      x2lo=XBD(1)

      CALL LOCATE(nbin_aer+1,XBF,x2lo,j1lo)

      DO js2=1,nbin_aer
         x2hi=XBD(js2+1)
         CALL LOCATE(nbin_aer+1,XBF,x2hi,j1hi)

         IF (q(js2).GT.TINYN) THEN
                                ! redistribute over fixed sections
            DO js1=MAX0(j1lo,1),MIN0(j1hi,nbin_aer)

               xmin=DMAX1(x2lo,XBF(js1))
               xmax=DMIN1(x2hi,XBF(js1+1))

               frac=(xmax-xmin)/HSD(js2)

               qnew(js1)=qnew(js1)+q(js2)*frac


               DO jesp=E1,E2
                  ji=IQ(jesp,js1)
                  jj=IQ(jesp,js2)

                  qnew(ji)=qnew(ji)+q(jj)*frac
               END DO
            END DO
         ENDIF

         x2lo=x2hi
         j1lo=j1hi
      END DO

C     ******check if quantities are not too small
      icpt=1                    ! to enter in the loop
      icpt = 0
      DO WHILE (icpt.GT.0)
         icpt = 0

         DO js1=1,nbin_aer
            DO jesp=E1,E2
               ji=IQ(jesp,js1)
               qtnew(js1)=qtnew(js1)+qnew(ji)
            END DO

            tmp1= qnew(js1)*( qnew(js1)-TINYN)
            tmp2=qtnew(js1)*(qtnew(js1)-TINYM)

            IF ( tmp1.LT.0.D0.OR.
     &           tmp2.LT.0.D0 ) THEN

               IF (qnew(js1).EQ.0.D0) THEN
                  WRITE(6,*)'SIREAM (redist.f): (1)Q<0 ',qnew(js1)
                  STOP
               ENDIF
               IF ((qtnew(js1)/qnew(js1)).LE.0.D0) THEN
                  WRITE(6,*)'SIREAM (redist.f): (2)Q<0 ',
     &                 qtnew(js1)/qnew(js1)
                  STOP
               ENDIF

               xx=DLOG(qtnew(js1)/qnew(js1))
               frac=(xx-XBD(js1))/HSD(js1)

               js2=js1+1
               IF (frac.LT.5.D-01) js2=js1-1

               IF (js1.EQ.1)   js2=js1+1
               IF (js1.EQ.nbin_aer) js2=js1-1

               qnew(js2)=qnew(js2)+qnew(js1)
               qtnew(js2)=qtnew(js2)+qtnew(js1)
               qnew(js1)=0.d0
               qtnew(js1)=0.d0

               DO jesp=E1,E2
                  ji=IQ(jesp,js1)
                  jj=IQ(jesp,js2)

                  qnew(jj)=qnew(jj)+qnew(ji)
                  qnew(ji)=0.d0
               END DO
               icpt=icpt+1
            ENDIF
         END DO

      END DO

C     ******turn back to concentration vector

      DO jj=1,IQ(E2,nbin_aer)
         q(jj)=qnew(jj)
      END DO

                                ! tot mass conc
      DO js1=1,nbin_aer
         QT(js1)=qtnew(js1)
      END DO

C     ******aerosol bound concentration : return to fixed sizes
      DO js1=1,nbin_aer+1
         MBD(js1)=MBF(js1)
         DBD(js1)=DBF(js1)
         XBD(js1)=XBF(js1)
      END DO

      DO js1=1,nbin_aer
         HSD(js1)=HSF(js1)
      END DO

      END
