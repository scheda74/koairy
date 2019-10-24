C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey and Bruno Sportisse
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

      SUBROUTINE ebiconc(neq,nbin_aer,q1,q,iq,QT,
     s      XSF,MSF,DSF,XSD,MSD,DSD)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine solves the lagrangian condensation/evaporation
C     equations together with nucleation with a simplified version of
C     the EBI (Euler Backward Iterative scheme).
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     XSF: neperian logarithm of fixed aerosol bin mass ([adim]).
C     MSF: fixed aerosol bin dry mass ([\mu g]).
C     DSF: fixed aerosol bin dry diameter ([\mu m]).
C     XSD: neperian logarithm of moving aerosol bin mass ([adim]).
C     MSD: moving aerosol bin dry mass ([\mu g]).
C     DSD: moving aerosol bin dry diameter ([\mu m]).
C     IQ: aerosol pointers.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q : gas/aerosol concentrations ([µg.m-3]).
C     QT: total aerosol concentration per bin ([\mu g.m^-3]).
C
C     -- OUTPUT VARIABLES
C
C     Q1: first-order arbitrary aerosol concentation ([µg.m-3])
C     # (needed for adaptive timestepping, not to be used in this case).
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C     2005/11/21: Add nucleation (Karine Sartelet, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     2004: Kathleen Fahey and Bruno Sportisse, CEREA.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'time.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'dynaero.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q(neq),q1(neq)
      INTEGER iq(NEXT,nbin_aer)

      INTEGER jj,ii,jesp,js,GLOC
      DOUBLE PRECISION kcpos(neq),kcneg(neq)

      DOUBLE PRECISION qaero(NEXT)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      DOUBLE PRECISION VSW(nbin_aer),MSW(nbin_aer)
      DOUBLE PRECISION QINT(NINTIS,nbin_aer),DSW(nbin_aer)
      DOUBLE PRECISION QGEQ(NEXT,nbin_aer),RHOW(nbin_aer)

      DOUBLE PRECISION KERCD(NEXT,nbin_aer),CKV(NEXT,nbin_aer)
      DOUBLE PRECISION AA(NEXT,nbin_aer)

      IF(ISULFCOND.EQ.1) THEN
         GLOC = ENH3            ! Do not take sulfate into account
      ELSE
         GLOC = G1
      ENDIF

C     compute inital aerosol mass per species
      DO jesp=G1,G2
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         IF (jesp.NE.ECl) THEN
#endif
            qaero(jesp)=0.D0
            DO js=1,nbin_aer
               jj=IQ(jesp,js)
               qaero(jesp)=qaero(jesp)+q(jj)
            END DO
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         ENDIF
#endif
      END DO

      DO ii=1,NITEBI

C     Initialise derivatives to zero
         DO js=1,neq
            kcpos(js) = 0.D0
            kcneg(js) = 0.D0
         ENDDO

         IF(ICUT2.NE.0.AND.INU2.EQ.1) THEN
C     Nucleation should be solved alone
C     Condensation is solved by bulkequi
            call FNUCL(neq,nbin_aer,q,iq,kcpos)

C     Number
            jj = 1
            q(jj)=(q(jj)+kcpos(jj)*dt2)
            q1(jj)=2.d0*q(jj)

C    Sulfate
            jj = IQ(ESO4,1)
            q(jj)=(q(jj)+kcpos(jj)*dt2)
            q1(jj)=2.d0*q(jj)

C     Ammonium
            jj = IQ(ENH3,1)
            q(jj)=(q(jj)+kcpos(jj)*dt2)
            q1(jj)=2.d0*q(jj)


         ELSE

            DO js=nbin_aer,ICUT2+1,-1

C     compute c/e kernels in dynamic sections
               CALL EQPART(nbin_aer,js,js,neq,q,iq,QT,
     s                  XSF,MSF,DSF,XSD,MSD,DSD,
     s                  QINT,QGEQ,VSW,MSW,DSW,RHOW)

               IKV2=IKELV
               CALL KERCD2(nbin_aer,js,js,neq,q,iq,QT,
     s                  QINT,QGEQ,VSW,DSW,AA,CKV,KERCD)

               DO jesp=GLOC,G2
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                  IF (jesp.NE.ECl) THEN
#endif
                     jj=IQ(jesp,js)

C     replace q(IG(jesp)) by qtot(jesp)-qaero(jesp)+q(jj) in kcpos
C     replace qgeq/q(jj) by qgeq/q(jj)+1 in kcneg
                     kcpos(jj) = aa(jesp,js) *
     &                    (qtot(jesp)-qaero(jesp)+q(jj))*q(js)
                     IF (q(jj).GT.tinym) THEN
                        kcneg(jj) = aa(jesp,js)*ckv(jesp,js)
     &                       *(qgeq(jesp,js)/q(jj)+1.D0)*q(js)
                     ENDIF
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                  ENDIF
#endif
               ENDDO

               IF(js.EQ.1.AND.INU2.EQ.1) then
                  call FNUCL(neq,nbin_aer,q,iq,kcpos)
C     The number of aerosols is modified in the first bin, because of
C     nucleation
                  q(1)=(q(1)+kcpos(1)*dt2)/(1.d0+kcneg(1)*dt2)
               ENDIF

               DO jesp=G1,G2    !Start the loop at G1 for nucleation
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                  IF (jesp.NE.ECl) THEN
#endif
C     prepare update aerosol mass per species
                     qaero(jesp)=qaero(jesp)-q(jj)
                     q(jj)=(q(jj)+kcpos(jj)*dt2)/(1.d0+kcneg(jj)*dt2)
                     q1(jj)=2.d0*q(jj)
C     update aerosol mass per species
                     qaero(jesp)=qaero(jesp)+q(jj)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                  ENDIF
#endif
               ENDDO

            ENDDO

         ENDIF

C     Compute mass conservation
         CALL MASSCNSRV(neq,nbin_aer,q,iq)

      ENDDO

      TIN2=TIN2+DT2

      END
