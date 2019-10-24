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

      SUBROUTINE INITSTEP(neq,nbin_aer,q,iq,couples_coag,
     s      first_index_coag,second_index_coag,
     s      coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine performs the initialization of the timestep
C     for the integration of the GDE. The criterion is related to
C     the timescales of the aerosol processes.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NEQ : number of equations.
C     nbin_aer: number of aerosol bins.
C     Q   : gas/aerosol concentrations ([\mu.g.m^-3]).
C     iq: index of aerosol species in q(*) vector.
C     couples_coag: coagulation couples for each bin.
C     first_index_coag: first index of coagulation couples.
C     second_index_coag: second index of coagulation couples.
C     coefficient_coag: coagulation partition coefficient ([adim]).
C     XSF: neperian logarithm of fixed aerosol bin mass ([adim]).
C     MSF: fixed aerosol bin dry mass ([\mu g]).
C     DSF: fixed aerosol bin dry diameter ([\mu m]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Q   : gas/aerosol concentrations ([\mu.g.m^-3]).
C     XSD: neperian logarithm of moving aerosol bin mass ([adim]).
C     MSD: moving aerosol bin dry mass ([\mu g]).
C     DSD: moving aerosol bin dry diameter ([\mu m]).
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
      INCLUDE 'dynaero.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'time.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nbin_aer
      DOUBLE PRECISION q(neq)
      INTEGER iq(NEXT,nbin_aer)
      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      INTEGER jj
      DOUBLE PRECISION tscale,tmp,dqdt(neq)

      DOUBLE PRECISION QT(nbin_aer)

      DOUBLE PRECISION AA(NEXT,nbin_aer)

C     Time step set to output time.

      TCG=TOUT
      TCE=TOUT

C     Coagulation time step.

      IF (ICOAG.EQ.1) THEN
         ICG2=1
         ICE2=0
         INU2=0
         CALL FGDE(neq,nbin_aer,q,iq,dqdt,couples_coag,
     s      first_index_coag,second_index_coag,
     s      coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA)

         DO jj=1,neq
            tmp=q(jj)*dqdt(jj)
            IF (tmp.NE.0.D0) THEN
               tscale=q(jj)/DABS(dqdt(jj))
               TCG=DMIN1(TCG,tscale)
            ENDIF
         END DO
      ENDIF

C     Cond/evap time step.

      IF (ICOND.EQ.1) THEN
         ICG2=0
         ICE2=1
         INU2=INUCL
         ICUT2=ICUT             ! current cutting size
         CALL FGDE(neq,nbin_aer,q,iq,dqdt,couples_coag,
     s      first_index_coag,second_index_coag,
     s      coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA)

         DO jj=1,neq
            tmp=q(jj)*dqdt(jj)
            IF (tmp.NE.0.D0) THEN

               tscale=q(jj)/DABS(dqdt(jj))
               TCE=DMIN1(TCE,tscale)
            ENDIF
         END DO
      ENDIF

C     Minimal time step.
                                ! it is assumed that coagulation is
                                ! always the slowest process
      DT=TCG
      DT=DMIN1(DT,TOUT-TIN)


                                ! if only one process then set
                                ! splitting step to whole timestep

      DT = DMAX1(DT,DTAEROMIN)
      TCE = DMIN1(DT,DMAX1(TCE,DTAEROMIN))
      TCG = DMIN1(DT,DMAX1(TCG,DTAEROMIN))
      IF(ICOAG.EQ.0.OR.ICOND.EQ.0)  DT=TOUT-TIN
      IF(ICUT.EQ.nbin_aer) DT=TOUT-TIN

      END

