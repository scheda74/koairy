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

      SUBROUTINE PROCESSAERO(neq,nesp_aer,nbin_aer,q,iq,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD, bin_density)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine integrates the different processes for the GDE.
C
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
C     2005/9/30: Add possibility to compute sulfate cond. explicitely
C     + Remove full eq. solver from FGDE (Karine Sartelet, CEREA).
C     2013/11/27: Added bin_density (Stephanie Deschamps, CEREA).
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
      INCLUDE 'time.inc'
      INCLUDE 'varp.inc'

      INTEGER neq,nesp_aer,nbin_aer
      DOUBLE PRECISION q(neq),q1(neq)
      INTEGER iq(nesp_aer,nbin_aer)

      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer),DSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer),DSD(nbin_aer)

      DOUBLE PRECISION QT(nbin_aer)
      double precision bin_density(nbin_aer)

      integer js,jesp

C     ******time loop
      DO WHILE ( TIN2 .LT. TOUT2 )

         IF (ICE2.EQ.1) THEN
C     Compute gas mass conservation.
            CALL MASSCNSRV(neq,nesp_aer,nbin_aer,q,iq)

C     solve sulfate dynamically
            IF (ISULFCOND.EQ.1) CALL SULFDYN(neq,nesp_aer,nbin_aer,
     &           q,iq,QT,XSF,MSF,DSF,XSD,MSD,DSD, bin_density)


C     solve equilibrium with
C     so-called bulk approach
            IF (ICUT2.GT.0) CALL BULKEQUI(neq,nesp_aer,nbin_aer,q,iq,QT,
     s           XSF,MSF,DSF,XSD,MSD,DSD,bin_density)

         ENDIF

         IF(ICUT2.EQ.nbin_aer.AND.
     s        INU2.EQ.0.AND.
     s        ICG2.EQ.0) THEN
            DT2 = TOUT2 - TIN2
            TIN2 = TOUT2
         ElSEIF (KDSLV2.EQ.1) THEN
                                ! etr dynamic solver

            CALL ETRCONC(neq,nesp_aer,nbin_aer,q1,q,iq,couples_coag,
     s           first_index_coag,second_index_coag,
     s           coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,
     s           bin_density)

         ELSEIF (KDSLV2.EQ.2) THEN
                                ! ros2 dynamic solver

            CALL ROS2CONC(neq,nesp_aer,nbin_aer,q1,q,iq,couples_coag,
     s           first_index_coag,second_index_coag,
     s           coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,
     s           bin_density)

         ELSEIF (KDSLV2.EQ.3) THEN
                                ! ebi dynamic solver

            CALL EBICONC(neq,nesp_aer,nbin_aer,q1,q,iq,QT,
     s           XSF,MSF,DSF,XSD,MSD,DSD)
         ENDIF

                                ! adaptive time step
C     check the need to call ADAPTSTEP (not for the last timestep)
         IF (TIN2.LT.TOUT2 ) THEN
            CALL ADAPTSTEP(neq,nesp_aer,nbin_aer,q1,q,iq)
         ENDIF

      END DO

      END
