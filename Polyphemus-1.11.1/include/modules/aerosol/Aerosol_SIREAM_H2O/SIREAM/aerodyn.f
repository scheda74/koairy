C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
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

      SUBROUTINE AERODYN(NGAS,nesp_aer,nbin_aer,
     s     T,P,HUMID,T0,T1,ZA,couples_coag,
     s     first_index_coag,second_index_coag,
     s     coefficient_coag,XSF,MSF,DSF,XBF,MBF,DBF,HSF,iq,ZNA,
     s     section_pass, DQLIMIT)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This is the main subroutine for solving the GDE for aerosol
C     with a size-resolved method in a given grid cell.
C     The initialization is performed for the aerosol variables, the GDE
C     is solved and the variables are modified for the CTM.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NGAS: number of gas species.
C     nesp_aer: number of aerosol species.
C     nbin_aer: number of aerosol bins.
C     T: temperature ([K]).
C     P: pressure ([Pa]).
C     HUMID: specific humidity ([%]).
C     T0: input timestep ([s]).
C     T1: output timestep ([s]).
C     ZA: vector of gas ans aerosol species ([\mu g.m^-3]).
C     couples_coag: coagulation couples for each bin.
C     first_index_coag: first index of coagulation couples.
C     second_index_coag: second index of coagulation couples.
C     coefficient_coag: coagulation partition coefficient ([adim]).
C     XSF: neperian logarithm of fixed aerosol bin mass ([adim]).
C     MSF: fixed aerosol bin dry mass ([\mu g]).
C     DSF: fixed aerosol bin dry diameter ([\mu m]).
C     XBF: neperian logarithm fixed aerosol bound dry mass ([]).
C     MBF: fixed aerosol bound dry mass ([\mu g]).
C     DBF: fixed aerosol bound dry diameter [\mu m]).
C     HSF: fixed width of aerosol bins ([adim]).
C     iq: index of aerosol species in q(*) vector.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     ZA: vector of gas and aerosol concentrations ([\mu g.m^-3]).
C     ZNA: vector of aerosol number concentrations ([m^-3]).
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
C     2004: add aqueous-phase module (Kathleen Fahey, CEREA).
C     2004: add clipping test at the end (Kathleen Fahey, CEREA).
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Kathleen Fahey, CEREA, 2004.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'param.inc'
      INCLUDE 'dynaero.inc'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'time.inc'
      INCLUDE 'meteo.inc'
      INCLUDE 'varq.inc'
      INCLUDE 'varp.inc'
      INCLUDE 'varg.inc'
      INCLUDE 'vara.inc'

      INTEGER NGAS,nesp_aer,nbin_aer
      INTEGER i,js,jesp,j,jk, ii
      INTEGER neq,NESP
      DOUBLE PRECISION q(nbin_aer*(1+nesp_aer)+nesp_aer)
      DOUBLE PRECISION T0,T1,T,P,HUMID
      DOUBLE PRECISION ZA(NGAS+nbin_aer*nesp_aer)
      DOUBLE PRECISION ZNA(nbin_aer)
      DOUBLE PRECISION totalwat,VIM
      DOUBLE PRECISION emw_tmp
      INTEGER couples_coag(nbin_aer)
      INTEGER first_index_coag(nbin_aer, 4 * nbin_aer)
      INTEGER second_index_coag(nbin_aer, 4 * nbin_aer)
      DOUBLE PRECISION coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

      INTEGER iq(nesp_aer,nbin_aer)

      DOUBLE PRECISION XSF(nbin_aer),MSF(nbin_aer)
      DOUBLE PRECISION DSF(nbin_aer),XBF(nbin_aer+1)
      DOUBLE PRECISION MBF(nbin_aer+1),DBF(nbin_aer+1)
      DOUBLE PRECISION HSF(nbin_aer)
      DOUBLE PRECISION XSD(nbin_aer),MSD(nbin_aer)
      DOUBLE PRECISION DSD(nbin_aer),XBD(nbin_aer+1)
      DOUBLE PRECISION MBD(nbin_aer+1),DBD(nbin_aer+1)
      DOUBLE PRECISION HSD(nbin_aer)

      DOUBLE PRECISION bin_density(nbin_aer) !! bin mass density = rho

      DOUBLE PRECISION QT(nbin_aer)

      DOUBLE PRECISION dsum,kcorr,tmp
      DOUBLE PRECISION dqdt(nbin_aer*(1+nesp_aer)+nesp_aer)
      DOUBLE PRECISION AA(nesp_aer,nbin_aer)
      DOUBLE PRECISION conc(nbin_aer, nesp_aer)
      DOUBLE PRECISION numconc(nbin_aer)
      DOUBLE PRECISION fixed_diameter(nbin_aer), dnew(nbin_aer)
      INTEGER section_pass
      DOUBLE PRECISION DQLIMIT

C     Set local dimensions.
      NESP = NGAS+nbin_aer*nesp_aer
      neq = nbin_aer*(1+nesp_aer)+nesp_aer

C     Set discretizations.
      DO js=1,nbin_aer+1
         MBD(js)=MBF(js)
         DBD(js)=DBF(js)
         XBD(js)=XBF(js)
      ENDDO

      DO js=1,nbin_aer
         MSD(js)=MSF(js)
         DSD(js)=DSF(js)
         XSD(js)=XSF(js)
         HSD(js)=HSF(js)
      ENDDO

C     Update numerical parameters.
      TIN = 0.D0
      TOUT = T1-T0

c
      DT=0.D0
      TIN2=0.D0
      TOUT2=0.D0
      DT2=0.D0
      TCG=0.D0
      TCE=0.D0
c
      KDSLV2=0
      IKV2=0
      ICUT2=0
c
      ICG2=0
      ICE2=0
      INU2=0

C     *************************************************

C     Update thermo and microphysical parameters.

      TEMP = T
      PRES = P

      CALL COMPUTE_RELATIVE_HUMIDITY(HUMID,TEMP,PRES,RH)
      RH = DMIN1(DMAX1(RH, Threshold_RH_inf), Threshold_RH_sup)

!     air free mean path
      IF (ICOAG.EQ.1) THEN
         CALL COMPUTE_AIR_FREE_MEAN_PATH(TEMP,
     &        PRES,AIRFMP,VIM)

      ENDIF

      IF (ICOND.EQ.1) THEN
         DO i=1,NEXT
            DIFFG(i)=0.D0
            VQMG(i)=0.D0
            QSAT(i)=0.D0
            KPART(i)=0.D0
         ENDDO

         DO i=1,nesp_aer
            IF (aerosol_species_interact(i).GT.0) THEN
               emw_tmp = EMW(i) * 1.D-6 ! g/mol
               CALL COMPUTE_GAS_DIFFUSIVITY(TEMP,PRES,
     &              SIGM(i),emw_tmp,PARM(i),DIFFG(i) ) ! gas diff coef in air

               CALL COMPUTE_QUADRATIC_MEAN_VELOCITY(TEMP,
     $              emw_tmp, VQMG(i) ) ! gas quad mean speed in air
            ENDIF
         ENDDO

         DO i=1,nesp_pankow
            j=pankow_species(i)
            emw_tmp = EMW(j) * 1.D-6 ! g/mol
            CALL COMPUTE_SATURATION_CONCENTRATION(TEMP,
     $           emw_tmp, DHVAP(j), PSATREF(j), QSAT(j) )
         ENDDO

         DO i = 1, nesp_pom
            j = poa_species(i)
            emw_tmp = EMW(j) * 1.D-6 ! g/mol
            CALL COMPUTE_SATURATION_CONCENTRATION(TEMP,
     $           emw_tmp, DHVAP(j), PSATREF(j), QSAT(j))
         ENDDO

         tmp = 1.D0/RGAS * (1.D0 / TEMP - 1.D0 / 298.D0)
!     the formula is inversed compared
!     to that of vapore pressure because
!     partition coefficient are inversely
!     proportional to vapore pressure
         DO i=1,nesp_pankow
            j=pankow_species(i)
            KPART(j) = TEMP / 298.D0 ! temperature dependency of partition coefficient
     &           * KPARTREF(j) * DEXP(DHVAP(j) * tmp)
         ENDDO
!     Reference at 300K for prymary organic aerosol.
         tmp = 1.D0 / RGAS * (1.D0 / TEMP - 1.D0 / 300.D0)
         DO i = 1, nesp_pom
            j = poa_species(i)
!     Temperature dependency of partition coefficient.
            KPART(j) = TEMP / 300.D0
     &           * KPARTREF(j) * DEXP(DHVAP(j) * tmp)
         ENDDO
      ENDIF

C     **************************************************
C     INITIALIZE CURRENT VARIABLES
C     **************************************************

      DO i = 1,neq
         q(i) = 0.D0
      ENDDO

C     Init gas/aerosol mass and number density
C     semi_volatile gas-phase concentrations

      DO jesp = 1,nesp_aer
         if (aerosol_species_interact(jesp).GT.0)
     $        q(IG(jesp)) = ZA(aerosol_species_interact(jesp))
      ENDDO

C     external dry composition (notice that the aerosol module
C     does not need H20 concentration)
C     + volumic sources

      DO js = 1,nbin_aer
         DO jesp = E1,E2
            q(IQ(jesp,js)) = ZA(NGAS+(jesp-1)*nbin_aer+js)
            conc(js, jesp) = q(IQ(jesp,js))
         ENDDO
         numconc(js) = q(js)
      ENDDO

C     Given mass volumic emission, compute volumic emission for number

      DO js=1,nbin_aer
         QT(js) = 0.D0

C     compute the total dry mass
         DO jesp = E1,E2
            QT(js) = QT(js)+q(IQ(jesp,js))
         ENDDO

C     compute the particle number (mass/geometric mean diameter)
C     Number concentration 

         bin_density(js) = RHOA
         IF (IDENS.EQ.1) THEN
           CALL compute_density(nbin_aer,nesp_aer,EH2O, TINYM,conc,
     &                            LMD,js,bin_density(js))
         ENDIF
         MSD(js)= cst_pi6*bin_density(js)*DSD(js)**3.D0
         XSD(js)=DLOG(MSD(js))


         IF(INUM.EQ.1) THEN
            q(js)=ZNA(js)       ! real value
         ELSE
            q(js)=QT(js)/MSD(js) ! approximate value 
         ENDIF
 
      ENDDO

C     compute total gas/aerosol mass in cell.
      DO jesp=1,nesp_aer
         QTOT(jesp)=0.D0
      ENDDO

      DO jesp=E1,E2
         DO js=1,nbin_aer
            j=IQ(jesp,js)
            QTOT(jesp)=QTOT(jesp)+q(j)
         END DO
      END DO

      DO jesp=E1,E2
         IF (aerosol_species_interact(jesp).GT.0) THEN
            j=IG(jesp)
            QTOT(jesp)=QTOT(jesp)+q(j)
         ENDIF
      END DO


      DO jk = 1, Nbin_aer
         fixed_diameter(jk) = DSQRT(DBF(jk) * DBF(jk+1))
      ENDDO


C     **************************************************
C     SOLVE GAS/AEROSOL GDE LAGRANGIAN EQS
C     **************************************************

      DO WHILE (TIN.LT.TOUT)
C     Compute timestep for each process
         CALL INITSTEP(neq,nesp_aer,nbin_aer,q,iq,couples_coag,
     s        first_index_coag,second_index_coag,
     s        coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,
     s        bin_density)

         conc = 0.d0
         DO js = 1,nbin_aer 
            DO jesp = E1,E2
               conc(js, jesp) = q(IQ(jesp,js))
            ENDDO
            numconc(js) = q(js)
         ENDDO
         
         DO jk=1, nbin_aer
            IF (IDENS.EQ.1) THEN
               CALL compute_density(nbin_aer,nesp_aer,EH2O, TINYM,conc,
     &              LMD,jk,bin_density(jk))
            ELSE
               bin_density(jk) = RHOA
            ENDIF            
         ENDDO

C     Begin with the slowest process (coagulation).
         IF (ICOAG.EQ.1) THEN

            TIN2=TIN
            TOUT2=TIN2+DT
            ICG2=ICOAG
            ICE2=0
            INU2=0
            DT2=TCG
            ICUT2 = 0
            KDSLV2=1            ! only etr for coagulation

            CALL PROCESSAERO(neq,nesp_aer,nbin_aer,q,iq,couples_coag,
     s           first_index_coag,second_index_coag,
     s           coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,
     s           bin_density)

c     Check mass conservation in order not to impact gas-phase
c     # through mass conservation for C/E

            do jesp=E1,E2
               dsum=0.D0
               do js=1,nbin_aer
                  dsum=dsum+q(IQ(jesp,js))
               enddo
               kcorr=1.D0
               if (dsum.gt.0.D0) kcorr=(qtot(jesp)-q(IG(jesp)))/dsum
               do js=1,nbin_aer
                  q(IQ(jesp,js))=q(IQ(jesp,js))*kcorr
               ENDDO
            ENDDO

            DO js = 1,nbin_aer
               QT(js) = 0.d0 
               DO jesp = E1,E2
                  QT(js) = QT(js)+q(IQ(jesp,js)) 
               ENDDO
            ENDDO

         ENDIF

C     Continue with the fastest process (condensation+nucleation).
         IF (ICOND+INUCL.GT.0) THEN

            TIN2   = TIN
            TOUT2  = TIN2 + DT
            ICG2   = 0
            ICE2   = ICOND
            INU2   = INUCL
            DT2    = TCE
            KDSLV2 = KDSLV
            ICUT2  = ICUT

            CALL PROCESSAERO(neq,nesp_aer,nbin_aer,q,iq,couples_coag,
     s           first_index_coag,second_index_coag,
     s           coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,
     s           bin_density)

C     Update all variables after last c/e step
            CALL FGDE(neq,nesp_aer,nbin_aer,q,iq,dqdt,couples_coag,
     s           first_index_coag,second_index_coag,
     s           coefficient_coag,QT,XSF,MSF,DSF,XSD,MSD,DSD,AA,
     s           bin_density)

C     If necessary, set to zero tiny concentrations
C     To do similarly as in aerosol.f

            IF((IREDIST.EQ.1).OR.(IREDIST.EQ.2)) THEN
C     Compute size bounds
            CALL SIZEBND(nbin_aer,XSF,XBF,XSD,XBD,MBD,DBD,HSD)

C     Redistribution onto the fixed grid
            
! array of concentration at the same length of the code
                          
               IF (IREDIST.EQ.1) CALL REDIST_NEW(neq,nesp_aer,
     s           nbin_aer,q,iq,QT,MSD,MSF,DSF,XBF,
     s           MBF,DBF,HSF,XBD,MBD,DBD,HSD)
            
                        
               IF (IREDIST.EQ.2) CALL REDIST(neq,nesp_aer,nbin_aer,q,
     s           iq,QT,XBF,MBF,DBF,HSF,XBD,MBD,DBD,HSD)

            ENDIF

C     Eulerian redistributions

            IF (IREDIST.EQ.3 .OR. IREDIST.EQ.4  
     s           .OR. IREDIST.EQ.5 .OR. IREDIST.EQ.6
     s           .OR. IREDIST.EQ.7 .OR. IREDIST.EQ.8
     s           .OR. IREDIST.EQ.9 .OR. IREDIST.EQ.10) THEN

! Concentration vectors for redist_euler :
               conc = 0.d0
               numconc = 0.d0
               DO js=1,nbin_aer
                  DO jesp= E1, E2
                     conc(js,jesp) = q(IQ(jesp,js))
                  ENDDO
                  numconc(js)=q(js)
               ENDDO
               
               CALL REDISTRIBUTION(nbin_aer,nesp_aer,EH2O,DBF,
     &              fixed_diameter,RHOA,IDENS,IREDIST,section_pass,LMD,
     &              DQLIMIT,conc,numconc,QT)

               DO js=1,nbin_aer
                  DO jesp= E1, E2
                     q(IQ(jesp,js)) = conc(js, jesp)
                  ENDDO
                  q(js)=numconc(js)
               ENDDO

            ENDIF

         ENDIF

         TIN=TIN+DT
      END DO

C     **************************************************
C     ******print all ouptut files

C     update variables for dispersion equation

      DO jesp = 1,nesp_aer
         if (aerosol_species_interact(jesp).GT.0)
     $        ZA(aerosol_species_interact(jesp)) = q(IG(jesp))
      ENDDO

      DO js = 1,nbin_aer
         DO jesp = 1,nesp_aer   ! all species including water
            ZA(NGAS+(jesp-1)*nbin_aer+js)=q(IQ(jesp,js))
         ENDDO
      ENDDO

      IF(INUM.EQ.1) THEN
         DO js = 1,nbin_aer
            ZNA(js) = q(js)
         ENDDO
      ENDIF
      totalwat = 0.d0
      DO js =1,nbin_aer
         totalwat = totalwat + q(IQ(EH2O,js))
      ENDDO

      IF (totalwat .GE. 1.D06)
     &     write(6,*) totalwat,'SIREAM (aerodyn.f): total water>1e6'

      DO i = 1, NESP
         if(ZA(i) .lt. 0.d0) ZA(i) = 0.d0
      ENDDO

      END
