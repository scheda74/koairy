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

      SUBROUTINE compute_coagulation_coefficient 
     &     (Nbin_aer, bin_bound,
     &     couple, first_index, second_index,
     &     partition_coefficient)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes coagulation partition coefficients.
C     All the variables (I/O) are in common blocks.
C     
C     couple(jk): number of bin couples whose coagulation fall into bin jk.
C     (first_index(jk,jl),second_index(jk,jl)): jl-th couple of bins whose coagulation 
C     # falls into bin jk.
C     partition_coefficient(j1,j2,jk): fraction of (j1,j2) coagulation  in bin jk.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     Nbin_aer: number of aerosol bins.
C     bin_bound: bound of aerosol bins ([\mu m]).
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     couple: coagulation couples for each bin.
C     first_index: first index of coagulation couples.
C     second_index: second index of coagulation couples.
C     partition_coefficient: coagulation partition coefficient ([adim]).
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

      INTEGER Nbin_aer
      DOUBLE PRECISION bin_bound(Nbin_aer + 1)
      INTEGER couple(Nbin_aer)
      INTEGER first_Index(Nbin_aer, 4 * Nbin_aer)
      INTEGER second_index(Nbin_aer, 4 * Nbin_aer)
      DOUBLE PRECISION partition_coefficient(Nbin_aer, 
     &     Nbin_aer, Nbin_aer)

C Local variables
      DOUBLE PRECISION cst_PI6
      PARAMETER(cst_PI6=0.52359877559829887307D0)
    
      ! fake aerosol density (\mu g / m^3), needed for computation
      DOUBLE PRECISION RHOA
      PARAMETER(RHOA=1.4D12)

      INTEGER Nquad
      PARAMETER (Nquad=9)

      DOUBLE PRECISION CS(Nquad), WE(Nquad)

      DATA CS/0.1591988024618696D-01, 0.8198444633668212D-01,
     &        0.1933142836497048D+00, 0.3378732882980955D+00,
     &        0.5000000000000000D+00, 0.6621267117019045D+00,
     &        0.8066857163502952D+00, 0.9180155536633179D+00,
     &        0.9840801197538130D+00/

      DATA WE/0.4063719418078723D-01, 0.9032408034742870D-01,
     &        0.1303053482014677D+00, 0.1561735385200014D+00,
     &        0.1651196775006299D+00, 0.1561735385200014D+00,
     &        0.1303053482014677D+00, 0.9032408034742870D-01,
     &        0.4063719418078723D-01/
 
      DOUBLE PRECISION bound_mass(Nbin_aer+1)
      DOUBLE PRECISION bound_mass_log(Nbin_aer+1)
      DOUBLE PRECISION bin_width(Nbin_aer)

      INTEGER jk,j1,j2,jklo,jkhi,jl,jq
      DOUBLE PRECISION xx,fx,x12lo,x12hi
      DOUBLE PRECISION fsum,hx,xmin,xmax
      INTEGER js


C     zero initialization
      DO js=1,Nbin_aer+1
         bound_mass(js)=0.D0
         bound_mass_log(js)=0.D0
      ENDDO

      DO js=1,Nbin_aer
         bin_width(js)=0.D0
      ENDDO

      DO jk=1,Nbin_aer
         couple(jk)=0

         DO jl=1,4*Nbin_aer
            first_index(jk,jl)=0
            second_index(jk,jl)=0
         ENDDO

         DO j1=1,Nbin_aer
            DO j2=1,Nbin_aer
               partition_coefficient(j1,j2,jk)=0.D0
            END DO
         END DO
      END DO

C     fixed discretization
      DO js=1,Nbin_aer+1
                                ! bound_mass in \mu g
         bound_mass(js)=RHOA*cst_PI6*bin_bound(js)
     &        *bin_bound(js)*bin_bound(js)
         bound_mass_log(js)=DLOG(bound_mass(js))
      ENDDO

      DO js=1,Nbin_aer
         bin_width(js)=bound_mass_log(js+1)-bound_mass_log(js)
      END DO

C     Compute coagulation features: couple of bins (j1,j2) whose
C     # coagulation falls into bin jk.
      DO j2=1,Nbin_aer
         DO j1=1,j2
C     
            x12lo=DLOG(bound_mass(j1)+bound_mass(j2))
            x12hi=DLOG(bound_mass(j1+1)+bound_mass(j2+1))

                                ! find location of (j1,j2) coagulation
            CALL LOCATE(Nbin_aer+1,bound_mass_log,x12lo,jklo)
            CALL LOCATE(Nbin_aer+1,bound_mass_log,x12hi,jkhi)

            DO jk=jklo,MIN0(jkhi,Nbin_aer)
                                ! if not inside size grid then
                                ! jklo>Nbin_aer and loop not performed

                                ! count each couple
               couple(jk)=couple(jk)+1

               IF (couple(jk).GE.4*Nbin_aer) THEN
                  WRITE(6,*)"SIREAM ",
     &                 "(compute_coagulation_coefficient.f): ",
     &                 "couple(jk) >= 4*Nbin_aer",
     &                 couple(jk),jk,4*Nbin_aer
                  STOP
               ENDIF 
               
                                ! store each couple
               jl=couple(jk)
               first_index(jk,jl)=j1
               second_index(jk,jl)=j2

                                ! bounds of integration
               xmin=DMAX1(x12lo,bound_mass_log(jk))
               xmax=DMIN1(x12hi,bound_mass_log(jk+1))
               
                                ! compute partition coef
               hx=xmax-xmin
               fsum=0.D0   

               DO jq=1,Nquad
                  xx=xmin+hx*CS(jq)
                  CALL FPART(Nbin_aer,bound_mass,bound_mass_log,
     &                 bin_width,j1,j2,xx,fx)
                  fsum=fsum+fx*WE(jq)
               END DO

               partition_coefficient(j1,j2,jk)= hx*fsum
     &              /bin_width(j1)
     &              /bin_width(j2)

            END DO
         END DO
      END DO

      END
