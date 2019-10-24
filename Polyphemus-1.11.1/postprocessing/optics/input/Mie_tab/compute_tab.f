      PROGRAM MAIN

      IMPLICIT NONE


      integer i,j,k
      integer Nre, Nimag, Ndiam
      PARAMETER(Nre=88,Nimag=100,Ndiam=100)
      double precision index_re_class(Nre+1),
     $     index_imag_class(Nimag+1)
      double precision diameters(Ndiam+1)

      double precision pi
      PARAMETER (pi=3.1416)

c     Change lambda for an other wavelength.
      double precision lambda0
      integer lambda0_nm
      PARAMETER (lambda0=5.50D-1)
      character*3 filen

      real xx
      complex ind_cplx

      logical perfct, anyang, prnt(2)

      real qext, qsca, gqsc, xmu(7), angle(7)
      real spike
      complex tforw(2), tback(2), s1(7), s2(7), sforw, sback
      real mimcut

C     Variables used in the self-test of Mie computations.
      real pmom(0:2, 4)
      integer numang, nmom, ipolzn, momdim

      real Qsca_save(Nre+1,Nimag+1,Ndiam+1),
     $     Qext_save(Nre+1,Nimag+1,Ndiam+1)

      integer NUC, NUC_grid

      real hx

      NUC_grid = 3
      NUC = 2
      lambda0_nm = lambda0*1.D3
      write(filen,'(I3)') lambda0_nm

      OPEN(NUC,FILE='efficiency_factors_tab_' // filen // '.dat')
      OPEN(NUC_grid,FILE='Grid_Mie.dat')

      perfct = .false.
      anyang = .false.
      mimcut = 1.D-9
      PRNT( 1 ) = .FALSE.
      PRNT( 2 ) = .FALSE.
      DO i = 1, 7
         ANGLE( i ) = ( i - 1 )*180. / ( 7 - 1 )
         xmu( i ) = COS( ANGLE( I ) * pi / 180. )
      enddo

C     Tabulation:
      do i=1,Nre+1
         index_re_class(i) = 1.11D0 + (i-1)*1D-2
      enddo
      do i=1,Nimag+1
         index_imag_class(i) = (i-1) * 4.4D-3
      enddo

      hx = log(20/1D-2) / Ndiam
      do i=1,Ndiam+1
         diameters(i) = 1D-2 * exp((i-1)*hx)
      enddo

      WRITE(NUC_grid,*) '#index_re_class= '
      WRITE(NUC_grid,*) (index_re_class(i),i=1,Nre+1)
      WRITE(NUC_grid,*) '#index_imag_class= '
      WRITE(NUC_grid,*) (index_imag_class(i),i=1,Nimag+1)
      WRITE(NUC_grid,*) '#diameters= '
      WRITE(NUC_grid,*) (diameters(i),i=1,Ndiam+1)

      xx = REAL( 2 * pi * 1.0D-2 / lambda0 )

      do i=1,Nre+1
         do j=1,Nimag+1
            do k=1,Ndiam+1
               xx = REAL( 2 * pi * diameters(k) / lambda0 )
               ind_cplx = CMPLX(index_re_class(i), index_imag_class(j))

C     Compute Qabs et Qext

C     Variables used in self-test of Mie computation.
               numang = 0
               nmom = 0
               ipolzn = 0
               momdim = 2

               call MIEV0(xx, ind_cplx, perfct, mimcut, anyang, numang,
     $              xmu, nmom, ipolzn, momdim, PRNT, Qext, Qsca, gqsc,
     $              pmom, sforw, sback, s1, s2, tforw, tback, spike)
               Qsca_save(i,j,k)=Qsca
               Qext_save(i,j,k)=Qext
            enddo
         enddo
      enddo
      WRITE(NUC,*) '# Qsca'
      WRITE(NUC,*) (((Qsca_save(i,j,k),k=1,Ndiam+1),j=1,Nimag+1)
     $     ,i=1,Nre+1)
      WRITE(NUC,*) '# Qext'
      WRITE(NUC,*) (((Qext_save(i,j,k),k=1,Ndiam+1),j=1,Nimag+1)
     $     ,i=1,Nre+1)

      END
