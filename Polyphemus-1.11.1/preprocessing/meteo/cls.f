C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Luc Musson-Genon
C
C     This file is part of the air quality modeling system Polyphemus.
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


C This file implements parameterizations to compute several quantities
C in the surface layer, such as the friction wind or the Monin-Obukhov
C length.

      subroutine cls(u2,v2,t2,t1,q2,q1,imer,zh,zb,zbt,z,lmo,ustar,
     $     qo, eo, fu, ft)
c
      complex uv
      real lmo,lmo1,k,fu,ft,cch
      logical stab
      integer cpte
      integer imer
      complex rot,fuvzo,duvdz
c
      common/sort/uz,vz,ffz,ddz,tetaz,qz
      common/grad/dudz,dvdz,dudzcls,dodz,dqdz
c
      common/phi/phiu,phit
      common/cbul/cka
c
c     ce programme calcule les flux uaster,qo,eo a partir des gradients
c     verticaux de u,t,q entre le sol et le dernier niveau du modele
c
c constantes physiques
c
c gravitation
       g=9.81
c karman
      cka=0.4
c hauteur de deplacement
      zd=0.
c inverse Prandtl turbulent
      at=1.
c constante de Charnock
      cch=0.0144
c
      difu=sqrt(u2*u2+v2*v2)
      difu=amax1(0.5,difu)
      dift=t2-t1
      difq=q2-q1
      tob=(t2+t1)/2.
      dtq=dift+0.608*difq*tob
      zd=0.
c
      if(dtq.gt.0.) then
c cas stable
      zth=(zh+zbt-zd+zbt)/2.
      dt=dift/(zh+zbt-zd-zbt)
      dq=difq/(zh+zbt-zd-zbt)
      dth=dtq/(zh+zbt-zd-zbt)
      lmo1=50.
      stab=.TRUE.
      else
c cas instable
      zth=sqrt((zh+zbt-zd)*zbt)
      dt=dift/(zth*alog((zh+zbt-zd)/zbt))
      dq=difq/(zth*alog((zh+zbt-zd)/zbt))
      dth=dtq/(zth*alog((zh+zbt-zd)/zbt))
      lmo1=-50.
      stab=.FALSE.
      endif
c
      du=difu
      tm=tob
      k=cka
      epsilon= 0.0001
      epsilon2=0.01
C      if(abs(dtq).gt.1.e-8) then
      if(dtq.ne.0) then
c         *************************************
c         BOUCLE METHODE ITERATIVE
      cpte=0.
c
c       print *,'cpte=',cpte
       ustar1=fustar(du,lmo1,zh,zb,cka)
c       print *,'ustar1=',ustar1
      do  icompt=1,15
       if(imer.eq.1) then
       zb=cch*ustar1*ustar1/g
       zbt=zb/10.
      if(dtq.gt.0.) then
c cas stable
      zth=(zh+zbt-zd+zbt)/2.
      dt=dift/(zh+zbt-zd-zbt)
      dq=difq/(zh+zbt-zd-zbt)
      dth=dtq/(zh+zbt-zd-zbt)
      else
c cas instable
      zth=sqrt((zh+zbt-zd)*zbt)
      dt=dift/(zth*alog((zh+zbt-zd)/zbt))
      dq=difq/(zth*alog((zh+zbt-zd)/zbt))
      dth=dtq/(zth*alog((zh+zbt-zd)/zbt))
      endif
      endif

200    lmo=flmo(tm,ustar1,dth,zth,stab,g,cka)
       ustar=fustar(du,lmo,zh,zb,cka)
C       print *,'ustar=',ustar,'lmo=',lmo
         ustar1=ustar
         lmo1=lmo
         if(ustar.le.5.e-3) go to 202
      enddo
 202  continue
      call  foncuo(ustar,lmo,zth,zb,zbt,zd,cka,fu,ft)
C      print *,'phit=',phit
      qo=-ustar*dt*cka*zth/phit
      eo=-ustar*dq*cka*zth/phit
      almo=-(ustar)**3./(cka*g/tm*(qo+0.608*tm*eo))
C      almo=lmo
      uaster=ustar
      uaster2=ustar*ustar
      taster=-qo/uaster
      qaster=-eo/uaster
      else
      qo=0.
      eo=0.
      lmo=1.e-30
      almo=lmo
      ustar=du*cka/log((zh+zb-zd)/zb)
      uaster=ustar
      uaster2=ustar*ustar
      taster=-qo/uaster
      qaster=-eo/uaster
      ft=0.
      fu=log((zh+zb-zd)/zb)
      endif
C      print 101,z,uaster,qo,eo
101   format(10x,'z=',f3.0,2x,'ustar=',e10.3,2x,'qo=',e10.3,2x,'eo='
     s,e10.3)
c
c***** calcul des flux de mvt a zo dans le repere meteo
c
      al=atan(v2/u2)
      if(u2.lt.0.) al=al+acos(-1.)
      rot=cexp((0.,1.)*al)
      fuvzo=-uaster2*rot
      fuzo=real(fuvzo)
      fvzo=aimag(fuvzo)
c
c     calcul du vent et de la temperature a zm

c
      call  foncuo(uaster,lmo,z,zb,zbt,zd,cka,fu,ft)
C      print *,'fu=',fu,'ft=',ft
      ucls=uaster*fu/cka
      uv=ucls*rot
      uz=real(uv)
      vz=aimag(uv)
      tetaz=t1+taster*ft/cka
      qz=q1+qaster*ft/cka
c
c     calcul du vent a z m en dd ff
c
      ffz=ucls
      ddz=al*180./(4.*atan(1.))
      ddz=270.-ddz
      if(ddz.gt.360.)ddz=ddz-360.
c
c***** calcul des gradients au niveau zh
c
      dodz=taster*phit/(cka*z)
      dqdz=qaster*phit/(cka*z)
      dudzcls=uaster*phiu/(cka*z)
      duvdz=dudzcls*rot
      dudz=real(duvdz)
      dvdz=aimag(duvdz)
C      print *,'dodz=',dodz,'dqdz=',dqdz,'dudz=',dudz,'dvdz=',dvdz
c
      return
      end
c
c
c ************************************************************
c
c fonction calculant U ( uaster , zb, zh, Lmo,zd)
c utilisant les fonctions psi et psi0 decrite plus haut
c
       subroutine foncuo(uaster,lmo,zh,zb,zbt,zd,k,fu,ft)
      implicit none
      external psiksi
      real psiksi,psi,psi0,du,zb,zh,lmo,k,ksib,ksih
      real t1,t2,t3,t4,uaster,phiu,phit,fu,ft,zbt
      real ksitb,ksith,zd

      common/phi/phiu,phit
c
      ksib=zb/lmo
      ksitb=zbt/lmo
      ksih=(zh-zd+zb)/lmo
      ksith=(zh-zd+zbt)/lmo
      if (ksih.le.(0.0)) then
        psi=psiksi(ksih)
        psi0=psiksi(ksib)
      endif

      if (ksih.le.0.0) then
        t1=log((zh-zd+zb)/zb)
        t2=-2.0*log((1+psi)/(1+psi0))
        t3=-log( (1+ (psi)**2 ) / (1+ (psi0)**2 ) )
        t4=2.0*(atan(psi)-atan(psi0))
        fu=t1+t2+t3+t4
        ft=0.74*(log((zh-zd+zbt)/zbt)-2.*(log(1.+(1.-9.*ksith)**(1./2.)
     s  )/2.)-log((1.+(1.-9.*ksith)**(1/2.))/2.))
        phiu=(1.-15.*ksih)**(-1./4.)
        phit=0.74*(1.-9.*ksith)**(-1./2.)
        goto 9000
       else
        if (ksih.lt.(0.5)) then
          t1=log((zh-zd+zb)/zb)
          t2=5.0*(ksih-ksib)
          fu=t1+t2
          ft=0.74*log((zh-zd+zbt)/zbt)+4.7*(ksith-ksitb)
          phiu=1.+5.*ksih
          phit=0.74+4.7*ksith
          goto 9000
         else
           if (ksih.lt.(10.0)) then
             t1=8.0*log(2.0*ksih)
             t2=4.25/ksih-0.5/( ksih**2 )
             t3=-log(2.0*ksib)-5.0*ksib-4.0
             fu=t1+t2+t3
             ft=0.74*log((zh-zd+zbt)/zbt)+4.7*(ksith-ksitb)
             phiu=8.-4.25/ksih+1./(ksih*ksih)
             phit=0.74+4.7*ksith
             goto 9000
            else
             t1=0.7585*ksih
             t2=8.0*log(20.0)-11.165
             t3=-log(2.0*ksib)-5.0*ksib
             fu=t1+t2+t3
             ft=0.74*log((zh-zd+zbt)/zbt)+4.7*(ksith-ksitb)
             phiu=0.7585*ksih
             phit=0.74+4.7*ksith
             goto 9000
           endif
        endif
      endif

 9000 continue
      return
      end
c
c
c ************************************************************
c
c fonction calculant uaster  (du, zb, zh, Lmo)
c utilisant les fonctions psi et psi0 decrite plus haut
c
      real function fustar(du,lmo,zh,zb,k)
      implicit none
      external psiksi
      real psiksi,psi,psi0,du,zb,zh,lmo,k,ksib,ksih
      real t1,t2,t3,t4
c
c
      ksib=zb/lmo
      ksih=(zh+zb)/lmo
c      print *,du,lmo,zh,zb,k
      if (ksih.le.(0.0)) then
        psi=psiksi(ksih)
        psi0=psiksi(ksib)
      endif

      if (ksih.le.0.0) then
        t1=log(zh/zb)
        t2=-2.0*log((1+psi)/(1+psi0))
        t3=-log( (1+ (psi)**2 ) / (1+ (psi0)**2 ) )
        t4=2.0*(atan(psi)-atan(psi0))
        fustar=(k*du)/(t1+t2+t3+t4)
        goto 9000
       else
        if (ksih.lt.(0.5)) then
          t1=log(zh/zb)
          t2=5.0*(ksih-ksib)
          fustar=(k*du)/(t1+t2)
          goto 9000
         else
           if (ksih.lt.(10.0)) then
             t1=8.0*log(2.0*ksih)
             t2=4.25/ksih-0.5/( ksih**2 )
             t3=-log(2.0*ksib)-5.0*ksib-4.0
             fustar=(k*du)/(t1+t2+t3)
             goto 9000
            else
             t1=0.7585*ksih
             t2=8.0*log(20.0)-11.165
             t3=-log(2.0*ksib)-5.0*ksib
             fustar=(k*du)/(t1+t2+t3)
             goto 9000
           endif
        endif
      endif

 9000 continue
      return
      end
c
c ************************************************************
c
      real function flmo(thm,ustar,dth,zth,stab,g,k)
      implicit none
      real thm,ustar,dth,zth
      logical stab
      real alpha,r,t1,t2
      real g,k
c
      alpha=(thm*ustar*ustar*0.74)/(g*k*k*dth*zth)
c      print *,thm,ustar,dth,zth,stab,g,k,alpha
c
      if (.not.stab) then
          t1=(9.0*zth)/(2.0*abs(alpha))
          t2=sqrt(1+ (t1*t1) )
          r=-(t1+t2)/abs(alpha)
          flmo=1/r
        else
          t1=0.74/(2.0*4.7*zth)
          t2=0.74/(alpha*4.7*zth)
          r=-t1+sqrt(t2+t1*t1)
          flmo=1/r
      endif
c      print *,t1,t2,r
c
      return
      end
c
c FONCTIONS ************************************************** FONCTIONS
c
c fonction calculant psi(z/Lmo)
c
      real function psiksi(ksi)
      implicit none
      real ksi,gm1
      parameter(gm1=15.0)
c
      psiksi=(1-gm1*ksi)**0.25
      return
      end
