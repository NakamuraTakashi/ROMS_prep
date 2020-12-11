      subroutine ll2utm(lat,lon,x,y,iizone,ispher)
c
c     universal transverse mercator conversion
c
c     rewritten 6-16-83 by j.f. waananen using formulation
c     by j.p. snyder in bulletin 1532, pages 68-69
c
c      rewrite utm.f for the case idir=0, which convert lat lon to utm
c       2/5/2001 j. klinck
c
      implicit none
c
      real*8 axis(19),bxis(19)
      real*8 lat,lon
      real*8 radsec
      real*8 ak0,a,b,es
      real*8 x,y,slat,slon,cm,phi,dlam,epri,en,t
      real*8 c,aa,s2,s4,s6,f1,f2,f3,f4,em,xx,yy,secs
c
      integer iizone,ispher,izone
      integer iiz,iutz
c
      data secs/3600.0d0/
c
      data axis/6378206.4d0,6378249.145d0,6377397.155d0,
     . 6378157.5d0,6378388.d0,6378135.d0,6377276.3452d0,
     . 6378145.d0,6378137.d0,6377563.396d0,6377304.063d0,
     . 6377341.89d0,6376896.0d0,6378155.0d0,6378160.d0,
     . 6378245.d0,6378270.d0,6378166.d0,6378150.d0/
c
      data bxis/6356583.8d0,6356514.86955d0,6356078.96284d0,
     . 6356772.2d0,6356911.94613d0,6356750.519915d0,6356075.4133d0,
     . 6356759.769356d0,6356752.31414d0,6356256.91d0,6356103.039d0,
     . 6356036.143d0,6355834.8467d0,6356773.3205d0,6356774.719d0,
     . 6356863.0188d0,6356794.343479d0,6356784.283666d0,
     . 6356768.337303d0/
c
      data ak0/0.9996d0/
c
      data radsec/206264.8062470964d0/
c
      slat = lat * secs
      slon = lon * secs
      slon = - slon
c
      a = axis(ispher)
      b = bxis(ispher)
      es= (a**2-b**2)/a**2
c
      izone = iizone
c
c     compute utm zone(izone) and central meridian in seconds for
c     geodetic to utm conversion where zone is not input.
c
      if(izone .eq. 0) then
         iiz=int(slon/2.16d4)
         izone = 30 - iiz
         if (slon.lt.0.0d0) izone=31- iiz
         if(iabs(izone).gt.30) then
            iutz = iabs(izone)-30
            cm=((iutz*6.0d0)-3.0d0)*(-3600.0d0)
         else
            iutz = 30-iabs(izone)
            cm=((iutz*6.0d0)+3.0d0)*3600.0d0
         endif
      else
         if(iabs(izone).gt.30) then
            iutz = iabs(izone)-30
            cm=((iutz*6.0d0)-3.0d0)*(-3600.0d0)
         else
            iutz = 30-iabs(izone)
            cm=((iutz*6.0d0)+3.0d0)*3600.0d0
         endif
      endif
      if (slat.lt.0.0d0) izone=-iabs(izone)
      if (izone.lt.0) slat = -dabs(slat)
      if (iizone.eq.0) iizone=izone
      if (slat.lt.0.0d0.and.iizone.gt.0) iizone=-iizone
c
c---------------------------------------------------------------
c     forward computation
c---------------------------------------------------------------
c
         phi = slat/radsec
         dlam = -(slon-cm)/radsec
c
         epri = es/(1.0d0-es)
         en = a/dsqrt(1.0d0-es*dsin(phi)**2)
         t = dtan(phi)**2
         c = epri*dcos(phi)**2
         aa = dlam*dcos(phi)
         s2 = dsin(2.0d0*phi)
         s4 = dsin(4.0d0*phi)
         s6 = dsin(6.0d0*phi)
         f1 = (1.d0-(es/4.d0)-(3.d0*es*es/6.4d1)
     .       -(5.d0*es*es*es/2.56d2))
         f2 = ((3.d0*es/8.d0)+(3.d0*es*es/3.2d1)
     .       +(4.5d1*es*es*es/1.024d3))
         f3 = ((1.5d1*es*es/2.56d2)+(4.5d1*es*es*es/1.024d3))
         f4 = (3.5d1*es*es*es/3.072d3)
         em = a*(f1*phi - f2*s2 + f3*s4 - f4*s6)
c
         xx = ak0*en*(aa + (1.d0-t+c)*aa**3/6.d0 + (5.d0-1.8d1*t+t*t+
     .        7.2d1*c-5.8d1*epri)*aa**5/1.2d2)
         xx = xx + 5.0d5 
c
         yy = ak0*(em + en*dtan(phi)*((aa*aa/2.d0) +
     .        (5.d0-t+9.d0*c+4.d0*c*c)*aa**4/2.4d1 +
     .        (6.1d1-5.8d1*t+t*t+6.d2*c-3.3d2*epri)*aa**6/7.2d2))
         if (izone.lt.0 .or. slat.lt.0.0d0) yy = yy + 1.0d7
c
c
c    test to see if within definition of utm projection
c
c  ---- latitude within 84 degrees
      if (dabs(slat).gt.302400) then
         write(99,*) ' *************   error exit from utm. '
         write(99,*) '   latitude value is poleward of 84 degrees'
         write(99,*) '   utm is not valid.'
         write(99,*) '   calculation has continued but values '
         write(99,*) '   may not be valid.'
         write(99,*) slat/3600.d0, slon/3600.d0
      endif
c  ---- delta longitude within 0.16 radians of central meridian
      if (dabs(dlam).gt.1.6d-1) then
         write(99,*) ' *************   error exit from utm. '
         write(99,*) '  d lon not within 0.16 radians of central merid.'
         write(99,*) '   calculation has continued but values '
         write(99,*) '   may not be valid.'
         write(99,*) slat/3600.d0, slon/3600.d0
      endif
c
      x = xx
      y = yy
c
c      write(*,*) '  ll2utm: x, y = ',x,y
      return
      end 
