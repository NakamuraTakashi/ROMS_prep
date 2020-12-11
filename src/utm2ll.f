      subroutine utm2ll(x,y,iizone,ispher,lat,lon,conv)
c
c     universal transverse mercator conversion
c
c     rewritten 6-16-83 by j.f. waananen using formulation
c     by j.p. snyder in bulletin 1532, pages 68-69
c
c     convert utm.f to convert from utm to lat lon (idir > 0)
c      2/5/2001  j. klinck
c 
      implicit none
c
      real*8 axis(19),bxis(19)
      real*8 lat,lon
      real*8 radsec
      real*8 ak0,a,b,es
      real*8 x,y,slat,slon,conv,cm,phi,dlam,epri,en,t
      real*8 c,em,xx,yy,um,e1
      real*8 phi1,r,d,alam,secs
c
      integer iizone,ispher,izone
      integer iutz

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
      a = axis(ispher)
      b = bxis(ispher)
      es= (a**2-b**2)/a**2
c
      izone = iizone
c
c     compute utm zone(izone) and central meridian in seconds for
c     geodetic to utm conversion where zone is not input.
c
      if (izone.eq.0 ) then
         write(*,*) ' *************   error exit from utm. '
         write(*,*) '  zone is not given for conversion from '
         write(*,*) '  x,y to lat,lon. '
         return
      endif
c
      if(iabs(izone).gt.30) then
         iutz = iabs(izone)-30
         cm=((iutz*6.0d0)-3.0d0)*(-3600.0d0)
      else
         iutz = 30-iabs(izone)
         cm=((iutz*6.0d0)+3.0d0)*3600.0d0
      endif
c
c---------------------------------------------------------------
c     inverse computation
c---------------------------------------------------------------
         yy = y
         if (izone.lt.0) yy = yy - 1.0d7
         xx = x - 5.0d5 
         em = yy/ak0
         um = em/(a*(1.d0-(es/4.d0)-(3.d0*es*es/6.4d1)-
     .        (5.d0*es*es*es/2.56d2)))
         e1 = (1.d0-dsqrt(1.d0-es))/(1.d0+dsqrt(1.d0-es))
c
         phi1 = um+((3.d0*e1/2.d0)-(2.7d1*e1**3/3.2d1))*dsin(2.d0*um) +
     .        ((2.1d1*e1*e1/1.6d1)-(5.5d1*e1**4/3.2d1))*dsin(4.d0*um) +
     .        (1.51d2*e1**3/9.6d1)*dsin(6.d0*um)
c
         en = a/dsqrt(1.0d0-es*dsin(phi1)**2)
         t  = dtan(phi1)**2
         epri = es/(1.d0-es)
         c  = epri*dcos(phi1)**2
         r  = (a*(1.d0-es))/((1.d0-es*dsin(phi1)**2)**1.5d0)
         d  = xx/(en*ak0)   
c
         phi = phi1 - (en*dtan(phi1)/r) * ((d*d/2.d0) -
     .        (5.d0+3.d0*t+10.d0*c-4.d0*c*c-9.d0*epri)*d**4/2.4d1 
     .      + (6.1d1+9.d1*t+2.98d2*c+4.5d1*t*t
     .        -2.52d2*epri-3.d0*c*c)*d**6/7.2d2)
c
         alam = (cm/radsec)-(d-(1.d0+2.d0*t+c)*d**3/6.d0 + 
     .        (5.d0-2.d0*c+2.8d1*t -3.d0*c*c+8.d0*epri+2.4d1*t*t)
     .        *d**5/1.2d2)/dcos(phi1)
c
         slat = phi*radsec
         slon = alam*radsec
         dlam = -(slon-cm)/radsec

         lat = slat / secs
         lon = -slon / secs
         if (abs(lon).gt.180.0) lon=360.0+lon

c---------------------------------------------------------------
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
c
c     compute convergence
c
          conv = dlam*(dsin(phi)+1.9587d-12*(dlam**2)
     .         *dsin(phi)*dcos(phi*phi))
          conv = conv*radsec
c     
      return
      end 
