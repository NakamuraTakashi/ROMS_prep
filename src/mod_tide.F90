
!!!=== Copyright (c) 2022 Takashi NAKAMURA  =====

!!!**** TIDE MODULE ************************************

MODULE mod_tide

  implicit none  

  CONTAINS

!============== SUBROUTINE for NAO99 ====================================

  SUBROUTINE readNAO99_head( iu, xmin, xmax, ymin, ymax, dx, dy         &
                     , mmax, nmax, ideff, fmt, aunit, punit )
  
    integer, intent(in ) :: iu            ! file unit number
    real(8), intent(out) :: xmin, xmax    ! 
    real(8), intent(out) :: ymin, ymax    ! 
    real(8), intent(out) :: dx, dy        ! 
    integer, intent(out) :: mmax, nmax    ! 
    integer, intent(out) :: ideff         ! 
    character(6), intent(out) :: fmt      !
    real(8), intent(out) :: aunit         ! 
    real(8), intent(out) :: punit         ! 

    character(50) :: ctmp
    integer :: idx,idy

    read(iu,'(13x,a50)') ctmp
    read(iu,'(13x,a3 )') ctmp
    read(iu,'(13x,f5.3,19x,f4.2)') aunit,punit
    read(iu,'(13x,a50)') ctmp
    read(iu,'(7x,f7.2,3(9x,f7.2))') xmin,xmax,ymin,ymax
    read(iu,'(12x,i2,14x,i2,2(9x,i7))') idx,idy,mmax,nmax
    read(iu,'(16x,i6,11x,a6)') ideff, fmt 
  
    if (idx.eq.50) then
       dx = 0.5d0
       dy = 0.5d0
    else
       dx = 1.d0/dfloat(idx)
       dy = 1.d0/dfloat(idy)
    endif

  END SUBROUTINE

!--------------------------------------------------------------------
  SUBROUTINE readNAO99_data( iu, mmax, nmax, fmt, aunit, punit  &
                           , amp, phs  )

    integer, intent(in ) :: iu            ! file unit number
    integer, intent(in ) :: mmax, nmax    ! 
    character(6), intent(in ) :: fmt      ! Model name
    real(8), intent(in ) :: aunit         ! 
    real(8), intent(in ) :: punit         ! 
    real(8), intent(out) :: amp(mmax,nmax)
    real(8), intent(out) :: phs(mmax,nmax)
  
    integer, parameter :: kc = 10
    integer :: iamp(mmax), iphs(mmax)
    integer :: kend, krem
    integer :: k,m,n
    integer :: m1,m2
   
  ! -----< Reading loop >-----
  
    kend = mmax/kc
    krem = mod(mmax,kc)
    do n = nmax,1,-1 !!!!!!!!!!!!!  Check 
       do k = 1,kend
          m1 = (k-1)*kc + 1
          m2 = m1 + kc - 1
          read(iu,fmt) (iamp(m),m=m1,m2)
       enddo
       if (krem.ne.0) then
          m1 = kend*kc + 1
          m2 = kend*kc + krem
          read(iu,fmt) (iamp(m),m=m1,m2)
       endif
       do k = 1,kend
          m1 = (k-1)*kc + 1
          m2 = m1 + kc - 1
          read(iu,fmt) (iphs(m),m=m1,m2)
       enddo
       if (krem.ne.0) then
          m1 = kend*kc + 1
          m2 = kend*kc + krem
          read(iu,fmt) (iphs(m),m=m1,m2)
       endif
       do m = 1,mmax
          amp(m,n) = dfloat(iamp(m))*aunit*0.01d0  ! in meters
          phs(m,n) = dfloat(iphs(m))*punit  ! in degrees
       enddo
    enddo
         
  END SUBROUTINE

!--------------------------------------------------------------------
  SUBROUTINE readNAO99_vel( iu, mmax, nmax, lat, lon  &
                           , Uamp, Uphs, Vamp, Vphs  )

    integer, intent(in ) :: iu            ! file unit number
    integer, intent(in ) :: mmax, nmax    ! 
    real(8), intent(in ) :: lat(mmax,nmax)
    real(8), intent(in ) :: lon(mmax,nmax)
    real(8), intent(out) :: Uamp(mmax,nmax)
    real(8), intent(out) :: Uphs(mmax,nmax)
    real(8), intent(out) :: Vamp(mmax,nmax)
    real(8), intent(out) :: Vphs(mmax,nmax)

    integer :: ios
    real(8) :: x,y,Ua,Up,Va,Vp
    integer :: m,n
    real(8) :: dx,dy
    
    Uamp(:,:) = 9999.0d0
    Uphs(:,:) = 9999.0d0
    Vamp(:,:) = 9999.0d0
    Vphs(:,:) = 9999.0d0

    LOOP: DO
      read(iu,'(2f8.3,4d11.4)',iostat=ios) x,y,Ua,Up,Va,Vp
      if(ios==-1) exit

      if(Ua==0.0d0 .or. Va==0.0d0) cycle LOOP

      do n = nmax,1,-1
        do m = 1,mmax
          dx = abs(lon(m,n)-x)
          dy = abs(lat(m,n)-y)
          if(dx<=0.01d0 .and. dy<=0.01d0) then
            Uamp(m,n) = Ua
            Uphs(m,n) = Up
            Vamp(m,n) = Va
            Vphs(m,n) = Vp
            cycle LOOP
          endif
        end do
      end do 
    END DO LOOP
         
  END SUBROUTINE
!--------------------------------------------------------------------
  SUBROUTINE Tvel_nao2roms( Au, Pu, Av, Pv, Pini               &
                          , Cphase, Cangle, Cmax, Cmin )
  
    real(8), intent(in ) :: Au, Av ! cm/s
    real(8), intent(in ) :: Pu, Pv ! degree
    real(8), intent(in ) :: Pini     ! degree
    real(8), intent(out) :: Cphase ! degree
    real(8), intent(out) :: Cangle ! degree
    real(8), intent(out) :: Cmax   ! m/s
    real(8), intent(out) :: Cmin   ! m/s

    real(8), parameter :: PI = 3.14159265359d0
    real(8), parameter :: deg2rad = PI/180.0d0

    real(8) :: cff
    real(8) :: tan2A, dP
    real(8) :: a,b,c,cmax2,cmin2
    real(8) :: Px,Py

    Px = Pini-Pu
    Py = Pini-Pv

    dP = (Py-Px)*deg2rad
    dP = mod(dP+2.0d0*PI,2.0d0*PI) ! 0 < dP < 2*PI

    a = 1.0d0
    b = -(Au*Au+Av*Av)
    c = Au*Au*Av*Av*sin(dP)*sin(dP)

    cmax2 = (-b+sqrt(b*b-4.0d0*a*c))/(2.0d0*a)
    cmin2 = Au*Au+Av*Av - cmax2

    Cmax = sqrt(cmax2)*0.01d0
    Cmin = -sqrt(cmin2)*0.01d0
    if(dP>PI) then
      Cmin = -Cmin
    endif

    if(Au==Av) then
      Cangle = sign(45.0d0,cos(dP))
    else
      tan2A = 2*Au*Av/(Au*Au-Av*Av)*cos(dP)
      cff = atan(tan2A)
      Cangle = 0.5d0*cff/deg2rad
      if(Au<Av) then
        Cangle = Cangle + 90.0d0
      endif
      if(Cangle<0.0d0) then
        Cangle = 180.0d0 + Cangle
      endif
    endif

    Cphase = mod(-Py+360.0d0, 360.0d0)


  END SUBROUTINE

!--------------------------------------------------------------------
    
END MODULE mod_tide
      
! -------------------------------------------------------------------------
