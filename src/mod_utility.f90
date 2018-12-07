
!!!=== Copyright (c) 2018 Takashi NAKAMURA  =====

!!!**** UTILITY MODULE ****************

      MODULE mod_utility

      implicit none
      real(8), parameter :: PI = 3.14159265359d0
      real(8), parameter :: deg2rad = PI/180.0d0
          
      CONTAINS

! --- Compute coriolis parameter ---------------------
 
      real(8) function Coriolis(lat)
      real(8), intent(in) :: lat      ! Latitude (degree)
      real(8) :: omega
      
      omega = 2.0d0*PI*366.25d0/(24.0d0*3600.0d0*365.25d0)
      Coriolis = 2.0d0*omega * sin( deg2rad*lat )

      end function Coriolis
!
! --- Compute distance ---------------------
 
      real(8) function distance(lat1, lon1, lat2, lon2)
      real(8), intent(in) :: lat1      ! Latitude (degree) at point1
      real(8), intent(in) :: lon1      ! Longitude (degree) at point1
      real(8), intent(in) :: lat2      ! Latitude (degree) at point2
      real(8), intent(in) :: lon2      ! Longitude (degree) at point2

      real(8), parameter :: Re = 6378137.0d0  ! equatorial radius
      real(8), parameter :: e2 = 0.00669437999019758  ! eccentricity^2
        
      real(8) :: radlat1, radlon1, radlat2, radlon2
      real(8) :: d_radlat, d_radlon, m_radlat
      real(8) :: W, M, N, t1, t2
      
      radlat1 = deg2rad*lat1
      radlon1 = deg2rad*lon1
      radlat2 = deg2rad*lat2
      radlon2 = deg2rad*lon2
      
      d_radlat = radlat1 - radlat2
      d_radlon = radlon1 - radlon2
      m_radlat = 0.5d0*(radlat1 + radlat2)
      
      W = sqrt(1.0d0-e2*sin(m_radlat)*sin(m_radlat))
      M = Re*(1-e2)/(W**3.0d0)
      N = Re/W
      t1 = M * d_radlat
      t2 = N * cos(m_radlat)*d_radlon

      distance = sqrt(t1*t1+t2*t2)

      end function distance
!
! --- Seek min max values in 2D array ---------------------
 
      SUBROUTINE min_max_2D( Im, Jm, data, data_min, data_max)
      integer, intent(in ) :: Im
      integer, intent(in ) :: Jm
      real(8), intent(in ) :: data(Im,Jm)
      real(8), intent(out) :: data_min  ! minimum value
      real(8), intent(out) :: data_max  ! maximum value
      integer :: i,j,k

      data_max = data(1,1)
      data_min = data(1,1)
      do i=1, Im
        do j=1, Jm
          if(data_max < data(i,j)) data_max = data(i,j)
          if(data_min > data(i,j)) data_min = data(i,j)
        enddo
      enddo

      END SUBROUTINE min_max_2D
!
! --- Seek index ---------------------
 
      SUBROUTINE seek_id_range( Im, Xd, Xmin, Xmax, Irmin, Irmax)
      integer, intent(in ) :: Im
      real(8), intent(in ) :: Xd(Im)
      real(8), intent(in ) :: Xmin
      real(8), intent(in ) :: Xmax
      integer, intent(out) :: Irmin  ! Index of minimum
      integer, intent(out) :: Irmax  ! Index of minimum
      integer :: i,j,k
      
      if(Xd(1) < Xd(Im)) then
        Irmin=1
        do i=1,Im
          if(xd(i) > Xmin) exit
          Irmin=i
        end do
        do i=Irmin,Im
          Irmax=i
          if(Xd(i) > Xmax) exit
        end do
      else
        Irmin=Im
        do i=Im,1,-1
          if(xd(i) > Xmin) exit
          Irmin=i
        end do
        do i=Irmin,1,-1
          Irmax=i
          if(Xd(i) > Xmax) exit
        end do
      endif

      END SUBROUTINE seek_id_range

! --- Compute landmask mask_rho ---------------------
 
      SUBROUTINE land_masking(N_xi_rho, N_eta_rho, h, hmin, mask_rho)
      integer, intent(in ) :: N_xi_rho
      integer, intent(in ) :: N_eta_rho
      real(8), intent(in ) :: h(N_xi_rho,N_eta_rho) ! depth (m)
      real(8), intent(in ) :: hmin                  ! minimum depth (m)
      real(8), intent(out) :: mask_rho(N_xi_rho,N_eta_rho) ! Land/Sea mask on RHO-points.
      integer :: i,j,k
      integer :: ip,im,jp,jm
      integer :: count
      real(8) :: rmsk(N_xi_rho,N_eta_rho)
      
      mask_rho(:,:) = 1.0d0
      do i=1, N_xi_rho
        do j=1, N_eta_rho
          if(h(i,j) < hmin) then
            mask_rho(i,j) = 0.0d0
          endif
        enddo
      enddo
      ! Remove isolated water area
      rmsk(:,:) = 0.0d0
      rmsk(1:2,1:2) = 1.0d0  !!! seed
      do k=1,1000
        count = 0
        do i=1, N_xi_rho
          do j=1, N_eta_rho
            if(mask_rho(i,j)==1.0d0 .and. rmsk(i,j)==0.0d0)then 
              ip=i+1
              im=i-1
              jp=j+1
              jm=j-1
              if(ip==N_xi_rho+1)  ip=im
              if(im==0)           im=ip
              if(jp==N_eta_rho+1) jp=jm
              if(jm==0)           jm=jp
              if(rmsk(im,j)==1.0d0.or.rmsk(ip,j)==1.0d0.or.   &
                 rmsk(i,jm)==1.0d0.or.rmsk(i,jp)==1.0d0   )then 
                rmsk(i,j)=1.0d0
                count = count+1
              endif
            endif
          enddo
        enddo
        if(count==0) exit
      enddo
      
      ! Remove one grid bay
      do k=1,1000
        count = 0
        do i=1, N_xi_rho
          do j=1, N_eta_rho
            if(rmsk(i,j)==1.0d0)then 
              ip=i+1
              im=i-1
              jp=j+1
              jm=j-1
              if(ip==N_xi_rho+1)  ip=im
              if(im==0)           im=ip
              if(jp==N_eta_rho+1) jp=jm
              if(jm==0)           jm=jp
              if(rmsk(im,j)+rmsk(ip,j)+rmsk(i,jm)+rmsk(i,jp)==1.0d0)then 
                rmsk(i,j)=0.0d0
                count = count+1
              endif
            endif
          enddo
        enddo
        if(count==0) exit
      enddo
      
      mask_rho(:,:) = rmsk(:,:)

      END SUBROUTINE land_masking

! --- Compute mask_u, mask_v, mask_psi ---------------------
 
      SUBROUTINE uvp_masks(N_xi_rho, N_eta_rho, mask_rho, mask_u, mask_v, mask_psi)
      integer, intent(in ) :: N_xi_rho
      integer, intent(in ) :: N_eta_rho
      real(8), intent(in ) :: mask_rho(N_xi_rho,N_eta_rho) ! Land/Sea mask on PHO-points.
      real(8), intent(out) :: mask_u(N_xi_rho-1,N_eta_rho) ! Land/Sea mask on U-points.
      real(8), intent(out) :: mask_v(N_xi_rho,N_eta_rho-1) ! Land/Sea mask on V-points.
      real(8), intent(out) :: mask_psi(N_xi_rho-1,N_eta_rho-1) ! Land/Sea mask on PSI-points.
      integer :: i,j
      real(8) :: tm
      
      do i=1, N_xi_rho-1
        do j=1, N_eta_rho
          mask_u(i,j) = mask_rho(i+1,j)*mask_rho(i,j)
        enddo
      enddo
      
      do i=1, N_xi_rho
        do j=1, N_eta_rho-1
          mask_v(i,j) = mask_rho(i,j+1)*mask_rho(i,j)
        enddo
      enddo
      
      do i=1, N_xi_rho-1
        do j=1, N_eta_rho
          tm = (mask_rho(i+1,j+1)+mask_rho(i,j))   &
     &        *(mask_rho(i+1,j)+mask_rho(i,j+1))
          if(tm>1.5d0) then
            mask_psi(i,j) = 1.0d0
          elseif(tm>0.5d0 ) then 
            mask_psi(i,j) = 2.0d0
          else 
            mask_psi(i,j) = 0.0d0
          endif
        enddo
      enddo

      END SUBROUTINE uvp_masks
      
! --- Compute bathymetry smoothing ---------------------
      SUBROUTINE bathy_smooth(N_xi_rho, N_eta_rho, rx0max, mask_rho, h)
      integer, intent(in ) :: N_xi_rho
      integer, intent(in ) :: N_eta_rho
      real(8), intent(in ) :: rx0max
      real(8), intent(in ) :: mask_rho(N_xi_rho,N_eta_rho) ! Land/Sea mask on PHO-points.
      real(8), intent(inout) :: h(N_xi_rho,N_eta_rho) ! depth (m).
      real(8), parameter :: c1 = 1.0d-2
      integer :: i,j,k
      integer :: count
      integer :: ip,im,jp,jm
      real(8) :: rx0, MAX_rx0

      do k=1,1000
        MAX_rx0 = 0.0d0
        do i=1, N_xi_rho-1
          do j=1, N_eta_rho-1
            if(mask_rho(i,j)*mask_rho(i,j+1)==1.0d0) then
              rx0 = abs(h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1))
              if(rx0>rx0max) then
                h(i,j)=h(i,j)-c1*(h(i,j)-h(i,j+1))
                h(i,j+1)=h(i,j+1)+c1*(h(i,j)-h(i,j+1))
                if(rx0>MAX_rx0) MAX_rx0=rx0
              endif
            endif
            if(mask_rho(i,j)*mask_rho(i+1,j)==1.0d0) then
              rx0 = abs(h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j))
              if(rx0>rx0max) then
                h(i,j)=h(i,j)-c1*(h(i,j)-h(i+1,j))
                h(i+1,j)=h(i+1,j)+c1*(h(i,j)-h(i+1,j))
                if(rx0>MAX_rx0) MAX_rx0=rx0
              endif
            endif
          enddo
        enddo
        if(MAX_rx0 == 0.0d0) exit
        write(*,*) k,'maximun rx0:', MAX_rx0
      enddo
    
      END SUBROUTINE bathy_smooth

      END MODULE mod_utility

