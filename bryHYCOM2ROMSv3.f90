
!!!=== ver 2016/03/22   Copyright (c) 2014-2016 Takashi NAKAMURA  =====

    PROGRAM bryHYCOM2ROMS
      use netcdf
      use mod_calendar
      use mod_interpolation
     
      implicit none
      
! -------------------------------------------------------------------------
!      integer, parameter :: N_s_rho = 30
      integer, parameter :: N_s_rho = 15
      
      integer, parameter :: spherical  = 1   
      
! Set vertical, terrain-following coordinates transformation equation and
! stretching function (see below for details), [1:Ngrids].

      integer, parameter :: Vtransform  = 2   !!! transformation equation
      integer, parameter :: Vstretching = 4   !!! stretching function

! Vertical S-coordinates parameters (see below for details), [1:Ngrids].

      real(8), parameter :: THETA_S = 7.0d0   !!! 1.0d0                      ! surface stretching parameter
      real(8), parameter :: THETA_B = 0.1d0   !!! 1.0d0                      ! bottom  stretching parameter
      real(8), parameter :: TCLINE  = 200.0d0  !!! 0.0d0                     ! critical depth (m)
      real(8), parameter :: DCRIT   = 0.10d0  !!! 0.0d0                      ! critical depth (m)
! -------------------------------------------------------------------------
      
      integer, parameter :: Ryear  = 2000   ! Reference year
      integer, parameter :: Rmonth = 1      ! Reference month
      integer, parameter :: Rday   = 1      ! Reference day
      
      integer, parameter :: mode = 1        ! mode=1, linear intrtporation
      logical, parameter :: write_ini = .FALSE.   ! .TRUE. : write INI_FILE; .FALSE.: not write
      
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "D:/ROMS/Yaeyama/Data/Yaeyama1_grd_v9.nc"
      character(len=*), parameter :: OCEAN_FILE  =                      &
     &     "D:/ROMS/Yaeyama/Data/HYCOM_extracted_data/Yaeyama1_HYCOM_extracted_1509.nc"
      character(len=*), parameter :: BRY_FILE  =                        &
!     &     "D:/ROMS/Yaeyama/Data/Yaeyama1_bry_Nz30_HYCOM_0709.nc"
     &     "D:/ROMS/Yaeyama/Data/Yaeyama1_bry_Nz15_HYCOM_1509.nc"
      character(len=*), parameter :: INI_FILE  =                        &
     &     "D:/ROMS/Yaeyama/Data/Yaeyama1_ini_Nz30_HYCOM_1504.nc"
      
! -------------------------------------------------------------------------

      character(33) :: TIME_ATT  = "seconds since 1992-01-01 00:00:00"
      
      real(8), parameter :: grav = 9.80665d0
      real(8), parameter :: dynz = 0.0d0
      real(8), parameter :: ref_press = 0.0d0
      real(8), parameter :: PI = 3.14159265359d0
      real(8) :: press, potemp

      real(8), allocatable :: time_all(:), lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:), depth(:),time(:),tau(:)
      real(8), allocatable :: time2(:)
      real(8), allocatable :: water_u(:, :, :, :), water_v(:, :, :, :)
      real(8), allocatable :: salinity(:, :, :, :),water_temp(:, :, :, :)
      real(8), allocatable :: surf_el(:, :, :)
      integer :: start1D(1), count1D(1)
      integer :: start2D(2), count2D(2)
      integer :: start3D(3), count3D(3)
      integer :: start4D(4), count4D(4)
      
      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      real(8), allocatable :: lat_u(:, :)
      real(8), allocatable :: lon_u(:, :)
      real(8), allocatable :: lat_v(:, :)
      real(8), allocatable :: lon_v(:, :)
      real(8), allocatable :: h(:,:)        ! depth (meter)
      real(8), allocatable :: h_u(:,:)        ! depth (meter)
      real(8), allocatable :: h_v(:,:)        ! depth (meter)
      real(8), allocatable :: z_w(:,:,:)  ! vertical depths z_w(x,y,0:N(ng)) at  W-points (top/bottom cell) (meter)
      real(8), allocatable :: z_r(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
      real(8), allocatable :: z_r_u(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
      real(8), allocatable :: z_r_v(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
     
      real(8) :: ocean_time(1)            ! Ocean time
      real(8), allocatable :: zeta(:,:)   ! free-surface (meter)
      real(8), allocatable :: ubar(:,:)   ! vertically integrated u-momentum component (meter second-1)
      real(8), allocatable :: vbar(:,:)   ! vertically integrated v-momentum component (milibar=hPa)
      real(8), allocatable :: u(:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: v(:,:,:)    ! v-momentum component (meter second-1)
      real(8), allocatable :: temp(:,:,:) ! potential temperature (Celsius)
      real(8), allocatable :: salt(:,:,:) ! salinity (psu)    

      real(8), allocatable :: cosAu(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: sinAu(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: cosAv(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: sinAv(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: ull_u(:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: vll_u(:,:,:)    ! v-momentum component (meter second-1)
      real(8), allocatable :: ull_v(:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: vll_v(:,:,:)    ! v-momentum component (meter second-1)
      
      real(8) :: hc       
      real(8) :: sc_w(0:N_s_rho)       
      real(8) :: sc_r(1:N_s_rho)  
      real(8) :: Cs_w(0:N_s_rho)       
      real(8) :: Cs_r(1:N_s_rho)  
      real(8) :: d_lat, d_lon
      
      integer :: i,j,k
      integer :: itime

      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_xi_u, N_eta_u
      integer :: N_xi_v, N_eta_v
      
      integer :: ncid,var_id
      integer :: ncid2,var_id2
      
      real(8) :: sf, off
      integer :: Im, Jm, Nz, Nt
      integer :: IL, IR, JB, JT
      integer :: ItS, ItE

      
!      integer :: time_varid, Uwind_varid
      
      
      write (YYYY, "(I4.4)") Ryear
      write (MM, "(I2.2)") Rmonth
      write (DD, "(I2.2)") Rday
      
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(15:18)=YYYY
      TIME_ATT(20:21)=MM
      TIME_ATT(23:24)=DD
      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      call get_dimension(ncid, 'xi_u',  N_xi_u)
      call get_dimension(ncid, 'eta_u', N_eta_u)
      call get_dimension(ncid, 'xi_v',  N_xi_v)
      call get_dimension(ncid, 'eta_v', N_eta_v)
      
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(lat_u(N_xi_u, N_eta_u))
      allocate(lon_u(N_xi_u, N_eta_u))
      allocate(lat_v(N_xi_v, N_eta_v))
      allocate(lon_v(N_xi_v, N_eta_v))
      allocate(h(N_xi_rho, N_eta_rho))
      allocate(h_u(N_xi_u, N_eta_u))
      allocate(h_v(N_xi_v, N_eta_v))
      allocate(z_w(N_xi_rho, N_eta_rho, 0:N_s_rho))
      allocate(z_r(N_xi_rho, N_eta_rho, 1:N_s_rho))
      allocate(z_r_u(N_xi_u, N_eta_u, 1:N_s_rho))
      allocate(z_r_v(N_xi_v, N_eta_v, 1:N_s_rho))
      
      allocate(zeta(N_xi_rho, N_eta_rho))
      allocate(ubar(N_xi_u, N_eta_u))
      allocate(vbar(N_xi_v, N_eta_v))
      allocate(u(N_xi_u, N_eta_u, N_s_rho))
      allocate(v(N_xi_v, N_eta_v, N_s_rho))
      allocate(temp(N_xi_rho, N_eta_rho, N_s_rho))
      allocate(salt(N_xi_rho, N_eta_rho, N_s_rho))
      
      allocate(cosAu(N_xi_u, N_eta_u))
      allocate(sinAu(N_xi_u, N_eta_u))
      allocate(cosAv(N_xi_v, N_eta_v))
      allocate(sinAv(N_xi_v, N_eta_v))
      allocate(ull_u(N_xi_u, N_eta_u, N_s_rho))
      allocate(vll_u(N_xi_u, N_eta_u, N_s_rho))
      allocate(vll_v(N_xi_v, N_eta_v, N_s_rho))
      allocate(ull_v(N_xi_v, N_eta_v, N_s_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      call check( nf90_inq_varid(ncid, 'lat_u', var_id) ) ! latitude at U-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_u) )
      call check( nf90_inq_varid(ncid, 'lon_u', var_id) ) ! longitude at U-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_u) )
      call check( nf90_inq_varid(ncid, 'lat_v', var_id) ) ! latitude at V-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_v) )
      call check( nf90_inq_varid(ncid, 'lon_v', var_id) ) ! longitude at V-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_v) )
      call check( nf90_inq_varid(ncid, 'h', var_id) ) 
      call check( nf90_get_var(ncid, var_id, h) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", GRID_FILE
      
      
      
      call set_scoord (             &
!        input parameters
     &              N_s_rho         &
     &            , Vtransform      &
     &            , Vstretching     &
     &            , THETA_S         &
     &            , THETA_B         &
     &            , TCLINE          &
     &            , DCRIT           &
!        output parameters
     &            , hc              &
     &            , sc_w            &
     &            , sc_r            &
     &            , Cs_w            &
     &            , Cs_r            &
     &           )

      
!      write(*,*) hc
!      open(50, file='check1.dat')
!      do i=0, N_s_rho
!        write(50,*) sc_w(i), Cs_w(i)
!      enddo
!      close(50)
!      open(50, file='check2.dat')
!      do i=1, N_s_rho
!        write(50,*) sc_r(i), Cs_r(i)
!      enddo
!      close(50)
      
      do i=1,N_xi_u
        do j=1,N_eta_u
          h_u(i,j)= (h(i,j)+h(i+1,j))*0.5d0
        enddo
      enddo
      do i=1,N_xi_v
        do j=1,N_eta_v
          h_v(i,j)= (h(i,j)+h(i,j+1))*0.5d0
        enddo
      enddo
      
      do i=1,N_xi_rho
        do j=1,N_eta_rho
          do k=1,N_s_rho
            z_r(i,j,k)= -h(i,j)*Cs_r(k)
          enddo
        enddo
      enddo
      do i=1,N_xi_u
        do j=1,N_eta_u
          do k=1,N_s_rho
            z_r_u(i,j,k)= -h_u(i,j)*Cs_r(k)
          enddo
        enddo
      enddo
      do i=1,N_xi_v
        do j=1,N_eta_v
          do k=1,N_s_rho
            z_r_v(i,j,k)= -h_v(i,j)*Cs_r(k)
          enddo
        enddo
      enddo
      do i=1,N_xi_rho
        do j=1,N_eta_rho
          do k=0,N_s_rho
            z_w(i,j,k)= -h(i,j)*Cs_w(k)
          enddo
        enddo
      enddo
      
      do i=1,N_xi_u
        do j=1,N_eta_u
          d_lat=lat_rho(i+1,j)-lat_rho(i,j)
          d_lon=lon_rho(i+1,j)-lon_rho(i,j)
          d_lon=d_lon*cos(lat_rho(i,j)/180.0d0*PI)
          cosAu(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
          sinAu(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
        enddo
      enddo
      do i=1,N_xi_v
        do j=1,N_eta_v
          d_lat=lat_rho(i,j+1)-lat_rho(i,j)
          d_lon=lon_rho(i,j)-lon_rho(i,j+1)
          d_lon=d_lon*cos(lat_rho(i,j)/180.0d0*PI)
          cosAv(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
          sinAv(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
        enddo         
      enddo

      
      if ( write_ini ) then
      
!---- Create the ROMS initial conditions netCDF file --------------------------------
          
        call createNetCDFini(           &
!            input parameters
     &          INI_FILE               &
     &        , TIME_ATT               &  
     &        , N_xi_rho, N_eta_rho, N_s_rho, 1          &   
     &        )
        
        call check( nf90_open(INI_FILE, NF90_WRITE, ncid) )
        call check( nf90_inq_varid(ncid, 'spherical', var_id) )
        call check( nf90_put_var(ncid, var_id, spherical) )
        call check( nf90_inq_varid(ncid, 'Vtransform', var_id) )
        call check( nf90_put_var(ncid, var_id, Vtransform ) )
        call check( nf90_inq_varid(ncid, 'Vstretching', var_id) )
        call check( nf90_put_var(ncid, var_id, Vstretching ) )
        call check( nf90_inq_varid(ncid, 'theta_s', var_id) )
        call check( nf90_put_var(ncid, var_id, THETA_S ) )
        call check( nf90_inq_varid(ncid, 'theta_b', var_id) )
        call check( nf90_put_var(ncid, var_id, THETA_B) )
        call check( nf90_inq_varid(ncid, 'Tcline', var_id) )
        call check( nf90_put_var(ncid, var_id, TCLINE) )
        call check( nf90_inq_varid(ncid, 'hc', var_id) )
        call check( nf90_put_var(ncid, var_id, hc ) )
        call check( nf90_close(ncid) )
      
!---- Write ROMS initial conditions netCDF file --------------------------------
      
        start1D = (/ 1 /)
        count1D = (/ N_s_rho /)
        call writeNetCDF_1d(           &
!            input parameters
     &          's_rho'                &
     &        , INI_FILE               &
     &        , N_s_rho                &
     &        , sc_r                   &
     &        , start1D, count1D       &
     &        )
        call writeNetCDF_1d(           &
!            input parameters
     &          'Cs_r'                 &
     &        , INI_FILE               &
     &        , N_s_rho                &
     &        , Cs_r                   &
     &        , start1D, count1D       &
     &        )
     
        start1D = (/ 1 /)
        count1D = (/ N_s_rho+1 /)
        call writeNetCDF_1d(           &
!            input parameters
     &          's_w'                  &
     &        , INI_FILE               &
     &        , N_s_rho+1              &
     &        , sc_w                   &
     &        , start1D, count1D       &
     &        )
        call writeNetCDF_1d(           &
!            input parameters
     &          'Cs_w'                 &
     &        , INI_FILE               &
     &        , N_s_rho+1              &
     &        , Cs_w                   &
     &        , start1D, count1D       &
     &        )
     
      end if
      
!---- Create ROMS boundary conditions netCDF file --------------------------------

      call createNetCDFbry(           &
!          input parameters
     &        BRY_FILE               &
     &      , TIME_ATT               &  
     &      , N_xi_rho, N_eta_rho, N_s_rho, 1          &   
     &      )
     
      call check( nf90_open(BRY_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, 'spherical', var_id) )
      call check( nf90_put_var(ncid, var_id, spherical) )
      call check( nf90_inq_varid(ncid, 'Vtransform', var_id) )
      call check( nf90_put_var(ncid, var_id, Vtransform ) )
      call check( nf90_inq_varid(ncid, 'Vstretching', var_id) )
      call check( nf90_put_var(ncid, var_id, Vstretching ) )
      call check( nf90_inq_varid(ncid, 'theta_s', var_id) )
      call check( nf90_put_var(ncid, var_id, THETA_S ) )
      call check( nf90_inq_varid(ncid, 'theta_b', var_id) )
      call check( nf90_put_var(ncid, var_id, THETA_B) )
      call check( nf90_inq_varid(ncid, 'Tcline', var_id) )
      call check( nf90_put_var(ncid, var_id, TCLINE) )
      call check( nf90_inq_varid(ncid, 'hc', var_id) )
      call check( nf90_put_var(ncid, var_id, hc ) )
      call check( nf90_close(ncid) )
      
!---- Write ROMS boundary conditions netCDF file --------------------------------
      
      start1D = (/ 1 /)
      count1D = (/ N_s_rho /)
      call writeNetCDF_1d(           &
!          input parameters
     &        's_rho'                &
     &      , BRY_FILE               &
     &      , N_s_rho                &
     &      , sc_r                   &
     &      , start1D, count1D       &
     &      )
      call writeNetCDF_1d(           &
!          input parameters
     &        'Cs_r'                 &
     &      , BRY_FILE               &
     &      , N_s_rho                &
     &      , Cs_r                   &
     &      , start1D, count1D       &
     &      )
     
      start1D = (/ 1 /)
      count1D = (/ N_s_rho+1 /)
      call writeNetCDF_1d(           &
!          input parameters
     &        's_w'                  &
     &      , BRY_FILE               &
     &      , N_s_rho+1              &
     &      , sc_w                   &
     &      , start1D, count1D       &
     &      )
      call writeNetCDF_1d(           &
!          input parameters
     &        'Cs_w'                 &
     &      , BRY_FILE               &
     &      , N_s_rho+1              &
     &      , Cs_w                   &
     &      , start1D, count1D       &
     &      )
     
      
!---- Read extracted HYCOM netCDF file --------------------------------

      ! Open NetCDF file
      write(*,*) "OPEN: ", OCEAN_FILE
      call check( nf90_open(OCEAN_FILE, nf90_nowrite, ncid) )
      call get_dimension(ncid, 'time', Nt)
      allocate(time_all(Nt))
      call check( nf90_inq_varid(ncid, 'time', var_id) )
      call check( nf90_get_var(ncid, var_id, time_all) )
      ! Get dimension data
      call get_dimension(ncid, 'lat', Jm)
      call get_dimension(ncid, 'lon', Im)
      call get_dimension(ncid, 'depth', Nz)
      
!      ! Allocate variable
      allocate(lat(Jm))
      allocate(lon(Im))
      allocate(depth(Nz))
      allocate(surf_el(Im,Jm,1))
      allocate(water_temp(Im,Jm,Nz,1))
      allocate(salinity(Im,Jm,Nz,1))
      allocate(water_u(Im,Jm,Nz,1))
      allocate(water_v(Im,Jm,Nz,1))

      call check( nf90_inq_varid(ncid, 'lat', var_id) )
      call check( nf90_get_var(ncid, var_id, lat) )
      call check( nf90_inq_varid(ncid, 'lon', var_id) )
      call check( nf90_get_var(ncid, var_id, lon) )
      call check( nf90_inq_varid(ncid, 'depth', var_id) )
      call check( nf90_get_var(ncid, var_id, depth) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      write(*,*) "CLOSE: ", OCEAN_FILE
      write(*,*) "*********************************************"
      
!==== LOOP START ======================================================

      do itime=1, Nt
      
!---- Read extracted HYCOM netCDF file --------------------------------

        write(*,*) 'time (days) = ', time_all(itime)

        write(*,*) "OPEN: ", OCEAN_FILE
        call check( nf90_open(OCEAN_FILE, nf90_nowrite, ncid) )
        
        start3D = (/ 1,  1,  itime /)
        count3D = (/ Im, Jm, 1     /)
        call check( nf90_inq_varid(ncid, 'surf_el', var_id) )
        call check( nf90_get_var(ncid, var_id, surf_el, start=start3D, count=count3D)  )

        start4D = (/ 1,  1,  1,  itime /)
        count4D = (/ Im, Jm, Nz, 1     /)
        call check( nf90_inq_varid(ncid, 'water_temp', var_id) )
        call check( nf90_get_var(ncid, var_id, water_temp, start=start4D, count=count4D)  )
        call check( nf90_inq_varid(ncid, 'salinity', var_id) )
        call check( nf90_get_var(ncid, var_id, salinity, start=start4D, count=count4D)  )
        call check( nf90_inq_varid(ncid, 'water_u', var_id) )
        call check( nf90_get_var(ncid, var_id, water_u, start=start4D, count=count4D)  )
        call check( nf90_inq_varid(ncid, 'water_v', var_id) )
        call check( nf90_get_var(ncid, var_id, water_v, start=start4D, count=count4D)  )
        
        ! Close NetCDF file
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", OCEAN_FILE
        
! --------------------------------------------------------

        ocean_time(1) = time_all(itime)*24.0d0*60.0d0*60.0d0  !! day -> second

! --- Linear interporation --------------------------------
      
        if (mode==1) then
        
          write(*,*) 'Linear Interporation'
        
          call LinearInterpolation2D_grid2(Im, Jm, lon, lat               &
     &                    , surf_el(:,:,1), -5.0d0, 5.0d0                 &
     &                    , N_xi_rho, N_eta_rho, lon_rho, lat_rho, zeta )
        
          call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                    , water_temp(:,:,:,1),  0.0d0, 40.0d0           &
     &                    , N_xi_rho, N_eta_rho, N_s_rho                  &
     &                    , lon_rho, lat_rho, z_r                         &
     &                    , temp )
     
          call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                    , salinity(:,:,:,1),  20.0d0, 36.0d0            &
     &                    , N_xi_rho, N_eta_rho, N_s_rho                  &
     &                    , lon_rho, lat_rho, z_r                         &
     &                    , salt )

        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , water_u(:,:,:,1),  -5.0d0, 5.0d0              &
     &                  , N_xi_u, N_eta_u, N_s_rho                      &
     &                  , lon_u, lat_u, z_r_u                           &
     &                  , ull_u(:,:,:) )
        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , water_v(:,:,:,1),  -5.0d0, 5.0d0              &
     &                  , N_xi_u, N_eta_u, N_s_rho                      &
     &                  , lon_u, lat_u, z_r_u                           &
     &                  , vll_u(:,:,:) )
        do i=1,N_xi_u
          do j=1,N_eta_u
            do k=1,N_s_rho
              u(i,j,k) = ull_u(i,j,k)*cosAu(i,j)                        &
                        +vll_u(i,j,k)*sinAu(i,j)
            enddo
          enddo
        enddo
        

        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , water_v(:,:,:,1),  -5.0d0, 5.0d0              &
     &                  , N_xi_v, N_eta_v, N_s_rho                      &
     &                  , lon_v, lat_v, z_r_v                           &
     &                  , vll_v(:,:,:) )
        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , water_u(:,:,:,1),  -5.0d0, 5.0d0              &
     &                  , N_xi_v, N_eta_v, N_s_rho                      &
     &                  , lon_v, lat_v, z_r_v                           &
     &                  , ull_v(:,:,:) )
        do i=1,N_xi_v
          do j=1,N_eta_v
            do k=1,N_s_rho
              v(i,j,k) =-ull_v(i,j,k)*sinAv(i,j)                        &
                        +vll_v(i,j,k)*cosAv(i,j)
            enddo
          enddo
        enddo
          
          write(*,*) '*** SUCCESS Linear Interporation'
          
        else if (mode ==2) then
        
        end if
!        write(*,*) N_xi_rho, N_eta_rho,N_xi_u, N_eta_u, N_xi_v, N_eta_v

!        open(50, file='check3.dat')
!        do i=1, N_eta_v
!          write(50,*) v(:,i,14)
!        enddo
!        close(50)
!        open(50, file='check4.dat')
!        do i=1, Jm
!          write(50,*) water_u(:,i,20,1)
!        enddo
!        close(50)

!---- Convert temperature to potential temperature --------------------------------

        do i=1,N_xi_u
          do j=1,N_eta_u
            do k=1,N_s_rho
              CALL pzcon(1,grav,dynz,press,z_r(i,j,k))
              CALL POTMP(press,temp(i,j,k),salt(i,j,k),ref_press,potemp)
              temp(i,j,k) = potemp
            enddo
          enddo
        enddo

! --------------------------------------------------------

        ubar(:,:)=0.0d0
        do i=1,N_xi_u
          do j=1,N_eta_u
            do k=1,N_s_rho
              ubar(i,j)= ubar(i,j)+u(i,j,k)*abs(Cs_w(k-1)-Cs_w(k))
            enddo
          enddo
        enddo
        vbar(:,:)=0.0d0
        do i=1,N_xi_v
          do j=1,N_eta_v
            do k=1,N_s_rho
              vbar(i,j)= vbar(i,j)+v(i,j,k)*abs(Cs_w(k-1)-Cs_w(k))
            enddo
          enddo
        enddo
        
        if ( write_ini ) then
        
!---- Write ROMS initial conditions netCDF file --------------------------------

          start1D = (/ itime /)
          count1D = (/ 1 /)
          call writeNetCDF_1d(           &
!              input parameters
     &            'ocean_time'           &
     &          , INI_FILE               &
     &          , 1                      &
     &          , ocean_time             &
     &          , start1D, count1D       &
     &          )

     
          start3D = (/ 1,  1,  itime /)
          count3D = (/ N_xi_rho, N_eta_rho, 1 /)
          call writeNetCDF_3d(                      &
!              input parameters
     &            'zeta'                            &
     &          , INI_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , zeta                              &
     &          , start3D, count3D                  &
     &          )

          start3D = (/ 1,  1,  itime /)
          count3D = (/ N_xi_u, N_eta_u, 1 /)
          call writeNetCDF_3d(                      &
!              input parameters
     &            'ubar'                            &
     &          , INI_FILE                          &
     &          , N_xi_u, N_eta_u, 1                &
     &          , ubar                              &
     &          , start3D, count3D                  &
     &          )
          start3D = (/ 1,  1,  itime /)
          count3D = (/ N_xi_v, N_eta_v, 1 /)
          call writeNetCDF_3d(                      &
!              input parameters
     &            'vbar'                            &
     &          , INI_FILE                          &
     &          , N_xi_v, N_eta_v, 1                &
     &          , vbar                              &
     &          , start3D, count3D                  &
     &          )
              
          start4D = (/ 1,  1,  1,  itime /)
          count4D = (/ N_xi_rho, N_eta_rho, N_s_rho, 1 /)
          call writeNetCDF_4d(                  &
!              input parameters
     &            'temp'                            &
     &          , INI_FILE                          &
     &          , N_xi_rho, N_eta_rho, N_s_rho, 1   &
     &          , temp                              &
     &          , start4D, count4D                  &
     &          )

          call writeNetCDF_4d(                      &
!              input parameters
     &            'salt'                            &
     &          , INI_FILE                          &
     &          , N_xi_rho, N_eta_rho, N_s_rho, 1   &
     &          , salt                              &
     &          , start4D, count4D                  &
     &          )

          start4D = (/ 1,  1,  1,  itime /)
          count4D = (/ N_xi_u, N_eta_u, N_s_rho, 1 /)
          call writeNetCDF_4d(                      &
!              input parameters
     &            'u'                               &
     &          , INI_FILE                          &
     &          , N_xi_u, N_eta_u, N_s_rho, 1       &
     &          , u                                 &
     &          , start4D, count4D                  &
     &          )

          start4D = (/ 1,  1,  1,  itime /)
          count4D = (/ N_xi_v, N_eta_v, N_s_rho, 1 /)
          call writeNetCDF_4d(                      &
!              input parameters
     &            'v'                               &
     &          , INI_FILE                          &
     &          , N_xi_v, N_eta_v, N_s_rho, 1       &
     &          , v                                 &
     &          , start4D, count4D                  &
     &          )
     
        end if


!---- Write ROMS boundary conditions netCDF file --------------------------------

        start1D = (/ itime /)
        count1D = (/ 1 /)
        call writeNetCDF_1d(           &
!            input parameters
     &          'bry_time'             &
     &        , BRY_FILE               &
     &        , 1                      &
     &        , ocean_time             &
     &        , start1D, count1D       &
     &        )

        start2D = (/  1,        itime /)
        count2D = (/ N_eta_rho, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'zeta_west'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, 1                  &
     &        , zeta(1,:)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'zeta_east'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, 1                  &
     &        , zeta(N_xi_rho,:)              &
     &        , start2D, count2D              &
     &        )
        start2D = (/  1,       itime /)
        count2D = (/ N_xi_rho, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'zeta_south'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, 1                   &
     &        , zeta(:,1)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'zeta_north'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, 1                   &
     &        , zeta(:,N_eta_rho)             &
     &        , start2D, count2D              &
     &        )
     
        start2D = (/  1,      itime /)
        count2D = (/ N_eta_u, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'ubar_west'                   &
     &        , BRY_FILE                      &
     &        , N_eta_u, 1                    &
     &        , ubar(1,:)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'ubar_east'                   &
     &        , BRY_FILE                      &
     &        , N_eta_u, 1                    &
     &        , ubar(N_xi_u,:)                &
     &        , start2D, count2D              &
     &        )
        start2D = (/  1,     itime /)
        count2D = (/ N_xi_u, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'ubar_south'                  &
     &        , BRY_FILE                      &
     &        , N_xi_u, 1                     &
     &        , ubar(:,1)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'ubar_north'                  &
     &        , BRY_FILE                      &
     &        , N_xi_u, 1                     &
     &        , ubar(:,N_eta_u)               &
     &        , start2D, count2D              &
     &        )
     
        start2D = (/  1,      itime /)
        count2D = (/ N_eta_v, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'vbar_west'                   &
     &        , BRY_FILE                      &
     &        , N_eta_v, 1                    &
     &        , vbar(1,:)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'vbar_east'                   &
     &        , BRY_FILE                      &
     &        , N_eta_v, 1                    &
     &        , vbar(N_xi_v,:)                &
     &        , start2D, count2D              &
     &        )
        start2D = (/  1,     itime /)
        count2D = (/ N_xi_v, 1 /)
        call writeNetCDF_2d(                  &
!            input parameters
     &          'vbar_south'                  &
     &        , BRY_FILE                      &
     &        , N_xi_v, 1                     &
     &        , vbar(:,1)                     &
     &        , start2D, count2D              &
     &        )
        call writeNetCDF_2d(                  &
!            input parameters
     &          'vbar_north'                  &
     &        , BRY_FILE                      &
     &        , N_xi_v, 1                     &
     &        , vbar(:,N_eta_v)               &
     &        , start2D, count2D              &
     &        )
     
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_eta_rho, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'temp_west'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, N_s_rho, 1         &
     &        , temp(1,:,:)                   &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'temp_east'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, N_s_rho, 1         &
     &        , temp(N_xi_rho,:,:)            &
     &        , start3D, count3D              &
     &        )
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_xi_rho, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'temp_south'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, N_s_rho, 1          &
     &        , temp(:,1,:)                   &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'temp_north'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, N_s_rho, 1          &
     &        , temp(:,N_eta_rho,:)           &
     &        , start3D, count3D              &
     &        )
     
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_eta_rho, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'salt_west'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, N_s_rho, 1         &
     &        , salt(1,:,:)                   &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'salt_east'                   &
     &        , BRY_FILE                      &
     &        , N_eta_rho, N_s_rho, 1         &
     &        , salt(N_xi_rho,:,:)            &
     &        , start3D, count3D              &
     &        )
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_xi_rho, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'salt_south'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, N_s_rho, 1          &
     &        , salt(:,1,:)                   &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'salt_north'                  &
     &        , BRY_FILE                      &
     &        , N_xi_rho, N_s_rho, 1          &
     &        , salt(:,N_eta_rho,:)           &
     &        , start3D, count3D              &
     &        )
     
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_eta_u, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'u_west'                      &
     &        , BRY_FILE                      &
     &        , N_eta_u, N_s_rho, 1           &
     &        , u(1,:,:)                      &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'u_east'                      &
     &        , BRY_FILE                      &
     &        , N_eta_u, N_s_rho, 1           &
     &        , u(N_xi_u,:,:)                 &
     &        , start3D, count3D              &
     &        )
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_xi_u, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'u_south'                     &
     &        , BRY_FILE                      &
     &        , N_xi_u, N_s_rho, 1            &
     &        , u(:,1,:)                      &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'u_north'                     &
     &        , BRY_FILE                      &
     &        , N_xi_u, N_s_rho, 1            &
     &        , u(:,N_eta_u,:)             &
     &        , start3D, count3D              &
     &        )

        start3D = (/  1,  1,  itime /)
        count3D = (/ N_eta_v, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'v_west'                      &
     &        , BRY_FILE                      &
     &        , N_eta_v, N_s_rho, 1           &
     &        , v(1,:,:)                      &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'v_east'                      &
     &        , BRY_FILE                      &
     &        , N_eta_v, N_s_rho, 1           &
     &        , v(N_xi_v,:,:)                 &
     &        , start3D, count3D              &
     &        )
        start3D = (/  1,  1,  itime /)
        count3D = (/ N_xi_v, N_s_rho, 1 /)
        call writeNetCDF_3d(                  &
!            input parameters
     &          'v_south'                     &
     &        , BRY_FILE                      &
     &        , N_xi_v, N_s_rho, 1            &
     &        , v(:,1,:)                      &
     &        , start3D, count3D              &
     &        )
        call writeNetCDF_3d(                  &
!            input parameters
     &          'v_north'                     &
     &        , BRY_FILE                      &
     &        , N_xi_v, N_s_rho, 1            &
     &        , v(:,N_eta_v,:)                &
     &        , start3D, count3D              &
     &        )
     
        write(*,*) "*********************************************"
      
      enddo
      
!==== LOOP END ======================================================

      write(*,*) 'FINISH!!'
      

!**** End of Main program **********************************************
        
    CONTAINS
    
!**** create initial conditions NetCDF file **********************************************

      SUBROUTINE createNetCDFini(   &
!        input parameters
     &      OUT_FILE             &
     &    , TIME_ATT             &  
     &    , Im, Jm, Nz, Nt       &   
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Im, Jm, Nz, Nt
      
      integer :: ncid2,var_id2
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: ocean_time_dimid
      integer :: dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid2, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid2, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid2, 'ocean_time', NF90_UNLIMITED, ocean_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid2, 'spherical', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid2, 'Vtransform', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid2, 'Vstretching', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid2, 'theta_s', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid2, 'theta_b', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid2, 'Tcline', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid2, 'hc', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )


      call check( nf90_def_var(ncid2, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 's_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid2, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid2, 'ocean_time', NF90_DOUBLE, ocean_time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      dim3Dids = (/ xi_rho_dimid, eta_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'zeta', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_u_dimid, eta_u_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'ubar', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertically integrated u-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_v_dimid, eta_v_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'vbar', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertically integrated v-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_rho_dimid, eta_rho_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'temp', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid2, 'salt', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'psu') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_u_dimid, eta_u_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'u', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_v_dimid, eta_v_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'v', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )

  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFini
      
!**** create boundary conditions NetCDF file **********************************************

      SUBROUTINE createNetCDFbry(   &
!        input parameters
     &      OUT_FILE             &
     &    , TIME_ATT             &  
     &    , Im, Jm, Nz, Nt       &   
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Im, Jm, Nz, Nt
      
      integer :: ncid2,var_id2
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: bry_time_dimid
      integer :: dim2Dids(2), dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid2, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid2, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid2, 'bry_time', NF90_UNLIMITED, bry_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid2, 'spherical', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid2, 'Vtransform', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid2, 'Vstretching', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid2, 'theta_s', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid2, 'theta_b', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid2, 'Tcline', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid2, 'hc', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )


      call check( nf90_def_var(ncid2, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 's_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid2, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid2, 'bry_time', NF90_DOUBLE, bry_time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'open boundary conditions time') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      dim2Dids = (/ eta_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'zeta_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'zeta_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'zeta_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'zeta_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'ubar_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'ubar_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'ubar_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'ubar_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'vbar_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'vbar_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'vbar_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'vbar_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'u_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'u_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'u_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'u_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'v_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'v_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'v_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'v_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'temp_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'temp_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'temp_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'temp_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'salt_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'salt_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'salt_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'salt_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )


  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFbry
      
!**** writeNetCDF_1d **********************************************
      
      SUBROUTINE writeNetCDF_1d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im                     &
     &    , data                   &
     &    , start1D, count1D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im 
      real(8), intent( in) :: data(Im )
      integer, intent( in) :: start1D(1), count1D(1)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start1D, count = count1D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_1d
      
!**** writeNetCDF_2d **********************************************
      
      SUBROUTINE writeNetCDF_2d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm                 &
     &    , data                   &
     &    , start2D, count2D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm
      real(8), intent( in) :: data(Im, Jm)
      integer, intent( in) :: start2D(2), count2D(2)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start2D, count = count2D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_2d
      
!**** writeNetCDF_3d **********************************************
      
      SUBROUTINE writeNetCDF_3d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nt             &
     &    , data                   &
     &    , start3D, count3D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nt 
      real(8), intent( in) :: data(Im, Jm, Nt )
      integer, intent( in) :: start3D(3), count3D(3)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_3d
      
!**** writeNetCDF_4d **********************************************
      
      SUBROUTINE writeNetCDF_4d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nz, Nt         &
     &    , data                   &
     &    , start4D, count4D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nz, Nt
      real(8), intent( in) :: data(Im, Jm, Nz, Nt)
      integer, intent( in) :: start4D(4), count4D(4)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start4D, count = count4D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_4d
      
!**** readNetCDF_3d **********************************************
      
      SUBROUTINE readNetCDF_3d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nt             &
     &    , start3D, count3D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nt
      integer, intent( in) :: start3D(3), count3D(3)
      real(8), intent(out) :: data(Im, Jm, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start3D, count=count3D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do

      data(:,:,:)=data(:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_3d
      
!**** readNetCDF_4d **********************************************
      
      SUBROUTINE readNetCDF_4d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start4D, count4D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start4D(4), count4D(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start4D, count=count4D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do
      
      data(:,:,:,:)=data(:,:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_4d
      
!**** NetCDF utility **********************************************
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension

! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          err_flag = 1
!          stop "Stopped"
      end if
      
      END SUBROUTINE check2
      
! -------------------------------------------------------------------------
      
    END PROGRAM bryHYCOM2ROMS
      
