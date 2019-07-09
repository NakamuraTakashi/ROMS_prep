
!!!=== Copyright (c) 2014-2019 Takashi NAKAMURA  =====

    PROGRAM iniHYCOMpNAO2ROMS
      use netcdf
      use mod_roms_netcdf
      use mod_calendar
      use mod_interpolation
     
      implicit none
      
! -------------------------------------------------------------------------
      integer :: Syear, Smonth, Sday
      integer :: Ryear, Rmonth, Rday
      integer :: itime
      integer :: mode
      character(256) :: GRID_FILE
      character(256) :: HYCOM_prefix
      character(256) :: INI_prefix
      integer :: N_s_rho
      integer :: spherical
      integer :: Vtransform, Vstretching
      real(8) :: THETA_S , THETA_B, TCLINE, DCRIT

! ----- Configulation for NAOTIDE ---------------------------------------
! -----< Set mode, location, epoch, data interval, and output file >-----
!      
! Mode selection
      integer, parameter :: itmode   = 2
!      Set itmode = 1 to calculate geocentric tide
!                 = 2 to calculate pure ocean tide 
!                        with respect to ocean floor
!                 = 3 to calculate radial loading tide 
      integer, parameter :: lpmode  = 1
!      Set lpmode = 1 to use long-period ocean tide map of Takanezawa
!                 = 2 to use equilibrium tide (Valid for itmode = 1,2)
!
      
! -------------------------------------------------------------------------

      character(33) :: TIME_ATT  = "seconds since 2000-01-01 00:00:00"

      character(10) :: HYCOM_suffix = "_200001.nc"
      character(256) :: HYCOM_FILE
      character(15) :: INI_suffix   = "_20000101.00.nc"
      character(256) :: INI_FILE
      
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
      real(8), allocatable :: mask_rho(:,:)        ! land mask
      real(8), allocatable :: h_u(:,:)        ! depth (meter)
      real(8), allocatable :: h_v(:,:)        ! depth (meter)
      real(8), allocatable :: z_w(:,:,:)  ! vertical depths z_w(x,y,0:N(ng)) at  W-points (top/bottom cell) (meter)
      real(8), allocatable :: z_r(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
      real(8), allocatable :: z_r_u(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
      real(8), allocatable :: z_r_v(:,:,:)  ! vertical depths z_w(x,y,1:N(ng)) at  RHO-points (cell center) (meter)
     
      real(8) :: ocean_time(3)            ! Ocean time
      real(8), allocatable :: zeta(:,:,:)   ! free-surface (meter)
      real(8), allocatable :: ubar(:,:,:)   ! vertically integrated u-momentum component (meter second-1)
      real(8), allocatable :: vbar(:,:,:)   ! vertically integrated v-momentum component (milibar=hPa)
      real(8), allocatable :: u(:,:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: v(:,:,:,:)    ! v-momentum component (meter second-1)
      real(8), allocatable :: temp(:,:,:,:) ! potential temperature (Celsius)
      real(8), allocatable :: salt(:,:,:,:) ! salinity (psu)    

      real(8), allocatable :: cosAu(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: sinAu(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: cosAv(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: sinAv(:,:)  ! angle differece between HYCOM and ROMS coordinates (radian)
      real(8), allocatable :: ull_u(:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: vll_u(:,:,:)    ! v-momentum component (meter second-1)
      real(8), allocatable :: ull_v(:,:,:)    ! u-momentum component (meter second-1)
      real(8), allocatable :: vll_v(:,:,:)    ! v-momentum component (meter second-1)

      real(8), allocatable :: zeta_tide(:,:)   ! free-surface by tide (meter)

      
      real(8) :: hc       
      real(8), allocatable :: sc_w(:)       
      real(8), allocatable :: sc_r(:)  
      real(8), allocatable :: Cs_w(:)       
      real(8), allocatable :: Cs_r(:)  
      real(8) :: d_lat, d_lon
      
      integer :: i,j,k
      integer :: idt,incdf
      real(8) :: wf  
      real(8) :: Smjd
      
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      character(11) :: YYYYMMDDpHH
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_xi_u, N_eta_u
      integer :: N_xi_v, N_eta_v
      
      integer :: ncid,var_id
      integer :: ncid2,var_id2
      
      real(8) :: sf, off
      integer :: Im, Jm, Nz, Nt
      integer :: IL, IR, JB, JT
      integer :: ItS, ItE

			real(8) :: height, hsp, hlp, jtime
      real(8) :: ramp
      Logical :: Ldata
      !
      namelist/grd/GRID_FILE
      namelist/sdate/Syear, Smonth, Sday
      namelist/refdate/Ryear, Rmonth, Rday
      namelist/intpmode/mode
      namelist/hycom/HYCOM_prefix
      namelist/ini/INI_prefix
      namelist/ini/itime
      namelist/hcoord/spherical
      namelist/zcoord/N_s_rho
      namelist/zcoord/Vtransform, Vstretching
      namelist/zcoord/THETA_S, THETA_B, TCLINE, DCRIT

      ! Read parameters in namelist file
      
      read (*, nml=grd)
      read (*, nml=sdate)
      read (*, nml=refdate)
      read (*, nml=intpmode)
      read (*, nml=hycom)
      read (*, nml=ini)
      read (*, nml=hcoord)
      read (*, nml=zcoord)
      
      !---- HYCOM file name --------------------------------

      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth

      HYCOM_suffix(2:5)=YYYY
      HYCOM_suffix(6:7)=MM
      HYCOM_FILE = trim( HYCOM_prefix )//HYCOM_suffix
      
      call mjdymd(Smjd, Ryear, Rmonth, Rday , 0,    &
     &            0, 0     , 1                      )  ! from naotidej.f

      
      write (YYYY, "(I4.4)") Ryear
      write (MM, "(I2.2)") Rmonth
      write (DD, "(I2.2)") Rday
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(15:18)=YYYY
      TIME_ATT(20:21)=MM
      TIME_ATT(23:24)=DD
      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", trim( GRID_FILE )
      
      ! Open NetCDF grid file
      call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
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
      allocate(mask_rho(N_xi_rho, N_eta_rho))
      allocate(h_u(N_xi_u, N_eta_u))
      allocate(h_v(N_xi_v, N_eta_v))
      allocate(z_w(N_xi_rho, N_eta_rho, 0:N_s_rho))
      allocate(z_r(N_xi_rho, N_eta_rho, 1:N_s_rho))
      allocate(z_r_u(N_xi_u, N_eta_u, 1:N_s_rho))
      allocate(z_r_v(N_xi_v, N_eta_v, 1:N_s_rho))
      
      allocate(zeta(N_xi_rho, N_eta_rho, 3))
      allocate(ubar(N_xi_u, N_eta_u, 3))
      allocate(vbar(N_xi_v, N_eta_v, 3))
      allocate(u(N_xi_u, N_eta_u, N_s_rho, 3))
      allocate(v(N_xi_v, N_eta_v, N_s_rho, 3))
      allocate(temp(N_xi_rho, N_eta_rho, N_s_rho, 3))
      allocate(salt(N_xi_rho, N_eta_rho, N_s_rho, 3))
      
      allocate(cosAu(N_xi_u, N_eta_u))
      allocate(sinAu(N_xi_u, N_eta_u))
      allocate(cosAv(N_xi_v, N_eta_v))
      allocate(sinAv(N_xi_v, N_eta_v))
      allocate(ull_u(N_xi_u, N_eta_u, N_s_rho))
      allocate(vll_u(N_xi_u, N_eta_u, N_s_rho))
      allocate(vll_v(N_xi_v, N_eta_v, N_s_rho))
      allocate(ull_v(N_xi_v, N_eta_v, N_s_rho))
      
      allocate(zeta_tide(N_xi_rho, N_eta_rho))
      
      allocate( sc_w(0:N_s_rho) )
      allocate( sc_r(1:N_s_rho) )
      allocate( Cs_w(0:N_s_rho) )
      allocate( Cs_r(1:N_s_rho) )
      
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
      call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
      call check( nf90_get_var(ncid, var_id, mask_rho) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", trim( GRID_FILE )
      
      
      
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

      call set_depth (              &
     &              N_xi_rho, N_eta_rho, N_s_rho       &
     &            , h               &
     &            , Vtransform      &
     &            , hc              &
     &            , sc_w, sc_r      &
     &            , Cs_w, Cs_r      &
!        output parameters
     &            , z_r, z_w)
      z_r = -z_r
      z_w = -z_w
      
      do i=1,N_xi_u
        do j=1,N_eta_u
          do k=1,N_s_rho
            z_r_u(i,j,k)= (z_r(i,j,k)+z_r(i+1,j,k))*0.5d0
          enddo
        enddo
      enddo
      do i=1,N_xi_v
        do j=1,N_eta_v
          do k=1,N_s_rho
            z_r_v(i,j,k)= (z_r(i,j,k)+z_r(i,j+1,k))*0.5d0
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
      
!---- Read extracted HYCOM netCDF file --------------------------------

      ! Open NetCDF file
      write(*,*) "OPEN: ", trim( HYCOM_FILE )
      call check( nf90_open(trim( HYCOM_FILE ), nf90_nowrite, ncid) )
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
      write(*,*) "CLOSE: ", trim( HYCOM_FILE )
      write(*,*) "*********************************************"
      
!---- Read extracted HYCOM netCDF file --------------------------------

      CALL cdate2(time_all(itime)/24.0d0,Ryear,Rmonth,Rday,YYYYMMDDpHH)
      write(*,*) 'time = ', YYYYMMDDpHH

      write(*,*) "OPEN: ", trim( HYCOM_FILE )
      call check( nf90_open(trim( HYCOM_FILE ), nf90_nowrite, ncid) )
      
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
      write(*,*) "CLOSE: ", trim( HYCOM_FILE )
      
! --------------------------------------------------------

      ocean_time(1) = time_all(itime)*60.0d0*60.0d0  !! hour -> second

! --- Linear interporation --------------------------------
      
      if (mode==1) then
      
        write(*,*) 'Linear Interporation: zeta'
        call LinearInterpolation2D_grid2(Im, Jm, lon, lat               &
     &                  , surf_el(:,:,1), -5.0d0, 5.0d0                 &
     &                  , N_xi_rho, N_eta_rho, lon_rho, lat_rho, zeta(:,:,1) )
      
        write(*,*) 'Linear Interporation: temp'
        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , water_temp(:,:,:,1),  0.0d0, 40.0d0           &
     &                  , N_xi_rho, N_eta_rho, N_s_rho                  &
     &                  , lon_rho, lat_rho, z_r                         &
     &                  , temp(:,:,:,1) )
     
        write(*,*) 'Linear Interporation: salt'
        call LinearInterpolation3D_grid3(Im, Jm, Nz, lon, lat, depth    &
     &                  , salinity(:,:,:,1),  20.0d0, 36.0d0            &
     &                  , N_xi_rho, N_eta_rho, N_s_rho                  &
     &                  , lon_rho, lat_rho, z_r                         &
     &                  , salt(:,:,:,1) )

        write(*,*) 'Linear Interporation: u'
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
              u(i,j,k,1) = ull_u(i,j,k)*cosAu(i,j)                      &
                          +vll_u(i,j,k)*sinAu(i,j)
            enddo
          enddo
        enddo
        

        write(*,*) 'Linear Interporation: v'
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
              v(i,j,k,1) =-ull_v(i,j,k)*sinAv(i,j)                      &
                          +vll_v(i,j,k)*cosAv(i,j)
            enddo
          enddo
        enddo

        write(*,*) '*** SUCCESS Linear Interporation'
        
      else if (mode ==2) then
      
      end if

!---- Convert temperature to potential temperature --------------------------------

      do i=1,N_xi_u
        do j=1,N_eta_u
          do k=1,N_s_rho
            CALL pzcon(1,grav,dynz,press,z_r(i,j,k))
            CALL POTMP(press,temp(i,j,k,1),salt(i,j,k,1),ref_press,potemp)
            temp(i,j,k,1) = potemp
          enddo
        enddo
      enddo
      
          
#if defined NAOTIDE || defined NAOTIDEJ

!---- NAOTIDE calculation --------------------------------
          
      jtime = dble(Smjd) + ocean_time(1)/86400.0d0
!      write(*,*) jtime, Smjd
      do i=1,N_xi_rho
        do j=1,N_eta_rho
# if defined NAOTIDE
          call naotide (lon_rho(i,j), lat_rho(i,j), jtime  , itmode, lpmode,       &
     &               zeta_tide(i,j), hsp   , hlp   , Ldata          )
# elif defined NAOTIDEJ
          call naotidej(lon_rho(i,j), lat_rho(i,j), jtime  , itmode, lpmode,       &
     &               zeta_tide(i,j), hsp   , hlp   , Ldata          )
# endif          
          zeta_tide(i,j)=zeta_tide(i,j)*0.01d0

        enddo
      enddo
      
      call FillExtraPoints2D(N_xi_rho, N_eta_rho, zeta_tide, -5.0d0, 5.0d0)

      zeta(:,:,1)=zeta(:,:,1)+zeta_tide(:,:)
#endif
         
!---- Depth averaged velocity calculation --------------------------------
      
      ubar(:,:,1)=0.0d0
      do i=1,N_xi_u
        do j=1,N_eta_u
          do k=1,N_s_rho
            ubar(i,j,1)= ubar(i,j,1)+u(i,j,k,1)*abs(Cs_w(k-1)-Cs_w(k))
          enddo
        enddo
      enddo
      vbar(:,:,1)=0.0d0
      do i=1,N_xi_v
        do j=1,N_eta_v
          do k=1,N_s_rho
            vbar(i,j,1)= vbar(i,j,1)+v(i,j,k,1)*abs(Cs_w(k-1)-Cs_w(k))
          enddo
        enddo
      enddo
!
!---- Initial netcdf file name -------------------------------------------------
      INI_suffix(2:12)=YYYYMMDDpHH
      INI_FILE = trim( INI_prefix )//INI_suffix

!---- Create the ROMS initial conditions netCDF file --------------------------------
          
      call createNetCDFini(                            &
!          input parameters
     &        trim( INI_FILE )                         &
     &      , TIME_ATT                                 &
     &      , N_xi_rho, N_eta_rho, N_s_rho, 1          &
     &      )
      
      call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid) )
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
!          input parameters
     &        's_rho'                &
     &      , trim( INI_FILE )       &
     &      , N_s_rho                &
     &      , sc_r                   &
     &      , start1D, count1D       &
     &      )
      call writeNetCDF_1d(           &
!          input parameters
     &        'Cs_r'                 &
     &      , trim( INI_FILE )       &
     &      , N_s_rho                &
     &      , Cs_r                   &
     &      , start1D, count1D       &
     &      )
     
      start1D = (/ 1 /)
      count1D = (/ N_s_rho+1 /)
      call writeNetCDF_1d(           &
!          input parameters
     &        's_w'                  &
     &      , trim( INI_FILE )       &
     &      , N_s_rho+1              &
     &      , sc_w                   &
     &      , start1D, count1D       &
     &      )
      call writeNetCDF_1d(           &
!          input parameters
     &        'Cs_w'                 &
     &      , trim( INI_FILE )       &
     &      , N_s_rho+1              &
     &      , Cs_w                   &
     &      , start1D, count1D       &
     &      )
      
!---- Write ROMS initial conditions netCDF file --------------------------------

      start1D = (/ 1 /)
      count1D = (/ 1 /)
      call writeNetCDF_1d(                      &
!          input parameters
     &        'ocean_time'                      &
     &      , trim( INI_FILE )                  &
     &      , 1                                 &
     &      , ocean_time(1)                     &
     &      , start1D, count1D                  &
     &      )

     
      start3D = (/ 1,  1,  1 /)
      count3D = (/ N_xi_rho, N_eta_rho, 1 /)
      call writeNetCDF_3d(                      &
!          input parameters
     &        'zeta'                            &
     &      , trim( INI_FILE )                  &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , zeta(:,:,1)                       &
     &      , start3D, count3D                  &
     &      )

      start3D = (/ 1,  1,  1 /)
      count3D = (/ N_xi_u, N_eta_u, 1 /)
      call writeNetCDF_3d(                      &
!          input parameters
     &        'ubar'                            &
     &      , trim( INI_FILE )                  &
     &      , N_xi_u, N_eta_u, 1                &
     &      , ubar(:,:,1)                       &
     &      , start3D, count3D                  &
     &      )
      start3D = (/ 1,  1,  1 /)
      count3D = (/ N_xi_v, N_eta_v, 1 /)
      call writeNetCDF_3d(                      &
!          input parameters
     &        'vbar'                            &
     &      , trim( INI_FILE )                  &
     &      , N_xi_v, N_eta_v, 1                &
     &      , vbar(:,:,1)                       &
     &      , start3D, count3D                  &
     &      )
          
      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ N_xi_rho, N_eta_rho, N_s_rho, 1 /)
      call writeNetCDF_4d(                  &
!          input parameters
     &        'temp'                            &
     &      , trim( INI_FILE )                  &
     &      , N_xi_rho, N_eta_rho, N_s_rho, 1   &
     &      , temp(:,:,:,1)                     &
     &      , start4D, count4D                  &
     &      )

      call writeNetCDF_4d(                      &
!          input parameters
     &        'salt'                            &
     &      , trim( INI_FILE )                  &
     &      , N_xi_rho, N_eta_rho, N_s_rho, 1   &
     &      , salt(:,:,:,1)                     &
     &      , start4D, count4D                  &
     &      )

      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ N_xi_u, N_eta_u, N_s_rho, 1 /)
      call writeNetCDF_4d(                      &
!          input parameters
     &        'u'                               &
     &      , trim( INI_FILE )                  &
     &      , N_xi_u, N_eta_u, N_s_rho, 1       &
     &      , u(:,:,:,1)                        &
     &      , start4D, count4D                  &
     &      )

      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ N_xi_v, N_eta_v, N_s_rho, 1 /)
      call writeNetCDF_4d(                      &
!          input parameters
     &        'v'                               &
     &      , trim( INI_FILE )                  &
     &      , N_xi_v, N_eta_v, N_s_rho, 1       &
     &      , v(:,:,:,1)                        &
     &      , start4D, count4D                  &
     &      )

      write(*,*) 'FINISH!!'
      
    END PROGRAM iniHYCOMpNAO2ROMS
      
