
!!!=== Copyright (c) 2017 Takashi NAKAMURA  =====

    PROGRAM frcJRA55
      use netcdf
      use wgrib2api
      use mod_interpolation
      use mod_calendar
     
      implicit none
      
! SETTINGS of JMA MSM Data reading area  -----------------------------------------------------------
! -------------------------------------------------------------------------
! NetCDF file for the box corners bounded by (Llon,Blat) and (Rlon,Tlat).
!
!  　　　　　　　 ______ (Rlon,Tlat)
!                |      |
!                |      |
!                |______|
!     (Llon,Blat)                     
!
! On Input:
!
!    Llon         Box left-edge   longitude (degrees, -180 - 180)
!    Rlon         Box right-edge  longitude (degrees, -180 - 180)
!    Blat         Box bottom-edge latitude  (degress,  -90 - 90 )
!    Tlat         Box top-edge    latitude  (degress,  -90 - 90 )
!
! Geographical and Grid parameters --------

      real(8), parameter :: Tlat = 27.5d0  ! 24.9d0   !27.5d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
      real(8), parameter :: Blat = 22.5d0  ! 23.7d0   !22.5d0    ! Latitude  (degrees) of the top-right corner of the grid.
      real(8), parameter :: Llon = 120.0d0 ! 123.5d0  !120.0d0   ! Longitude (degrees) of the bottom-left corner of the grid. 
      real(8), parameter :: Rlon = 127.0d0 ! 125.0d0  !127.0d0   ! Longitude (degrees) of the top-right corner of the grid. 

      integer, parameter :: Syear  = 2016   ! Starting year
      integer, parameter :: Smonth = 4      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
      
      integer, parameter :: Eyear  = 2016   ! Ending year   !!!4/30, 5/31, 6/30, 7/31, 8/31, 9/30
      integer, parameter :: Emonth = 5      ! Ending month
      integer, parameter :: Eday   = 1      ! Ending day
      
      integer, parameter :: Ryear  = 2000   ! Reference year
      integer, parameter :: Rmonth = 1      ! Reference month
      integer, parameter :: Rday   = 1      ! Reference day
      
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "D:/ROMS/Data/Coral_Triangle/CT_0.04_grd_v1.nc"

      character(len=*), parameter :: OUT_FILE  = "D:/ROMS/Data/Coral_Triangle/CT_frc_JRA55_1604.nc"

      character(len=*) :: G_TEMP_FILE = "D:/JRA-55/fcst_surf125.011_tmp.2016040100_2016043021"
      character(len=*) :: G_U10_FILE  = "D:/JRA-55/fcst_surf125.033_ugrd.2016040100_2016043021"
      character(len=*) :: G_V10_FILE  = "D:/JRA-55/fcst_surf125.034_vgrd.2016040100_2016043021"
      character(len=*) :: G_RH_FILE   = "D:/JRA-55/fcst_surf125.052_rh.2016040100_2016043021"
      character(len=*) :: G_PRES_FILE = "D:/JRA-55/fcst_surf125.001_pres.2016040100_2016043021"
      character(len=*) :: G_TCDC_FILE = "D:/JRA-55/fcst_surf125.071_tcdc.2016040100_2016043021"
      character(len=*) :: G_PREC_FILE = "D:/JRA-55/fcst_column125.054_pwat.2016040100_2016043021"
     
      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8) :: time2(1) ! Ocean time
      real(8), allocatable :: Uwind(:,:) ! surface u-wind component (meter second-1)
      real(8), allocatable :: Vwind(:,:) ! surface v-wind component (meter second-1)
      real(8), allocatable :: Pair (:,:) ! surface air pressure (milibar=hPa)
      real(8), allocatable :: Tair (:,:) ! surface air temperature (Celsius)
      real(8), allocatable :: Qair (:,:) ! surface air relative humidity (percentage)
      real(8), allocatable :: rain (:,:) ! rain fall rate (kilogram meter-2 second-1)
      real(8), allocatable :: cloud(:,:) ! cloud fraction (0 to 1)    
           
      real(8), allocatable :: lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:), time(:)
      real(8), allocatable :: u(:, :, :), v(:, :, :),temp(:, :, :)
      real(8), allocatable :: psea(:, :, :), rh(:, :, :),ncld(:, :, :)
      real(8), allocatable :: r1h(:, :, :)
      integer :: start1D(1), count1D(1)
      integer :: start3D(3), count3D(3)
      real(8), allocatable :: u2(:, :)
      
      integer :: N_days
      integer :: N_ref_days
      integer :: julian_date
      integer :: ref_julian_date
      integer :: iyear, imonth, iday
      integer :: i,j,k
      integer :: idays
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      real(8) :: sf, off
      integer :: Im, Jm, Nt
      integer :: IL, IR, JB, JT
      
      integer :: ncid2,var_id2
      integer :: start1Dw(1), count1Dw(1)
      integer :: start3Dw(3), count3Dw(3)
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
!      integer :: time_varid, Uwind_varid
      


      call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
      call ndays(Smonth, Sday, Syear, Rmonth, Rday, Ryear, N_ref_days)
      call jd(Syear, Smonth, Sday, julian_date)
      call jd(Ryear, Rmonth, Rday, ref_julian_date)
      
!      write(*,*) N_days
      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(Uwind(N_xi_rho, N_eta_rho))
      allocate(Vwind(N_xi_rho, N_eta_rho))
      allocate(Pair (N_xi_rho, N_eta_rho))
      allocate(Tair (N_xi_rho, N_eta_rho))
      allocate(Qair (N_xi_rho, N_eta_rho))
      allocate(rain (N_xi_rho, N_eta_rho))
      allocate(cloud(N_xi_rho, N_eta_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
!---- Modify time-unit description ---------------------------------
      
      call cdate(ref_julian_date, iyear, imonth, iday)
      
      write (YYYY, "(I4.4)") iyear
      write (MM, "(I2.2)") imonth
      write (DD, "(I2.2)") iday
      
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      
!---- Update file name  -----------------------------------------------
      
      call cdate(julian_date, iyear, imonth, iday)
      
      write (YYYY, "(I4.4)") iyear
      write (MM, "(I2.2)") imonth
      write (DD, "(I2.2)") iday
      
       
      NC_FILE(74:77)=YYYY
      NC_FILE(79:80)=MM
      NC_FILE(81:82)=DD
      
      NC_FILE2(78:81)=YYYY
      NC_FILE2(83:84)=MM
      NC_FILE2(85:86)=DD
      
      write(*,*) "OPEN: ", NC_FILE
            
!---- Read JRA-55 GRIB2 file --------------------------------

      !Open NetCDF file
      call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )

      ! Get dimension data
      call get_dimension(ncid, 'lat', Jm)
      call get_dimension(ncid, 'lon', Im)
      call get_dimension(ncid, 'time', Nt)
          
!      write(*,*) Im, Jm, Nt
      
!      ! Allocate variable
      allocate(lat_all(Jm))
      allocate(lon_all(Im))
      
      call check( nf90_inq_varid(ncid, 'lat', var_id) )
      call check( nf90_get_var(ncid, var_id, lat_all) )
      call check( nf90_inq_varid(ncid, 'lon', var_id) )
      call check( nf90_get_var(ncid, var_id, lon_all) )
      
      JT=1
      do j=1,Jm
        if(lat_all(j)<Tlat) exit
        JT=j
      end do
      JB=JT
      do j=JT,Jm
        if(lat_all(j)<=Blat) exit
        JB=j
      end do
      IL=1
      do i=IL,Im
        if(lon_all(i)>=Llon) exit
        IL=i
      end do
      IR=IL
      do i=IL,Im
        if(lon_all(i)>Rlon) exit
        IR=i
      end do
      write(*,*) JB, lat_all(JB), JT, lat_all(JT)
      write(*,*) IL, lon_all(IL), IR, lon_all(IR)
      
      Jm=JB-JT-1
      Im=IR-IL-1
      
      ! Initialize count and start for reading data
      start3D = (/ IL, JT, 1 /)
      count3D = (/ Im, Jm, Nt /)
      
      allocate(lat(Jm))
      allocate(lon(Im))
      
      allocate(time(Nt))
      allocate(u(Im, Jm, Nt))
      allocate(u2(Im, Jm))
      allocate(v(Im, Jm, Nt))
      allocate(temp(Im, Jm, Nt))
      allocate(psea(Im, Jm, Nt))
      allocate(rh(Im, Jm, Nt))
      allocate(ncld(Im, Jm, Nt))
      allocate(r1h(Im, Jm, Nt))
      
      write(*,*) "*** SUCCESS check dimension!"
      
!---- Create the forcing netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'xi_rho', N_xi_rho, xi_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_rho',N_eta_rho, eta_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'time', NF90_UNLIMITED, time_dimid) )
!      call check( nf90_def_dim(ncid2, 'time', 48, time_dimid) )

      dimids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)

    ! Define the netCDF variables for the pressure and temperature data.
      call check( nf90_def_var(ncid2, 'time', NF90_DOUBLE, time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'atmospheric forcing time') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      call check( nf90_def_var(ncid2, 'Uwind', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'surface u-wind component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'Vwind', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'surface v-wind component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'Pair', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'surface air pressure') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'millibar') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'Tair', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'surface air temperature') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'Qair', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'surface air relative humidity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'percentage') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'rain', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'rain fall rate') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'kilogram meter-2 second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'cloud', NF90_REAL, dimids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'cloud fraction') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     '0 to 1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )


  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )




      julian_date = julian_date - 2
      
!---- LOOP START --------------------------------

      do idays=1, N_days+1

! --- Update  -----------------------------------------------

        julian_date = julian_date + 1
        
        call cdate(julian_date, iyear, imonth, iday)
      
        write (YYYY, "(I4.4)") iyear
        write (MM, "(I2.2)") imonth
        write (DD, "(I2.2)") iday
         
        NC_FILE(74:77)=YYYY
        NC_FILE(79:80)=MM
        NC_FILE(81:82)=DD
        
        NC_FILE2(78:81)=YYYY
        NC_FILE2(83:84)=MM
        NC_FILE2(85:86)=DD
      
! --- Read JMA MSM netCDF file --------------------------------
      
        write(*,*) "OPEN: ", NC_FILE
      
        !Open NetCDF file
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )

!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'time', var_id) )  !!!  since 2014-MM-DD 00:00:00+00:00, not Japan time (00:00:00+09:00)
        call check( nf90_get_var(ncid, var_id, time) )
        
        start1D(1)=JT
        count1D(1)=Jm
        call check( nf90_inq_varid(ncid, 'lat', var_id) )
        call check( nf90_get_var(ncid, var_id, lat, start=start1D, count=count1D) )
        
        start1D(1)=IL
        count1D(1)=Im
        call check( nf90_inq_varid(ncid, 'lon', var_id) )
        call check( nf90_get_var(ncid, var_id, lon, start=start1D, count=count1D) )
        
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'u', var_id) ) ! Eastward component of wind (m/s)
        call check( nf90_get_var(ncid, var_id, u, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        u(:,:,:)=u(:,:,:)*sf+off
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'v', var_id) ) ! Northward component of wind (m/s)
        call check( nf90_get_var(ncid, var_id, v, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        v(:,:,:)=v(:,:,:)*sf+off
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'temp', var_id) ) ! Temperature (K)
        call check( nf90_get_var(ncid, var_id, temp, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        temp(:,:,:)=temp(:,:,:)*sf+off
        temp(:,:,:)=temp(:,:,:)-273.15d0   !!! convert K -> oC
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'psea', var_id) )  ! Sea level pressure (Pa)
        call check( nf90_get_var(ncid, var_id, psea, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        psea(:,:,:)=psea(:,:,:)*sf+off
        psea(:,:,:)=psea(:,:,:)*1.0d-2     !!! convert Pa -> millibar (= hPa)
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'rh', var_id) ) ! Relative humidity (%)
        call check( nf90_get_var(ncid, var_id, rh, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        rh(:,:,:)=rh(:,:,:)*sf+off
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'ncld', var_id) ) ! Cloud amount (%)
        call check( nf90_get_var(ncid, var_id, ncld, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        ncld(:,:,:)=ncld(:,:,:)*sf+off
        ncld(:,:,:)=ncld(:,:,:)*1.0d-2     !!! convert % -> 0 to 1
       ! Close NetCDF file
        call check( nf90_close(ncid) )
        
        write(*,*) "OPEN: ", NC_FILE2
      
        !Open NetCDF file
        call check( nf90_open(NC_FILE2, nf90_nowrite, ncid) )
!        ! Get variable id
        call check( nf90_inq_varid(ncid, 'r1h', var_id) ) ! Rainfall in 1 hour (mm/h)
        call check( nf90_get_var(ncid, var_id, r1h, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        r1h(:,:,:)=r1h(:,:,:)*sf+off
        r1h(:,:,:)=r1h(:,:,:)*1.0d0/3600.0d0  !!! convert mm/h -> kg/m2/s
                                              !!! 1mm x 1m x 1m = 0.1x100x100 cm3 = 1000 cm3 = 1L = 1kg
       ! Close NetCDF file
        call check( nf90_close(ncid) )
        
        write(*,*) "*** SUCCESS reading JMA MSM netcdf file!"
        
        do k=1,Nt
          
          start1Dw = (/ k+(idays-1)*Nt /)
          count1Dw = (/ 1 /)
          
          time2(1) = dble(idays)-2.0d0+dble(k)/24.0d0 + dble(N_ref_days)
          
          write(*,*) 'time (days) = ', time2(1)
          
          call check( nf90_open(OUT_FILE, NF90_WRITE, ncid2) )
          call check( nf90_inq_varid(ncid2, "time", var_id2) )
          call check( nf90_put_var(ncid2, var_id2, time2, start = start1Dw, count = count1Dw) )
          call check( nf90_close(ncid2) )
          
          call interp_writeNetCDF(   &
!          input parameters
     &      'Uwind'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , u                    &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
          call interp_writeNetCDF(   &
!          input parameters
     &      'Vwind'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , v                    &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
          call interp_writeNetCDF(   &
!          input parameters
     &      'Pair'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , psea                 &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
          call interp_writeNetCDF(   &
!          input parameters
     &      'Tair'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , temp                 &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
          call interp_writeNetCDF(   &
!          input parameters
     &      'Qair'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , rh                   &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
          call interp_writeNetCDF(   &
!          input parameters
     &      'rain'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , r1h                  &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
          call interp_writeNetCDF(   &
!          input parameters
     &      'cloud'              &  !!!
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , ncld                 &  !!!
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
     
     
          write(*,*) '*************************************************'
        end do
        

      end do
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
        
    CONTAINS
    
      SUBROUTINE interp_writeNetCDF(  &
!        input parameters
     &      NCNAME               &
     &    , mode                 &  ! mode=1, kirnging; mode=2, linear intrtporation
     &    , OUT_FILE             &
     &    , Im, Jm, Nt           &   
     &    , N_xi_rho, N_eta_rho  &   
     &    , idays, k             &   
     &    , lat                  &  
     &    , lon                  &  
     &    , Ydata                &  
     &    , lat_rho              &  
     &    , lon_rho              &  
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: NCNAME
      integer, intent( in) :: mode           
      character(len=*),  intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nt            
      integer, intent( in) :: N_xi_rho, N_eta_rho
      integer, intent( in) :: idays, k
      real(8), intent( in) :: lat(Jm)      
      real(8), intent( in) :: lon(Im)      
      real(8), intent( in) :: Ydata(Im, Jm, Nt)      
      real(8), intent( in) :: lat_rho(N_xi_rho, N_eta_rho)      
      real(8), intent( in) :: lon_rho(N_xi_rho, N_eta_rho)      


      real(8) :: data(N_xi_rho, N_eta_rho)
      real(8) :: time(1)
      integer :: ncid,var_id
      integer :: start3Dw(3), count3Dw(3)
      
! --- Linear interporation --------------------------------
      if (mode==1) then
      
        call LinearInterpolation2D_grid(Im, Jm, lon, lat, Ydata(:,:,k)  &
     &                  , N_xi_rho, N_eta_rho, lon_rho, lat_rho, data )
        
        write(*,*) '*** SUCCESS Linear Interporation of ', NCNAME
        
      else if (mode ==2) then
      
      end if
      
      
! --- Write NetCDF bulk forcing file ------------------------
      
      start3Dw = (/ 1,        1,         k+(idays-1)*Nt /)
      count3Dw = (/ N_xi_rho, N_eta_rho, 1 /)
      

      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3Dw, count = count3Dw) )
      call check( nf90_close(ncid) )

      END SUBROUTINE interp_writeNetCDF
      
      
!---- NetCDF utility -------------------------------------------------
     
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension


      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
    END PROGRAM frcJRA55
      
