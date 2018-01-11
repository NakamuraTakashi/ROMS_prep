
!!!=== Copyright (c) 2018 Takashi NAKAMURA  =====

    PROGRAM frcDSJRA55
      use netcdf
      use eccodes
      use mod_roms_netcdf
      use mod_interpolation
      use mod_calendar
     
      implicit none
      
! SETTINGS of input and output files  -----------------------------------

      integer, parameter :: Syear  = 1996   ! Sampling year
      integer, parameter :: Smonth = 9      ! Sampling month

     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "D:/ROMS/Data/Yaeyama/Yaeyama1_grd_v10.nc"

      integer, parameter :: N_Param = 7
      character(len=*), parameter :: GRIB_FCST_SURF_dir  =  &
     &   "D:/DSJRA-55/Hist/Daily/fcst_surf/199609/"
      character(len=*), parameter :: GRIB_FCST_PHY2M_dir =  &
     &   "D:/DSJRA-55/Hist/Daily/fcst_phy2m/199609/"
     
      character(len=*), parameter :: GRIB_LL = "D:/DSJRA-55/Consts/Lambert5km_latlon.dat"
      
      character(len=*), parameter :: OUT_prefix  = "output/Yaeyama1_frc_DSJRA55_"
      
      integer, parameter :: Ryear  = 2000   ! Reference year for NetCDF
      integer, parameter :: Rmonth = 1      ! Reference month for NetCDF
      integer, parameter :: Rday   = 1      ! Reference day for NetCDF
      
!----------------------------------------------------------------------

      character(len=*), parameter :: GRIB_FCST_SURF_prefix  = "fcst_surf"
      character(len=*), parameter :: GRIB_FCST_PHY2M_prefix = "fcst_phy2m"
      character(11) :: GRIB_suffix  = ".1958010100"
      character(10) :: OUT_suffix   = "_195801.nc"
      
      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"

      character(256) :: GRIB_FILE
      character(256) :: OUT_FILE(N_Param)
      
      integer, parameter :: GRIB_NAME(3,N_Param) = reshape ((/   &
!     productDefinitionTemplateNumber, parameterCategory, parameterNumber
     &    0, 2, 2          &  ! 10 metre U wind component
     &  , 0, 2, 3          &  ! 10 metre V wind component
     &  , 0, 0, 0          &  ! air_temperature
     &  , 0, 0, 7          &  ! dew-point deficit (K)
     &  , 0, 3, 0          &  ! Surface pressure
     &  , 0, 6, 1          &  ! Total Cloud Cover
     &  , 8, 1, 52         &  ! Total precipitation rate  (kg m-2 s-1)
     &  /), (/3, N_Param/))
     
      character(256) :: NC_NAME(N_Param) = (/     &
     &   "Uwind"                                  &
     &  ,"Vwind"                                  &
     &  ,"Tair "                                  &
     &  ,"Qair "                                  &
     &  ,"Pair "                                  &
     &  ,"cloud"                                  &
     &  ,"rain "                                  &
     &  /)
      character(256) :: NC_LNAME(N_Param) = (/    &
     &   "surface u-wind component     "          &
     &  ,"surface v-wind component     "          &
     &  ,"surface air temperature      "          &
     &  ,"surface air relative humidity"          &
     &  ,"surface air pressure         "          &
     &  ,"cloud fraction               "          &
     &  ,"rain fall rate               "          &
     &  /)
      character(256) :: NC_UNIT(N_Param) = (/     &
     &   "meter second-1           "              &
     &  ,"meter second-1           "              &
     &  ,"Celsius                  "              &
     &  ,"percentage               "              &
     &  ,"millibar                 "              &
     &  ,"0 to 1                   "              &
     &  ,"kilogram meter-2 second-1"              &
     &  /)
     
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8), allocatable :: out_data(:,:) ! output forcing data
           
      integer :: Im, Jm
!      integer, parameter :: Im = 721
!      integer, parameter :: Jm = 577
      real(8), allocatable :: lat(:, :), lon(:, :)
!      real(8), allocatable :: grib_data(:, :)
      real(8), allocatable :: grib_data(:,:, :)
!      real(8) :: grib_data(Im,Jm)
      real(8) :: time(1)
      integer :: start1D(1), count1D(1)
      integer :: start3D(3), count3D(3)
      
      integer :: iyear, imonth, iday
      integer :: ihour, imin
      integer :: i,j,k
      integer :: idays
      integer :: itime
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      character(2) :: hh
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
!      integer :: Im, Jm
      
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
      
      integer :: iparam
      integer :: ifile,idx,iret,igrib
      integer :: istart, iend
      integer :: YYYYMMDD, hhmm
      real(8), allocatable :: values(:)
      real :: bd
      
      integer, allocatable :: tmp(:)
      integer :: itmp

!---- Modify time-unit description ---------------------------------
      
      write (YYYY, "(I4.4)") Ryear
      write (MM, "(I2.2)") Rmonth
      write (DD, "(I2.2)") Rday
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(out_data(N_xi_rho, N_eta_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
!---- Read DSJRA-55 GRIB2 file --------------------------------

      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth
      GRIB_suffix(2:5)=YYYY
      GRIB_suffix(6:7)=MM

      GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
      !Open GRIB file
      call codes_grib_multi_support_on	(	iret	)	
      write(*,*) "OPEN: ", trim( GRIB_FILE )
      call codes_open_file(ifile, GRIB_FILE,'r', iret)
      call codes_grib_new_from_file(ifile,igrib, iret)
!
!      ! Get dimension data
      call codes_get(igrib,'Ny', Jm)
      call codes_get(igrib,'Nx', Im)
      write(*,*) Im,Jm
!          
      ! Allocate variable
      allocate(values(Im*Jm))
      allocate(lat(Im,Jm))
      allocate(lon(Im,Jm))
      allocate(grib_data(N_Param,Im, Jm))
      
!      open(unit=21, file='lat.txt')
!      open(unit=22, file='lon.txt')
      
!  ---- Get Lat Lon coordinates from GRIB2 file --------------
!
!      call codes_get(igrib,'latitudes', values)
!      write(*,*) values(1)
!      do i=1, Jm
!        istart = 1 + Im*(i-1)
!        iend   = Im*i
!        lat(:,i) = values(istart:iend)
!!        write(21,*) lat(:,i)
!      end do
!      call codes_get(igrib,'longitudes', values)
!      write(*,*) values(1)
!      do i=1, Jm
!        istart = 1 + Im*(i-1)
!        iend   = Im*i
!        lon(:,i) = values(istart:iend)
!!        write(22,*) lon(:,i)
!      end do
!
!      write(*,*) lat(1,1), lon(1,1)
!      write(*,*) lat(1,Jm), lon(1,Jm)
!      write(*,*) lat(Im,1), lon(Im,1)
!      write(*,*) lat(Im,Jm), lon(Im,Jm)
      
      call codes_release(igrib)
      call codes_close_file(ifile)

!  ---- Get Lat Lon coordinates from Binary file --------------
!
      open(unit=20, file=GRIB_LL, action='read',                 &
           & form='unformatted', access='direct', recl=4,        &
           & CONVERT='BIG_ENDIAN', status='old')
      do j=1, Jm
        do i=1,Im
          read(20, rec=i+Im*(j-1)) bd
          lat(i,j) =dble(bd)
!          write(*,*) i+Im*(j-1), lat(i,j)
!           write(*,*) lat(i,j)
        end do
!        write(21,*) lat(:,j)
      end do
      do j=1, Jm
        do i=1,Im
          read(20, rec=Im*Jm+i+Im*(j-1)) bd
          lon(i,j) = dble(bd)
!          write(*,*) i+Im*(j-1), lat(i,j)
        end do
!        write(22,*) lon(:,j)
      end do
      close(20)
!      close(21)
!      close(22)
      write(*,*) 'NW corner:', lat(1,1),  lon(1,1)
      write(*,*) 'SW corner:', lat(1,Jm), lon(1,Jm)
      write(*,*) 'SE corner:', lat(Im,1), lon(Im,1)
      write(*,*) 'NE corner:', lat(Im,Jm),lon(Im,Jm)
      
      
!---- Create the forcing netCDF file --------------------------------

      OUT_suffix(2:5)=YYYY
      OUT_suffix(6:7)=MM

      DO iparam=1,N_Param
      
        OUT_FILE(iparam) = OUT_prefix//trim(NC_NAME(iparam))//OUT_suffix
        
        write(*,*) "CREATE: ", trim( OUT_FILE(iparam) )

        call check( nf90_create(OUT_FILE(iparam), nf90_clobber, ncid) )

        call check( nf90_def_dim(ncid, 'xi_rho', N_xi_rho, xi_rho_dimid) )
        call check( nf90_def_dim(ncid, 'eta_rho',N_eta_rho, eta_rho_dimid) )
        call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid) )

        dimids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)

    !   Define the netCDF variables for the pressure and temperature data.
        call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, var_id) )
        call check( nf90_put_att(ncid, var_id, 'long_name', 'atmospheric forcing time') )
        call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

        call check( nf90_def_var(ncid, trim( NC_NAME(iparam) ), NF90_REAL, dimids, var_id) )
        call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(iparam) )) )
        call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(iparam) ) ) )
        call check( nf90_put_att(ncid, var_id, 'time',      'time') )

  ! End define mode.
        call check( nf90_enddef(ncid) )
        call check( nf90_close(ncid) )

      END DO
      
!---- LOOP1 START --------------------------------

      itime = 1
      
      DO
      
        write (DD, "(I2.2)") 1+int((itime-1)*1/24)
        write (hh, "(I2.2)") mod((itime-1)*1,24)
        GRIB_suffix(8:9)=DD
        GRIB_suffix(10:11)=hh
      
        GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
        
        !Open GRIB file
        write(*,*) "OPEN: ", trim( GRIB_FILE )
        call codes_open_file(ifile, GRIB_FILE,'r', iret)
        if (iret /= CODES_SUCCESS) then
          write(*,*) "CANNOT OPEN: ", trim( GRIB_FILE )
          exit
        end if
        call codes_grib_new_from_file(ifile,igrib, iret)
        write(*,*) "READ: ", trim( GRIB_FILE )
        call codes_get(igrib,'validityDate',YYYYMMDD)
        write(*,*) 'validityDate=', YYYYMMDD
        call codes_get(igrib,'validityTime',hhmm)
        write(*,*) 'validityTime=', hhmm
        call codes_release(igrib)
        call codes_close_file(ifile)
        
        iyear  = YYYYMMDD/10000
        imonth = (YYYYMMDD-iyear*10000)/100
        iday   = YYYYMMDD-iyear*10000-imonth*100
        ihour  = hhmm/100
        imin   = hhmm-100*ihour
        call ndays(imonth, iday, iyear, Rmonth, Rday, Ryear, idays)
        
        time(1) = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
        

!  ----   LOOP2 START --------------------------------
        DO iparam=1,N_Param
        
          if(iparam==7) then !!! for rain (rain fall rate)
            GRIB_FILE = trim(GRIB_FCST_PHY2M_dir)//GRIB_FCST_PHY2M_prefix//GRIB_suffix
          else
            GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
          end if
          
          write(*,*) "READ: ", trim( GRIB_FILE )
          call codes_index_create(idx,GRIB_FILE,'productDefinitionTemplateNumber:i,parameterCategory:i,parameterNumber:i')
      
          call codes_index_select(idx,'productDefinitionTemplateNumber',GRIB_NAME(1,iparam),iret)
          call codes_index_select(idx,'parameterCategory',GRIB_NAME(2,iparam),iret)
          call codes_index_select(idx,'parameterNumber',GRIB_NAME(3,iparam),iret)
          call codes_new_from_index(idx,igrib, iret)
          call codes_get(igrib,'values', values)
          do i=1, Jm
            istart = 1 + Im*(i-1)
            iend   = Im*i
            grib_data(:,i) = values(istart:iend)
          end do
          
          if(iparam==3) then  !!! for Tair
            grib_data = grib_data - 273.15d0  ! K -> degC
          end if
          if(iparam==5) then !!! for Pair (Pressure)
            grib_data = grib_data*0.01  ! Pa -> millibar (= hPa)
          end if
          if(iparam==6) then !!! for cloud (cloud fraction)
            grib_data = grib_data*0.01  ! percent -> ratio(0 to 1)
          end if
          if(iparam==7) then !!! for rain (rain fall rate)
            grib_data = grib_data   ! kg m-2 s-1
            time(1) = time(1)+0.5d0/24.0d0  !!! + 0.5 hours
          end if
          
          write(*,*) time(1),TIME_ATT
          
          write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )
!          call LinearInterpolation2D_grid2(Im, Jm, lon, lat               &
!     &                    , grib_data, -9999.0d0, 9999.0d0                &
!     &                    , N_xi_rho, N_eta_rho, lon_rho, lat_rho, out_data )
          
          start1D = (/ itime /)
          count1D = (/ 1 /)
          call writeNetCDF_1d(           &
!              input parameters
     &            'time'                 &
     &          , OUT_FILE(iparam)       &
     &          , 1                      &
     &          , time                   &
     &          , start1D, count1D       &
     &          )

          start3D = (/ 1,  1,  itime /)
          count3D = (/ N_xi_rho, N_eta_rho, 1 /)
          
          call writeNetCDF_3d(                    &
!            input parameters
     &          trim( NC_NAME(iparam) )           &
     &        , OUT_FILE(iparam)                  &
     &        , N_xi_rho, N_eta_rho, 1            &
     &        , out_data                          &
     &        , start3D, count3D                  &
     &        )
          
          
          call codes_index_release(idx)
          
        END DO
!  ---- LOOP2 END --------------------------------
        call codes_release(igrib)

        itime = itime + 1

      END DO
!---- LOOP1 END --------------------------------

      deallocate(values)
!      deallocate(lat)
!      deallocate(lon)
      deallocate(grib_data)
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
    END PROGRAM frcDSJRA55
      
