
!!!=== Copyright (c) 2017 Takashi NAKAMURA  =====

    PROGRAM frcJRA55
      use netcdf
      use grib_api
      use mod_roms_netcdf
      use mod_interpolation
      use mod_calendar
     
      implicit none
      
! SETTINGS of input and output files  -----------------------------------

      integer, parameter :: Ryear  = 2000   ! Reference year
      integer, parameter :: Rmonth = 1      ! Reference month
      integer, parameter :: Rday   = 1      ! Reference day

     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "D:/ROMS/Data/Coral_Triangle/CT_0.04_grd_v1.nc"

      integer, parameter :: N_Param = 7
      character(len=*), parameter :: GRIB_dir  = "D:/JRA-55/"
      character(256) :: GRIB_prefix(N_Param) = (/ &
     &   "fcst_surf125.011_tmp.   "               &
     &  ,"fcst_surf125.033_ugrd.  "               &
     &  ,"fcst_surf125.034_vgrd.  "               &
     &  ,"fcst_surf125.052_rh.    "               &
     &  ,"fcst_surf125.002_prmsl. "               &
     &  ,"fcst_surf125.071_tcdc.  "               &
     &  ,"fcst_phy2m125.061_tprat."               &
     &  /)
!      character(len=*), parameter :: GRIB_suffix  = "2016030100_2016033121"
      character(len=*), parameter :: GRIB_suffix  = "2016040100_2016043021"
      
      character(len=*), parameter :: OUT_prefix  = "output/CT_0.04_frc_JRA55_"
      character(len=*), parameter :: OUT_suffix  = "_1604.nc"
      
      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
      
!----------------------------------------------------------------------

      character(256) :: GRIB_FILE(N_Param)
      character(256) :: OUT_FILE(N_Param)
      
      character(256) :: NC_NAME(N_Param) = (/     &
     &   "Tair "                                  &
     &  ,"Uwind"                                  &
     &  ,"Vwind"                                  &
     &  ,"Qair "                                  &
     &  ,"Pair "                                  &
     &  ,"cloud"                                  &
     &  ,"rain "                                  &
     &  /)
      character(256) :: NC_LNAME(N_Param) = (/    &
     &   "surface air temperature      "          &
     &  ,"surface u-wind component     "          &
     &  ,"surface v-wind component     "          &
     &  ,"surface air relative humidity"          &
     &  ,"surface air pressure         "          &
     &  ,"cloud fraction               "          &
     &  ,"rain fall rate               "          &
     &  /)
      character(256) :: NC_UNIT(N_Param) = (/     &
     &   "Celsius                  "              &
     &  ,"meter second-1           "              &
     &  ,"meter second-1           "              &
     &  ,"percentage               "              &
     &  ,"millibar                 "              &
     &  ,"0 to 1                   "              &
     &  ,"kilogram meter-2 second-1"              &
     &  /)
     
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8), allocatable :: out_data(:,:) ! output forcing data
           
      real(8), allocatable :: lat(:), lon(:)
      real(8), allocatable :: grib_data(:, :)
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
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      integer :: Im, Jm
      
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
      
      integer :: iparam
      integer :: ifile,iret,igrib
      integer :: istart, iend
      integer :: YYYYMMDD, hhmm
      real(8), allocatable :: values(:)

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
      
!---- LOOP1 START --------------------------------

      DO iparam=7,N_Param
      
!---- Create the forcing netCDF file --------------------------------
      
        OUT_FILE(iparam) = OUT_prefix//trim(NC_NAME(iparam))//OUT_suffix
        
        write(*,*) "CREATE: ", OUT_FILE(iparam)

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

      
!---- Read JRA-55 GRIB1 file --------------------------------

        GRIB_FILE(iparam) = GRIB_dir//trim(GRIB_prefix(iparam))//GRIB_suffix
        !Open GRIB file
        write(*,*) "OPEN: ", trim( GRIB_FILE(iparam) )
        call grib_open_file(ifile, GRIB_FILE(iparam),'r')
        call grib_new_from_file(ifile,igrib, iret)

        ! Get dimension data
        call grib_get(igrib,'Nj', Jm)
        call grib_get(igrib,'Ni', Im)
            
        write(*,*) Im, Jm
        
!        ! Allocate variable
        allocate(lat(Jm))
        allocate(lon(Im))
        
        call grib_get(igrib,'distinctLatitudes',lat)
        call grib_get(igrib,'distinctLongitudes',lon)
        
        allocate(values(Im*Jm))
        allocate(grib_data(Im, Jm))

!----   LOOP2 START --------------------------------
        itime = 1
        
        DO WHILE (iret /= GRIB_END_OF_FILE)

          write(*,*) "READ: ", trim( GRIB_FILE(iparam) )
          call grib_get(igrib,'validityDate',YYYYMMDD)
          write(*,*) 'validityDate=', YYYYMMDD
          call grib_get(igrib,'validityTime',hhmm)
          write(*,*) 'validityTime=', hhmm
          
          iyear  = YYYYMMDD/10000
          imonth = (YYYYMMDD-iyear*10000)/100
          iday   = YYYYMMDD-iyear*10000-imonth*100
          ihour  = hhmm/100
          imin   = hhmm-100*ihour
          call ndays(imonth, iday, iyear, Rmonth, Rday, Ryear, idays)
          
          time(1) = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
          
          call grib_get(igrib,'values', values)
          do i=1, Jm
            istart = 1 + Im*(i-1)
            iend   = Im*i
            grib_data(:,i) = values(istart:iend)
          end do
          
          if(iparam==1) then  !!! for Tair
            grib_data = grib_data - 273.15d0  ! K -> degC
          end if
          if(iparam==5) then !!! for Pair (Pressure)
            grib_data = grib_data*0.01  ! Pa -> millibar (= hPa)
          end if
          if(iparam==6) then !!! for cloud (cloud fraction)
            grib_data = grib_data*0.01  ! percent -> ratio(0 to 1)
          end if
          if(iparam==7) then !!! for rain (rain fall rate)
            grib_data = grib_data*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1
                                 !!! 1mm x 1m x 1m = 0.1x100x100 cm3 = 1000 cm3 = 1L = 1kg
            time(1) = time(1)-1.5d0/24.0d0  !!! minus 1.5 hours
          end if
          
          write(*,*) time(1),TIME_ATT
          
          write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )
          call LinearInterpolation2D_grid2(Im, Jm, lon, lat               &
     &                    , grib_data, -9999.0d0, 9999.0d0                &
     &                    , N_xi_rho, N_eta_rho, lon_rho, lat_rho, out_data )
          
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
          
          
          call grib_new_from_file(ifile,igrib, iret)
          itime = itime + 1

        END DO
!---- LOOP2 END --------------------------------
      
        deallocate(values)
        deallocate(lat)
        deallocate(lon)
        deallocate(grib_data)
      
      END DO
!---- LOOP1 END --------------------------------
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
    END PROGRAM frcJRA55
      
