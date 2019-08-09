
!!!=== Copyright (c) 2018 Takashi NAKAMURA  =====

    PROGRAM frcJRA55
      use netcdf
      use eccodes
      use mod_roms_netcdf
      use mod_interpolation
      use mod_calendar
     
      implicit none
      
!-------------------------------------------------------------------------------
      integer :: Syear, Smonth, Sday
      character(256) :: GRID_FILE
      character(256) :: GRIB_FCST_SURF_dir
      character(256) :: GRIB_FCST_PHY2M_dir
      character(256) :: FRC_prefix
      integer :: Ryear, Rmonth, Rday
      integer :: mode
!----------------------------------------------------------------------
      integer, parameter :: N_Param = 7
      character(len=*), parameter :: GRIB_FCST_SURF_prefix  = "fcst_surf125"
      character(len=*), parameter :: GRIB_FCST_PHY2M_prefix = "fcst_phy2m125"
      character(11) :: GRIB_suffix  = ".1958010100"
      character(10) :: FRC_suffix   = "_195801.nc"
      
      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"

      character(256) :: GRIB_FILE
      character(256) :: FRC_FILE(N_Param)
      
      character(256) :: GRIB_NAME(N_Param) = (/   &
     &   "10u     "                               &
     &  ,"10v     "                               &
     &  ,"2t      "                               &
     &  ,"2r      "                               &
     &  ,"msl     "                               &
     &  ,"tcc     "                               &
     &  ,"tpratsfc"                               &
     &  /)
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

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8), allocatable :: out_data(:,:) ! output forcing data
           
      real(8), allocatable :: lat(:), lon(:)
      real(8), allocatable :: grib_data(:, :)
      real(8) :: time(1)
      integer, allocatable :: ID_cont(:,:)
      real(8), allocatable :: w_cont(:,:)
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
      character(11) :: YYYYMMDDpHH
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      integer :: Im, Jm
      
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
      
      integer :: iparam
      integer :: ifile,idx,iret,igrib
      integer :: istart, iend
      integer :: YYYYMMDD, hhmm
      real(8), allocatable :: values(:)

!-------------------------------------------------------------------------------
      namelist/grd/GRID_FILE
      namelist/sdate/Syear, Smonth, Sday
      namelist/refdate/Ryear, Rmonth, Rday
      namelist/intpmode/mode
      namelist/frc_jra55/GRIB_FCST_SURF_dir
      namelist/frc_jra55/GRIB_FCST_PHY2M_dir
      namelist/frc_jra55/FRC_prefix

      ! Read parameters in namelist file
      
      read (*, nml=grd)
      read (*, nml=sdate)
      read (*, nml=refdate)
      read (*, nml=intpmode)
      read (*, nml=frc_jra55)

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
      allocate(ID_cont(4, N_xi_rho*N_eta_rho))
      allocate(w_cont(4, N_xi_rho*N_eta_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
!---- Read JRA-55 GRIB1 file --------------------------------

      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth
      GRIB_suffix(2:5)=YYYY
      GRIB_suffix(6:7)=MM

      GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
      !Open GRIB file
      write(*,*) "OPEN: ", trim( GRIB_FILE )
      call codes_open_file(ifile, GRIB_FILE,'r')
      call codes_grib_new_from_file(ifile,igrib, iret)

      ! Get dimension data
      call codes_get(igrib,'Nj', Jm)
      call codes_get(igrib,'Ni', Im)
          
      write(*,*) Im, Jm
      
!      ! Allocate variable
!      allocate(lat(Jm))
!      allocate(lon(Im))
      
      call codes_get(igrib,'distinctLatitudes',lat)
      call codes_get(igrib,'distinctLongitudes',lon)
      
      call codes_release(igrib)
      call codes_close_file(ifile)
      
      allocate(values(Im*Jm))
      allocate(grib_data(Im, Jm))
      
      write(*,*) "CALC.: weight parameters for interpolation"
      call weight2D_grid(Im,Jm,lon,lat,N_xi_rho,N_eta_rho,lon_rho,lat_rho,ID_cont,w_cont)
      
!---- Create the forcing netCDF file --------------------------------

      FRC_suffix(2:5)=YYYY
      FRC_suffix(6:7)=MM

      DO iparam=1,N_Param
      
        FRC_FILE(iparam) = trim(FRC_prefix)//'_'//trim(NC_NAME(iparam))//FRC_suffix
        
        write(*,*) "CREATE: ", trim( FRC_FILE(iparam) )

        call check( nf90_create(trim( FRC_FILE(iparam) ), nf90_clobber, ncid) )

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
      
      DO WHILE (iret /= GRIB_END_OF_FILE)
      
        write (DD, "(I2.2)") 1+int((itime-1)*3/24)
        write (hh, "(I2.2)") mod((itime-1)*3,24)
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
        call codes_close_file(ifile) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        iyear  = YYYYMMDD/10000
        imonth = (YYYYMMDD-iyear*10000)/100
        iday   = YYYYMMDD-iyear*10000-imonth*100
        ihour  = hhmm/100
        imin   = hhmm-100*ihour
        call ndays(imonth, iday, iyear, Rmonth, Rday, Ryear, idays)
        
        time(1) = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
        

!----   LOOP2 START --------------------------------
!        DO iparam=1,N_Param
        DO iparam=1,N_Param
        
          if(iparam==7) then !!! for rain (rain fall rate)
            GRIB_FILE = trim(GRIB_FCST_PHY2M_dir)//GRIB_FCST_PHY2M_prefix//GRIB_suffix
          else
            GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
          end if
          
          write(*,*) "READ: ", trim( GRIB_FILE )
          call codes_index_create(idx,trim( GRIB_FILE ),'shortName')
      
          call codes_index_select(idx,'shortName',GRIB_NAME(iparam))
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
            grib_data = grib_data*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1
                                 !!! 1mm x 1m x 1m = 0.1x100x100 cm3 = 1000 cm3 = 1L = 1kg
            time(1) = time(1)+1.5d0/24.0d0  !!! + 1.5 hours
          end if
          
          CALL oceantime2cdate(time(1),Ryear,Rmonth,Rday,YYYYMMDDpHH)
          write(*,*) 'time = ', YYYYMMDDpHH
          
          write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )
          call interp2D_grid(Im, Jm, grib_data                             &
     &                       , N_xi_rho, N_eta_rho, out_data               &
     &                       , Id_cont, w_cont)
          
          start1D = (/ itime /)
          count1D = (/ 1 /)
          call writeNetCDF_1d(                &
!              input parameters
     &            'time'                      &
     &          , trim( FRC_FILE(iparam) )    &
     &          , 1                           &
     &          , time                        &
     &          , start1D, count1D            &
     &          )

          start3D = (/ 1,  1,  itime /)
          count3D = (/ N_xi_rho, N_eta_rho, 1 /)
          
          call writeNetCDF_3d(                    &
!            input parameters
     &          trim( NC_NAME(iparam) )           &
     &        , trim( FRC_FILE(iparam) )          &
     &        , N_xi_rho, N_eta_rho, 1            &
     &        , out_data                          &
     &        , start3D, count3D                  &
     &        )
          
          
          call codes_index_release(idx)
          
        END DO
!---- LOOP2 END --------------------------------
        call codes_release(igrib)

        itime = itime + 1

      END DO
!---- LOOP1 END --------------------------------

      deallocate(values)
      deallocate(lat)
      deallocate(lon)
      deallocate(grib_data)
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
    END PROGRAM frcJRA55
      
