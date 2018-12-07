
!!!=== Copyright (c) 2014-2018 Takashi NAKAMURA  =====

    PROGRAM frcJMAMSM
      use netcdf
      use mod_roms_netcdf
      use mod_interpolation
      use mod_calendar
      use mod_utility
           
      implicit none

!-------------------------------------------------------------------------------
      integer :: Syear, Smonth, Sday
      character(256) :: GRID_FILE
      character(256) :: MSM_SURF_dir
      character(256) :: MSM_R1H_dir
      character(256) :: FRC_prefix
      integer :: Ryear, Rmonth, Rday
      integer :: mode
!----------------------------------------------------------------------

      integer, parameter :: N_Param = 7
      
      character(12) :: MSM_suffix  = "2002/0101.nc"
      character(10) :: FRC_suffix   = "_195801.nc"
      
      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"

      character(256) :: MSM_FILE
      character(256) :: FRC_FILE(N_Param)
      
      character(256) :: MSM_NAME(N_Param) = (/    &
     &   "u   "                                   &
     &  ,"v   "                                   &
     &  ,"temp"                                   &
     &  ,"rh  "                                   &
     &  ,"psea"                                   &
     &  ,"ncld"                                   &
     &  ,"r1h "                                   &
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
      !
      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8), allocatable :: out_data(:,:,:) ! output forcing data
           
      real(8), allocatable :: lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:)
      real(8), allocatable :: in_data(:,:,:)
      real(8), allocatable :: time(:)
      integer, allocatable :: ID_cont(:,:)
      real(8), allocatable :: w_cont(:,:)
      integer :: start1Dr(1), count1Dr(1)
      integer :: start1Dw(1), count1Dw(1)
      integer :: start3Dr(3), count3Dr(3)
      integer :: start3Dw(3), count3Dw(3)
      
      integer :: iyear, imonth, iday
      integer :: ihour, imin
      integer :: i,j,k
      integer :: idays
      integer :: itime, days_since20000101
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      character(2) :: hh
      character(11) :: YYYYMMDDpHH
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      integer :: Im, Jm, Nt
      integer :: I1, I2, J1, J2
      real(8) :: Tlat, Blat, Llon, Rlon
      real(8) :: sf, off
      
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
      
      integer :: iparam
      integer :: ifile,idx,iret,igrib
      integer :: istart, iend
      integer :: YYYYMMDD, hhmm

!-------------------------------------------------------------------------------
      namelist/grd/GRID_FILE
      namelist/sdate/Syear, Smonth, Sday
      namelist/refdate/Ryear, Rmonth, Rday
      namelist/intpmode/mode
      namelist/frc_jmamsm/MSM_SURF_dir
      namelist/frc_jmamsm/MSM_R1H_dir
      namelist/frc_jmamsm/FRC_prefix

      ! Read parameters in namelist file
      
      read (*, nml=grd)
      read (*, nml=sdate)
      read (*, nml=refdate)
      read (*, nml=intpmode)
      read (*, nml=frc_jmamsm)

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
      MSM_suffix(1:4)=YYYY
      MSM_suffix(6:7)=MM

      MSM_FILE = trim(MSM_SURF_dir)//MSM_suffix
      !Open GRIB file
      write(*,*) "OPEN: ", trim( MSM_FILE )
      !Open NetCDF file
      call check( nf90_open(trim( MSM_FILE ), nf90_nowrite, ncid) )

      ! Get dimension data
      call get_dimension(ncid, 'lat', Jm)
      call get_dimension(ncid, 'lon', Im)
          
      write(*,*) Im, Jm, Nt
      
!      ! Allocate variable
      allocate(lat_all(Jm))
      allocate(lon_all(Im))
      
      call check( nf90_inq_varid(ncid, 'lat', var_id) )
      call check( nf90_get_var(ncid, var_id, lat_all) )
      call check( nf90_inq_varid(ncid, 'lon', var_id) )
      call check( nf90_get_var(ncid, var_id, lon_all) )
      call check( nf90_close(ncid) )
      
      call min_max_2D(N_xi_rho, N_eta_rho, lat_rho, Blat, Tlat)
      call min_max_2D(N_xi_rho, N_eta_rho, lon_rho, Llon, Rlon)
      write(*,*) Blat, Tlat
      write(*,*) Llon, Rlon
      
      call seek_id_range(Jm, lat_all, Blat, Tlat, J2, J1)
      call seek_id_range(Im, lon_all, Llon, Rlon, I1, I2)

      write(*,*) J1, lat_all(J1), J2, lat_all(J2)
      write(*,*) I1, lon_all(I1), I2, lon_all(I2)
      
      Jm=J2-J1+1
      Im=I2-I1+1
      
      
      
      allocate(lat(Jm))
      allocate(lon(Im))
      allocate(ID_cont(4, N_xi_rho*N_eta_rho))
      allocate(w_cont(4, N_xi_rho*N_eta_rho))
      
      lat(:) = lat_all(J1:J2)
      lon(:) = lon_all(I1:I2)
      write(*,*) Im,Jm, lon(1), lon(Im), lat(1), lat(Jm)

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
      
      
!==== LOOP1 START ==============================================================
      itime = 1
      
      DO
      
        write (DD, "(I2.2)") itime
        MSM_suffix(8:9)=DD

!----   LOOP2 START --------------------------------
        DO iparam=1,N_Param
        
          if(iparam==7) then !!! for rain (rain fall rate)
            MSM_FILE = trim(MSM_R1H_dir)//MSM_suffix
          else
            MSM_FILE = trim(MSM_SURF_dir)//MSM_suffix
          end if
          
          write(*,*) "READ: ", trim( MSM_FILE )
          !Open NetCDF file
          call check( nf90_open(trim( MSM_FILE ), nf90_nowrite, ncid) )
          call get_dimension(ncid, 'time', Nt)
          allocate(time(Nt))
          allocate(in_data(Im, Jm, Nt))
          allocate(out_data(N_xi_rho, N_eta_rho, Nt))

          ! Define the netCDF variables for the pressure and temperature data.
          call check( nf90_inq_varid(ncid, 'time', var_id) )  !!!  not Japan time (00:00:00+09:00)
          call check( nf90_get_var(ncid, var_id, time) )
          call ndays(Smonth, itime, Syear, 1, 1, 2000, days_since20000101)
          time = time/24.0d0 + dble(days_since20000101)

          ! Initialize count and start for reading data
          start3Dr = (/ I1, J1, 1 /)
          count3Dr = (/ Im, Jm, Nt /)

          ! Get variable id
          call check( nf90_inq_varid(ncid, trim( MSM_NAME(iparam) ), var_id) ) ! Rainfall in 1 hour (mm/h)
          call check( nf90_get_var(ncid, var_id, in_data, start=start3Dr, count=count3Dr) )
          call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
          call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
          in_data=in_data*sf+off
          ! Close NetCDF file
          call check( nf90_close(ncid) )
          if(iparam==3) then  !!! for Tair
            in_data = in_data - 273.15d0  ! K -> degC
          end if
          if(iparam==5) then !!! for Pair (Pressure)
            in_data = in_data*0.01  ! Pa -> millibar (= hPa)
          end if
          if(iparam==6) then !!! for cloud (cloud fraction)
            in_data = in_data*0.01  ! percent -> ratio(0 to 1)
          end if
          if(iparam==7) then !!! for rain (rain fall rate)
            in_data = in_data*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1
                                 !!! 1mm x 1m x 1m = 0.1x100x100 cm3 = 1000 cm3 = 1L = 1kg
            time = time-0.5d0/24.0d0  !!! -0.5 hours
          end if
          
          
          do k=1, Nt
          
            CALL cdate2(time(k),Ryear,Rmonth,Rday,YYYYMMDDpHH)
            write(*,*) 'time = ', YYYYMMDDpHH
          
            write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )

            call interp2D_grid(Im, Jm, in_data(:,:,k)                         &
     &                       , N_xi_rho, N_eta_rho, out_data(:,:,k)           &
     &                       , Id_cont, w_cont)
     
          end do
          
          start1Dw = (/ 1+Nt*(itime-1) /)
          count1Dw = (/ Nt /)
          call writeNetCDF_1d(                &
!              input parameters
     &            'time'                      &
     &          , trim( FRC_FILE(iparam) )    &
     &          , Nt                          &
     &          , time                        &
     &          , start1Dw, count1Dw          &
     &          )

          start3Dw = (/ 1,  1,  1+Nt*(itime-1) /)
          count3Dw = (/ N_xi_rho, N_eta_rho, Nt /)
          
          call writeNetCDF_3d(                    &
!            input parameters
     &          trim( NC_NAME(iparam) )           &
     &        , trim( FRC_FILE(iparam) )          &
     &        , N_xi_rho, N_eta_rho, Nt           &
     &        , out_data                          &
     &        , start3Dw, count3Dw                &
     &        )
          
          deallocate(time)
          deallocate(in_data)
          deallocate(out_data)
          
        END DO
!---- LOOP2 END --------------------------------

        itime = itime + 1

      END DO
!==== LOOP1 END ================================================================

      deallocate(lat)
      deallocate(lon)
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
    END PROGRAM frcJMAMSM
      
