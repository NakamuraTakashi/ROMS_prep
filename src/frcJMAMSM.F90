
!!!=== Copyright (c) 2014-2020 Takashi NAKAMURA  =====

#ifdef NETCDF_INPUT
# undef SWRAD
#endif

PROGRAM frcJMAMSM
  use netcdf
  use eccodes
  use mod_roms_netcdf
  use mod_interpolation
  use mod_calendar
 
  implicit none
  
!-------------------------------------------------------------------------------
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  character(256) :: GRID_FILE
  character(256) :: MSM_dir
  character(256) :: FRC_prefix
  integer :: Ryear, Rmonth, Rday
!----------------------------------------------------------------------
#if defined SWRAD
  integer, parameter :: N_Param = 11
#else
  integer, parameter :: N_Param = 10
#endif

#if defined NETCDF_INPUT

  character(256) :: NCIN_NAME(N_Param) = (/   &
     "psea      "                             &
    ,"u         "                             &
    ,"v         "                             &
    ,"temp      "                             &
    ,"rh        "                             &
    ,"ncld_low  "                             &
    ,"ncld_mid  "                             &
    ,"ncld_upper"                             &
    ,"ncld      "                             &
    ,"r1h       "                             &
    /)

    integer :: d_ref_days
    real(8) :: sf, off
#else

  character(len=*), parameter :: GRIB_prefix  = "Z__C_RJTD_"
  character(len=*), parameter :: GRIB_suffix  =   &
        "0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"
  character(10) :: GRIB_yyyymmddhh = "2002070100"

  integer :: ifile,idx,iret,igrib
  integer :: istart, iend
  integer :: YYYYMMDD(N_Param), hhmm(N_Param)
  real(8), allocatable :: values(:)
  integer :: p1
  character(256) :: p2

  character(256) :: GRIB_NAME(N_Param) = (/   &
     "Pressure reduced to MSL           "     &
    ,"u-component of wind               "     &
    ,"v-component of wind               "     &
    ,"Temperature                       "     &
    ,"Relative humidity                 "     &
    ,"Low cloud cover                   "     &
    ,"Medium cloud cover                "     &
    ,"High cloud cover                  "     &
    ,"Total cloud cover                 "     &
    ,"Total precipitation               "     &
# if defined SWRAD
    ,"Downward short-wave radiation flux"     &
# endif
    /)
  integer :: GRIB_STEP(15) = (/ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 /)

#endif

  character(9) :: FRC_yyyymmdd = "_20020701" 
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: FRC_FILE(2)

  character(256) :: NC_NAME(N_Param) = (/     &
     "Pair  "                                 &
    ,"Uwind "                                 &
    ,"Vwind "                                 &
    ,"Tair  "                                 &
    ,"Qair  "                                 &
    ,"lcloud"                                 &
    ,"mcloud"                                 &
    ,"hcloud"                                 &
    ,"cloud "                                 &
    ,"rain  "                                 &
#if defined SWRAD
    ,"swrad "                                 &
#endif
    /)
  character(256) :: NC_LNAME(N_Param) = (/    &
     "surface air pressure          "         &
    ,"surface u-wind component      "         &
    ,"surface v-wind component      "         &
    ,"surface air temperature       "         &
    ,"surface air relative humidity "         &
    ,"low cloud fraction            "         &
    ,"medium cloud fraction         "         &
    ,"high cloud fraction           "         &
    ,"cloud fraction                "         &
    ,"rain fall rate                "         &
#if defined SWRAD
    ,"solar shortwave radiation flux"         &
#endif
    /)
  character(256) :: NC_UNIT(N_Param) = (/     &
     "millibar                 "              &
    ,"meter second-1           "              &
    ,"meter second-1           "              &
    ,"Celsius                  "              &
    ,"percentage               "              &
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"kilogram meter-2 second-1"              &
#if defined SWRAD
    ,"watt meter-2             "              &
#endif
    /)

  real(8), parameter :: PI = 3.141592653589793d0
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: cosAu(:, :)
  real(8), allocatable :: sinAu(:, :)
  real(8), allocatable :: cosAv(:, :)
  real(8), allocatable :: sinAv(:, :)
  
  real(8), allocatable :: in_data(:,:, :)
  real(8), allocatable :: out_data(:,:,:) ! output forcing data
       
  integer :: Im, Jm
  real(8), allocatable :: lat(:), lon(:)
  real(8) :: t
  real(8) :: time(1)
  integer, allocatable :: ID_cont(:,:)
  real(8), allocatable :: w_cont(:,:)
  integer :: start1D(1), count1D(1)
  integer :: start3D(3), count3D(3)
  
  character(256) :: IN_FILE(2), IN_FILE2(2)
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: i,j,k
  integer :: idays,ihours,juliandate
  integer :: itime
  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh
  
  integer :: ncid,var_id
  integer :: Nxr, Nyr !, Nxu, Nyu, Nxv, Nyv
  integer :: L, M  
  integer :: dimids(3)
  integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
  
  integer :: iparam,ifc,inc,isp
  real(8) :: d_lat, d_lon
  real(8) :: u, v
!  real(8) :: dpT, sat_VP, VP
  logical :: file_exists
!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/frc_jmamsm/MSM_dir
  namelist/frc_jmamsm/FRC_prefix
  ! Read parameters in namelist file
  
  read (*, nml=grd)
  read (*, nml=sdate)
  read (*, nml=edate)
  read (*, nml=refdate)
  read (*, nml=frc_jmamsm)

!---- Modify time-unit description ---------------------------------
      
  write (YYYY, "(I4.4)") Ryear
  write (MM, "(I2.2)") Rmonth
  write (DD, "(I2.2)") Rday
  
  TIME_ATT(12:15)=YYYY
  TIME_ATT(17:18)=MM
  TIME_ATT(20:21)=DD
      
!---- Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
!  call get_dimension(ncid, 'xi_u',    Nxu)
!  call get_dimension(ncid, 'eta_u',   Nyu)
!  call get_dimension(ncid, 'xi_v',    Nxv)
!  call get_dimension(ncid, 'eta_v',   Nyv)
  L = Nxr-1
  M = Nyr-1
  allocate( latr(0:L, 0:M) )
  allocate( lonr(0:L, 0:M) )
  allocate( cosAu(0:L, 0:M) )
  allocate( sinAu(0:L, 0:M) )
  allocate( cosAv(0:L, 0:M) )
  allocate( sinAv(0:L, 0:M) )
  allocate(out_data(N_Param, 0:L, 0:M))
  allocate(ID_cont(4, Nxr*Nyr))
  allocate(w_cont(4, Nxr*Nyr))

  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr) )
  
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
  do i=0,L-1
    do j=0,M
      d_lat=latr(i+1,j)-latr(i,j)
      d_lon=lonr(i+1,j)-lonr(i,j)
      d_lon=d_lon*cos(latr(i,j)/180.0d0*PI)
      cosAu(i,j) = d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAu(i,j) = d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo
  enddo
  cosAu(L,:) = cosAu(L-1,:)
  sinAu(L,:) = sinAu(L-1,:)

  do i=0,L
    do j=0,M-1
      d_lat=latr(i,j+1)-latr(i,j)
      d_lon=lonr(i,j)-lonr(i,j+1)
      d_lon=d_lon*cos(latr(i,j)/180.0d0*PI)
      cosAv(i,j) = d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAv(i,j) = d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo         
  enddo
  cosAv(:,M) = cosAv(:,M-1)
  sinAv(:,M) = sinAv(:,M-1)

!---- Read JMA-MSM file --------------------------------
  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  write (DD, "(I2.2)") Sday

#if defined NETCDF_INPUT
!---- Read JMA-MSM NetCDF file --------------------------------
  IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//DD//".nc"
  IN_FILE(2) = trim(MSM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

  write(*,*) "OPEN: ", trim( IN_FILE(1) )
  !Open NetCDF file
  call check( nf90_open(trim( IN_FILE(1) ), nf90_nowrite, ncid) )

  ! Get dimension data
  call get_dimension(ncid, 'lat', Jm)
  call get_dimension(ncid, 'lon', Im)
  write(*,*) Im,Jm
  
  ! Allocate variable
  allocate(lat(Jm))
  allocate(lon(Im))
  
  call check( nf90_inq_varid(ncid, 'lat', var_id) )
  call check( nf90_get_var(ncid, var_id, lat) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) )
  call check( nf90_get_var(ncid, var_id, lon) )
  call check( nf90_close(ncid) )

#else
!---- Read JMA-MSM GRIB2 file --------------------------------

  GRIB_yyyymmddhh(1:4)=YYYY
  GRIB_yyyymmddhh(5:6)=MM
  GRIB_yyyymmddhh(7:8)=DD

  IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
              GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix
  !Open GRIB file
  call codes_grib_multi_support_on	(	iret	)	
  write(*,*) "OPEN: ", trim( IN_FILE(1) )
  call codes_open_file(ifile, trim( IN_FILE(1) ),'r', iret)
  call codes_grib_new_from_file(ifile,igrib, iret)
!
! ! Get dimension data
  call codes_get(igrib,'Ny', Jm)
  call codes_get(igrib,'Nx', Im)
  write(*,*) Im,Jm
!
  call codes_get(igrib,'distinctLatitudes',lat)
  call codes_get(igrib,'distinctLongitudes',lon)
      
  call codes_release(igrib)
  call codes_close_file(ifile)

  ! Allocate variable
  allocate(values(Im*Jm))
#endif

  allocate(in_data(N_Param, Im, Jm))
      
  write(*,*) "CALC.: weight parameters for interpolation"
  call weight2D_grid(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)
   
      
!---- Create the forcing netCDF file --------------------------------

  FRC_yyyymmdd(2:5)=YYYY
  FRC_yyyymmdd(6:7)=MM
  FRC_yyyymmdd(8:9)=DD  

  FRC_FILE(1) = trim(FRC_prefix)//FRC_yyyymmdd//'_1.nc'
  
  write(*,*) "CREATE: ", trim( FRC_FILE(1) )
  call check( nf90_create(trim( FRC_FILE(1) ), nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xi_rho_dimid) )
  call check( nf90_def_dim(ncid, 'eta_rho',Nyr, eta_rho_dimid) )
  call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid) )
  dimids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)
  call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, var_id) )
  call check( nf90_put_att(ncid, var_id, 'long_name', 'atmospheric forcing time') )
  call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

  DO iparam=1,9
    call check( nf90_def_var(ncid, trim( NC_NAME(iparam) ), NF90_REAL, dimids, var_id) )
    call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(iparam) )) )
    call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(iparam) ) ) )
    call check( nf90_put_att(ncid, var_id, 'time',      'time') )
  END DO
! End define mode.
  call check( nf90_enddef(ncid) )
  call check( nf90_close(ncid) ) 
  
  FRC_FILE(2) = trim(FRC_prefix)//FRC_yyyymmdd//'_2.nc'
  
  write(*,*) "CREATE: ", trim( FRC_FILE(2) )
  call check( nf90_create(trim( FRC_FILE(2) ), nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xi_rho_dimid) )
  call check( nf90_def_dim(ncid, 'eta_rho',Nyr, eta_rho_dimid) )
  call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid) )
  dimids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)
! Define the netCDF variables for the pressure and temperature data.
  call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, var_id) )
  call check( nf90_put_att(ncid, var_id, 'long_name', 'atmospheric forcing time') )
  call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

  DO iparam=10,N_Param
    call check( nf90_def_var(ncid, trim( NC_NAME(iparam) ), NF90_REAL, dimids, var_id) )
    call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(iparam) )) )
    call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(iparam) ) ) )
    call check( nf90_put_att(ncid, var_id, 'time',      'time') )
  END DO
! End define mode.
  call check( nf90_enddef(ncid) )
  call check( nf90_close(ncid) )

      
!---- LOOP set up --------------------------------
  itime = 1
  ihours = 0
  iyear = Syear
  imonth = Smonth
  iday = Sday
  call jd(iyear, imonth, iday, juliandate)

!---- LOOP1 START --------------------------------
  DO

    ihour = mod(ihours,24)
!    juliandate = juliandate + int(ihours/24)
    call cdate( juliandate + int(ihours/24), iyear, imonth, iday )
    ! Check end date
    if(iyear==Eyear .and. imonth==Emonth .and. iday==Eday) then
      write(*,*) "Completed!!!"
      STOP
    endif

    write (YYYY, "(I4.4)") iyear
    write (MM, "(I2.2)") imonth
    write (DD, "(I2.2)") iday ! 1+int((itime-1)*1/24)
    write (hh, "(I2.2)") ihour ! mod((itime-1)*1,24)

#if defined NETCDF_INPUT
    ihours = ihours + 24  !!! Files exist daily interval

    IN_FILE2(1) = IN_FILE(1)  
    IN_FILE2(2) = IN_FILE(2)  
    IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//DD//".nc"
    IN_FILE(2) = trim(MSM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

#else
    ihours = ihours + 3  !!! Files exist 3 hourly interval

    GRIB_yyyymmddhh(1:4)=YYYY
    GRIB_yyyymmddhh(5:6)=MM
    GRIB_yyyymmddhh(7:8)=DD
    GRIB_yyyymmddhh(9:10)=hh

    IN_FILE2(1) = IN_FILE(1)  
    IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
                  GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix

#endif
    ! Check GRIB/nc file
    write(*,*) "CHECK: ", trim( IN_FILE(1) )
    inquire(FILE=trim( IN_FILE(1) ), EXIST=file_exists)
    if( file_exists ) then
      isp=1
      write(*,*) "PASSED"
    else
      IN_FILE(1) = IN_FILE2(1)
      isp=isp+3
      write(*,*) "Not found..."
      if(isp>13) cycle
    endif

!  ---- LOOP2 START --------------------------------        
#if defined NETCDF_INPUT
    DO ifc=1,24
#else
    DO ifc=isp,isp+2
#endif
!  ---- LOOP3.1 START --------------------------------
      DO iparam=1,N_Param
#if defined NETCDF_INPUT
        if(iparam>=10) then !!! for rain (rain fall rate)
          inc=2
        else
          inc=1
        end if
        write(*,*) "OPEN: ", trim( IN_FILE(inc) )
        !Open NetCDF file
        call check( nf90_open(trim( IN_FILE(inc) ), nf90_nowrite, ncid) )

        start3D = (/ 1,  1,  ifc /)
        count3D = (/ Im, Jm, 1   /)
        ! Get variable id
        call check( nf90_inq_varid(ncid, trim( NCIN_NAME(iparam) ), var_id) ) 
        call check( nf90_get_var(ncid, var_id, in_data(iparam,:,:), start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        call check( nf90_close(ncid) )     
        in_data(iparam,:,:)=in_data(iparam,:,:)*sf+off
   
#else
!      ---- Seek message --------------------------------
        write(*,*) "OPEN: ", trim( IN_FILE(1) )
        call codes_open_file(ifile, trim( IN_FILE(1) ),'r', iret)
        call codes_grib_new_from_file(ifile,igrib, iret)
        DO WHILE (iret /= GRIB_END_OF_FILE)
          call codes_get(igrib,'forecastTime',p1)
          call codes_get(igrib,'parameterName', p2)
          if (p1==GRIB_STEP(ifc) .and.             &
              trim(p2)==trim(GRIB_NAME(iparam))  ) exit
          call codes_release(igrib)
          call codes_grib_new_from_file(ifile,igrib, iret)
        END DO
   
        write(*,*) "READ GRIB DATA: ", trim(p2)
        call codes_get(igrib,'validityDate',YYYYMMDD(iparam))
        write(*,*) 'validityDate = ', YYYYMMDD(iparam)
        call codes_get(igrib,'validityTime',hhmm(iparam))
        write(*,*) 'validityTime = ', hhmm(iparam)
  
!        call codes_get(igrib,'dataDate',YYYYMMDD)
        call codes_get(igrib,'values', values)
        call codes_release(igrib)
        call codes_close_file(ifile)      
        
        do i=1, Jm
          istart = 1 + Im*(i-1)
          iend   = Im*i
          in_data(iparam,:,i) = values(istart:iend)
        end do
#endif
      END DO
#if defined NETCDF_INPUT
    ! Set date & time
      call check( nf90_open(trim( IN_FILE(1) ), nf90_nowrite, ncid) )
      start1D = (/ ifc /)
      count1D = (/ 1 /)
      call check( nf90_inq_varid(ncid, 'time', var_id) )  !!!  not Japan time (00:00:00+09:00)
      call check( nf90_get_var(ncid, var_id, time, start=start1D, count=count1D) )
      call check( nf90_close(ncid) )
      call ndays(imonth, iday, iyear, 1, 1, 2000, d_ref_days)
     t = time(1)/24.0d0 + dble(d_ref_days)
#else
    ! Set date & time
      iyear  = YYYYMMDD(1)/10000
      imonth = (YYYYMMDD(1)-iyear*10000)/100
      iday   = YYYYMMDD(1)-iyear*10000-imonth*100
      ihour  = hhmm(1)/100
      imin   = hhmm(1)-100*ihour

    !  ihour  = ihour-1 ! since time for precipitation is set +1 hour

      call ndays(imonth, iday, iyear, Rmonth, Rday, Ryear, idays)
      
      t = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
#endif
!  ---- LOOP3.1 END --------------------------------
    ! for Pair (Pressure)
      in_data(1,:,:) = in_data(1,:,:)*0.01  ! Pa -> millibar (= hPa)    
    ! for Tair
      in_data(4,:,:) = in_data(4,:,:) - 273.15d0  ! K -> degC  
    ! for cloud (cloud fraction)
      in_data(6,:,:) = in_data(6,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data(7,:,:) = in_data(7,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data(8,:,:) = in_data(8,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data(9,:,:) = in_data(9,:,:)*0.01  ! percent -> ratio(0 to 1)    
    ! for rain (Total precipitation rate)
      in_data(10,:,:) = in_data(10,:,:)/3600.0d0  ! kg m-2 h-1 -> kg m-2 s-1  
          
!  ---- LOOP3.2 START --------------------------------
      DO iparam=1,N_Param
      
        write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )
        call interp2D_grid(Im, Jm, in_data(iparam,:,:)                  &
                           , Nxr, Nyr, out_data(iparam,:,:)    &
                           , Id_cont, w_cont)
      END DO
!  ---- LOOP3.2 END --------------------------------

  !!! for U V: change regular Lat Lon to ROMS grid coordinat vectors 
      do i=0,L
        do j=0,M
          u = out_data(2,i,j)*cosAu(i,j)+out_data(3,i,j)*sinAu(i,j)
          v =-out_data(2,i,j)*sinAv(i,j)+out_data(3,i,j)*cosAv(i,j)
          out_data(2,i,j) = u
          out_data(3,i,j) = v
        enddo
      enddo

!  ---- LOOP3.3 START --------------------------------
      time(1) = t
      write(*,*) time(1),TIME_ATT
      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d( 'time', trim( FRC_FILE(1) )                  &
            , 1, time, start1D, count1D )

      time(1) = t+0.5d0/24.0d0
      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d( 'time', trim( FRC_FILE(2) )                  &
            , 1, time, start1D, count1D )

      DO iparam=1,N_Param
        if(iparam>=10) then !!! for rain (rain fall rate)
          inc=2
        else
          inc=1
        end if
  
        start3D = (/ 1,  1,  itime /)
        count3D = (/ Nxr, Nyr, 1 /)     
        call writeNetCDF_3d(trim( NC_NAME(iparam) ), trim( FRC_FILE(inc))   &
            , Nxr, Nyr, 1, out_data(iparam,:,:), start3D, count3D )
        
      END DO
!  ---- LOOP3.3 END --------------------------------         
      itime = itime + 1
    END DO
! ---- LOOP2 END --------------------------------
  END DO
!---- LOOP1 END --------------------------------
  
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM frcJMAMSM
      
