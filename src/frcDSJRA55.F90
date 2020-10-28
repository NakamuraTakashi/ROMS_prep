
!!!=== Copyright (c) 2018-2020 Takashi NAKAMURA  =====

PROGRAM frcDSJRA55
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
  character(256) :: GRIB_LL
  character(256) :: FRC_prefix
  integer :: Ryear, Rmonth, Rday
  integer :: mode
!----------------------------------------------------------------------
  integer, parameter :: N_Param = 7
  
  character(len=*), parameter :: GRIB_FCST_SURF_prefix  = "fcst_surf"
  character(len=*), parameter :: GRIB_FCST_PHY2M_prefix = "fcst_phy2m"
  character(11) :: GRIB_suffix  = ".1958010100"
  character(10) :: FRC_suffix   = "_195801.nc"
  
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: GRIB_FILE
  character(256) :: FRC_FILE(N_Param)
  
  integer, parameter :: GRIB_NUM(2,N_Param) = reshape ((/   &
! parameter  parameter
! Category  ,Number
     2      ,2        &  ! 10 metre U wind component
    ,2      ,3        &  ! 10 metre V wind component
    ,0      ,0        &  ! Air temperature
    ,0      ,7        &  ! Dewpoint depression (or deficit) (K)
    ,3      ,0        &  ! Surface pressure
    ,6      ,1        &  ! Total Cloud Cover
    ,1      ,52       &  ! Total precipitation rate  (kg m-2 s-1)
    /), (/2, N_Param/))
  character(256) :: GRIB_NAME(N_Param) = (/   &
     "10u    "                                &
    ,"10v    "                                &
    ,"t      "                                &
    ,"unknown"                                &
    ,"msl    "                                &
    ,"tcc    "                                &
    ,"mtpf   "                                &
    /)
 
  character(256) :: NC_NAME(N_Param) = (/     &
     "Uwind"                                  &
    ,"Vwind"                                  &
    ,"Tair "                                  &
    ,"Qair "                                  &
    ,"Pair "                                  &
    ,"cloud"                                  &
    ,"rain "                                  &
    /)
  character(256) :: NC_LNAME(N_Param) = (/    &
     "surface u-wind component     "          &
    ,"surface v-wind component     "          &
    ,"surface air temperature      "          &
    ,"surface air relative humidity"          &
    ,"surface air pressure         "          &
    ,"cloud fraction               "          &
    ,"rain fall rate               "          &
    /)
  character(256) :: NC_UNIT(N_Param) = (/     &
     "meter second-1           "              &
    ,"meter second-1           "              &
    ,"Celsius                  "              &
    ,"percentage               "              &
    ,"millibar                 "              &
    ,"0 to 1                   "              &
    ,"kilogram meter-2 second-1"              &
    /)
  real(8), parameter :: pi = 3.141592653589793d0
  real(8), allocatable :: lat_rho(:, :)
  real(8), allocatable :: lon_rho(:, :)
  real(8), allocatable :: cosAu(:, :)
  real(8), allocatable :: sinAu(:, :)
  real(8), allocatable :: cosAv(:, :)
  real(8), allocatable :: sinAv(:, :)
  
  real(8), allocatable :: out_data(:,:,:) ! output forcing data
       
  integer :: Im, Jm
  real(8), allocatable :: lat(:, :), lon(:, :)
  real(8), allocatable :: cosAx(:, :), sinAx(:, :)
  real(8), allocatable :: cosAy(:, :), sinAy(:, :)
  real(8), allocatable :: grib_data(:,:, :)
  real(8) :: t
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
  integer :: N_xi_u, N_eta_u
  integer :: N_xi_v, N_eta_v
  
  integer :: dimids(3)
  integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
  
  integer :: iparam
  integer :: ifile,idx,iret,igrib
  integer :: istart, iend
  integer :: YYYYMMDD, hhmm
  real :: bd
  real(8), allocatable :: values(:)
  integer :: p1,p2
  character(10) :: p3
  real(8) :: d_lat, d_lon
  real(8) :: u, v
  real(8) :: dpT, sat_VP, VP
!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/intpmode/mode
  namelist/frc_dsjra55/GRIB_FCST_SURF_dir
  namelist/frc_dsjra55/GRIB_FCST_PHY2M_dir
  namelist/frc_dsjra55/GRIB_LL
  namelist/frc_dsjra55/FRC_prefix
  ! Read parameters in namelist file
  
  read (*, nml=grd)
  read (*, nml=sdate)
  read (*, nml=refdate)
  read (*, nml=intpmode)
  read (*, nml=frc_dsjra55)


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
  call get_dimension(ncid, 'xi_rho',  N_xi_rho)
  call get_dimension(ncid, 'eta_rho', N_eta_rho)
  call get_dimension(ncid, 'xi_u',    N_xi_u)
  call get_dimension(ncid, 'eta_u',   N_eta_u)
  call get_dimension(ncid, 'xi_v',    N_xi_v)
  call get_dimension(ncid, 'eta_v',   N_eta_v)
  
  allocate(lat_rho(N_xi_rho, N_eta_rho))
  allocate(lon_rho(N_xi_rho, N_eta_rho))
  allocate(cosAu(N_xi_u, N_eta_u))
  allocate(sinAu(N_xi_u, N_eta_u))
  allocate(cosAv(N_xi_v, N_eta_v))
  allocate(sinAv(N_xi_v, N_eta_v))
  allocate(out_data(N_Param, N_xi_rho, N_eta_rho))
  
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lat_rho) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lon_rho) )
  
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
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
      
!---- Read DSJRA-55 GRIB2 file --------------------------------

  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  GRIB_suffix(2:5)=YYYY
  GRIB_suffix(6:7)=MM
  GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
  !Open GRIB file
  call codes_grib_multi_support_on	(	iret	)	
  write(*,*) "OPEN: ", trim( GRIB_FILE )
  call codes_open_file(ifile, trim( GRIB_FILE ),'r', iret)
  call codes_grib_new_from_file(ifile,igrib, iret)
!
! ! Get dimension data
  call codes_get(igrib,'Ny', Jm)
  call codes_get(igrib,'Nx', Im)
  write(*,*) Im,Jm
!          
  ! Allocate variable
  allocate(values(Im*Jm))
  allocate(lat(Im,Jm))
  allocate(lon(Im,Jm))
  allocate(cosAx(Im,Jm))
  allocate(sinAx(Im,Jm))
  allocate(cosAy(Im,Jm))
  allocate(sinAy(Im,Jm))
  allocate(grib_data(N_Param, Im, Jm))
      
!  open(unit=21, file='lat.txt')
!  open(unit=22, file='lon.txt')
      
  call codes_release(igrib)
  call codes_close_file(ifile)

!  ---- Get Lat Lon coordinates from Binary file --------------
!
  open(unit=20, file=trim(GRIB_LL), action='read',                 &
       & form='unformatted', access='direct', recl=4,        &
       & CONVERT='BIG_ENDIAN', status='old')
  do j=1, Jm
    do i=1,Im
      read(20, rec=i+Im*(j-1)) bd
      lat(i,j) =dble(bd)
!      write(*,*) i+Im*(j-1), lat(i,j)
!       write(*,*) lat(i,j)
    end do
!    write(21,*) lat(:,j)
  end do
  do j=1, Jm
    do i=1,Im
      read(20, rec=Im*Jm+i+Im*(j-1)) bd
      lon(i,j) = dble(bd)
!      write(*,*) i+Im*(j-1), lat(i,j)
    end do
!    write(22,*) lon(:,j)
  end do
  close(20)
!  close(21)
!  close(22)
  write(*,*) 'NW corner:', lat(1,1),  lon(1,1)
  write(*,*) 'SW corner:', lat(1,Jm), lon(1,Jm)
  write(*,*) 'SE corner:', lat(Im,1), lon(Im,1)
  write(*,*) 'NE corner:', lat(Im,Jm),lon(Im,Jm)

  allocate(ID_cont(8, N_xi_rho*N_eta_rho))
  allocate(w_cont(3, N_xi_rho*N_eta_rho))
  
  write(*,*) "CALC.: weight parameters for interpolation"
  call weight2D_grid2(Im,Jm,lon,lat,N_xi_rho,N_eta_rho,lon_rho,lat_rho,ID_cont,w_cont)
  write(*,*) ID_cont(:,1)
  write(*,*) ID_cont(:,N_xi_rho*N_eta_rho)
  write(*,*) w_cont(:,1)
  write(*,*) w_cont(:,N_xi_rho*N_eta_rho)
      
!  ---- Calc. rotation angle --------------
  do j=1, Jm
    do i=1,Im-1
      d_lat=lat(i+1,j)-lat(i,j)
      d_lon=lon(i+1,j)-lon(i,j)
      d_lon=d_lon*cos(lat(i,j)/180.0d0*PI)
      cosAx(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAx(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
    end do
  end do
  cosAx(Im,:) = cosAx(Im-1,:)
  sinAx(Im,:) = sinAx(Im-1,:)
  
  do j=1, Jm-1
    do i=1,Im
      d_lat=lat(i,j)-lat(i,j+1)
      d_lon=lon(i,j+1)-lon(i,j)
      d_lon=d_lon*cos(lat(i,j)/180.0d0*PI)
      cosAy(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAy(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
    end do
  end do
  cosAy(:,Jm) = cosAy(:,Jm-1)
  sinAy(:,Jm) = sinAy(:,Jm-1)
      
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
  
  DO
  
    write (DD, "(I2.2)") 1+int((itime-1)*1/24)
    write (hh, "(I2.2)") mod((itime-1)*1,24)
    GRIB_suffix(8:9)=DD
    GRIB_suffix(10:11)=hh
  
    GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
    
    !Open GRIB file
    write(*,*) "OPEN: ", trim( GRIB_FILE )
    call codes_open_file(ifile, trim( GRIB_FILE ),'r', iret)
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
    
    t = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
        

!  ---- LOOP2.1 START --------------------------------
    DO iparam=1,N_Param
    
      if(iparam==7) then !!! for rain (rain fall rate)
        GRIB_FILE = trim(GRIB_FCST_PHY2M_dir)//GRIB_FCST_PHY2M_prefix//GRIB_suffix
      else
        GRIB_FILE = trim(GRIB_FCST_SURF_dir)//GRIB_FCST_SURF_prefix//GRIB_suffix
      end if
      
      !Open GRIB file
      write(*,*) "OPEN: ", trim( GRIB_FILE )
      call codes_open_file(ifile, trim( GRIB_FILE ),'r', iret)
      call codes_grib_new_from_file(ifile,igrib, iret)
          
!    ---- Seek message --------------------------------
      DO WHILE (iret /= GRIB_END_OF_FILE)
        call codes_get(igrib,'parameterCategory',p1)
        call codes_get(igrib,'parameterNumber',  p2)
        call codes_get(igrib,'shortName',        p3)
        if (p1==GRIB_NUM(1,iparam) .and.             &
            p2==GRIB_NUM(2,iparam) .and.             &
            p3==GRIB_NAME(iparam)        ) exit
        call codes_release(igrib)
        call codes_grib_new_from_file(ifile,igrib, iret)
      END DO

!    ---- Get GRIB data --------------------------------
      write(*,*) "READ GRIB DATA: ", trim( GRIB_NAME(iparam) )
      call codes_get(igrib,'values', values)
      
      call codes_release(igrib)
      call codes_close_file(ifile)
      
      do i=1, Jm
        istart = 1 + Im*(i-1)
        iend   = Im*i
        grib_data(iparam,:,i) = values(istart:iend)
      end do
    END DO
        

!  ---- LOOP2.1 END --------------------------------

  !!! for U V: change DSJRA55 Lambert conformal to regular Lat Lon coordinat vectors 
    do j=1, Jm
      do i=1,Im
        u = grib_data(1,i,j)*cosAx(i,j)-grib_data(2,i,j)*sinAy(i,j)
        v = grib_data(1,i,j)*sinAx(i,j)+grib_data(2,i,j)*cosAy(i,j)
        grib_data(1,i,j) = u
        grib_data(2,i,j) = v
      end do
    end do
      
  !!! for Tair
    grib_data(3,:,:) = grib_data(3,:,:) - 273.15d0  ! K -> degC
    
  !!! for Qair: convert Dewpoint depression (K) to Relative humidity (%)
    do j=1, Jm
      do i=1,Im
        ! Dewpoint (oC)
        dpT = grib_data(3,i,j)-grib_data(4,i,j) ! (oC)
        ! Saturation vapor pressure (hPa)
        sat_VP=6.1078d0*10.0d0**(7.5d0*grib_data(3,i,j)/(237.3d0+grib_data(3,i,j)))
        ! Vapor pressure (hPa)
        VP    =6.1078d0*10.0d0**(7.5d0*dpT/(237.3d0+dpT))
        ! Relative humidity (%)
        grib_data(4,i,j) = VP/sat_VP*100.0d0 ! (%)
      end do
    end do
  !!! for Pair (Pressure)
    grib_data(5,:,:) = grib_data(5,:,:)*0.01  ! Pa -> millibar (= hPa)
    
  !!! for cloud (cloud fraction)
    grib_data(6,:,:) = grib_data(6,:,:)*0.01  ! percent -> ratio(0 to 1)
     
          
!  ---- LOOP2.2 START --------------------------------
    DO iparam=1,N_Param
    
      write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )
      call interp2D_grid2(Im, Jm, grib_data(iparam,:,:)                       &
                        , N_xi_rho, N_eta_rho, out_data(iparam,:,:)           &
                        , Id_cont, w_cont)
    END DO
!  ---- LOOP2.2 END --------------------------------

  !!! for U V: change regular Lat Lon to ROMS grid coordinat vectors 
      do i=1,N_xi_rho
        do j=1,N_eta_rho
          u = out_data(1,i,j)*cosAu(i,j)+out_data(2,i,j)*sinAu(i,j)
          v =-out_data(1,i,j)*sinAv(i,j)+out_data(2,i,j)*cosAv(i,j)
          out_data(1,i,j) = u
          out_data(2,i,j) = v
        enddo
      enddo

!  ---- LOOP2.3 START --------------------------------
    DO iparam=1,N_Param
      if(iparam==7) then !!! for rain (rain fall rate)
        time(1) = t+0.5d0/24.0d0  !!! + 0.5 hours
      else
        time(1) = t
      end if
      write(*,*) time(1),TIME_ATT
      
      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d(                    &
!           input parameters
              'time'                          &
            , trim( FRC_FILE(iparam) )        &
            , 1                               &
            , time                            &
            , start1D, count1D                &
            )
      start3D = (/ 1,  1,  itime /)
      count3D = (/ N_xi_rho, N_eta_rho, 1 /)
      
      call writeNetCDF_3d(                    &
!         input parameters
            trim( NC_NAME(iparam) )           &
          , trim( FRC_FILE(iparam))           &
          , N_xi_rho, N_eta_rho, 1            &
          , out_data(iparam,:,:)              &
          , start3D, count3D                  &
          )
      
    END DO
!  ---- LOOP2.3 END --------------------------------
        
    itime = itime + 1
  END DO
!---- LOOP1 END --------------------------------

  deallocate(values)
!  deallocate(lat)
!  deallocate(lon)
  deallocate(grib_data)
  
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM frcDSJRA55
      
