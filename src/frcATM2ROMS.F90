
!!!=== Copyright (c) 2014-2021 Takashi NAKAMURA  =====

#if defined JMA_MSM
# undef BULK_FLUX
# if defined NETCDF_INPUT
#  undef SWRAD
# endif
#elif defined DSJRA55 || defined JRA55
# undef SWRAD
#endif

PROGRAM frcATM2ROMS
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
  integer :: Ryear, Rmonth, Rday
!----------------------------------------------------------------------
#if defined JMA_MSM
!--- JMA_MSM parameter setting -----------------
  character(256) :: MSM_dir
  character(256) :: FRC_prefix
  integer :: jd_msmnew

# if defined SWRAD
  integer, parameter :: N_InPar  = 11
  integer, parameter :: N_OutPar = 11
# else
  integer, parameter :: N_InPar  = 10
  integer, parameter :: N_OutPar = 10
# endif

# if defined NETCDF_INPUT

  character(256) :: NCIN_NAME(N_InPar) = (/   &
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
# else

  character(len=*), parameter :: GRIB_prefix  = "Z__C_RJTD_"
  character(len=*), parameter :: GRIB_suffix  =   &
        "0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"

  character(256) :: GRIB_NAME(N_InPar) = (/   &
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
# endif

#elif defined DSJRA55
!--- DSJRA55 parameter setting -----------------
  character(256) :: DSJRA55_dir
  character(256) :: FRC_prefix
  character(256) :: LL_FILE
  real :: bd
# if defined BULK_FLUX
  integer, parameter :: N_InPar  = 16
  integer, parameter :: N_OutPar = 15
# else
  integer, parameter :: N_InPar  = 10
  integer, parameter :: N_OutPar = 10
# endif
  integer, parameter :: GRIB_NUM(3,N_InPar) = reshape ((/   &
! parameter  parameter  Type of First
! Category   Number     Fixed Surface           
     3         ,0         ,101        &  ! Mean sea level pressure
    ,2         ,2         ,103        &  ! 10 metre U wind component
    ,2         ,3         ,103        &  ! 10 metre V wind component
    ,0         ,0         ,103        &  ! Air temperature
    ,0         ,7         ,103        &  ! Dewpoint depression (or deficit) (K)
    ,6         ,3         ,1          &  ! Low Cloud Cover
    ,6         ,4         ,1          &  ! Medium Cloud Cover
    ,6         ,5         ,1          &  ! High Cloud Cover
    ,6         ,1         ,1          &  ! Total Cloud Cover
    ,1         ,52        ,1          &  ! Total precipitation rate  (kg m-2 s-1)
# if defined BULK_FLUX
    ,0         ,10        ,1          &  ! Latent heat net flux  (W m-2)
    ,0         ,11        ,1          &  ! Sensible heat net flux  (W m-2)
    ,4         ,7         ,1          &  ! Downward short-wave radiation flux  (W m-2)
    ,4         ,8         ,1          &  ! Upward short-wave radiation flux  (W m-2)
    ,5         ,3         ,1          &  ! Downward long-wave radiation flux  (W m-2)
    ,5         ,4         ,1          &  ! Upward long-wave radiation flux  (W m-2)
# endif
    /), (/3, N_InPar/))

#elif defined JRA55
!--- JRA55 parameter setting -----------------
  character(256) :: JRA55_dir
  character(256) :: FRC_prefix

# if defined BULK_FLUX
  integer, parameter :: N_InPar  = 16
  integer, parameter :: N_OutPar = 15
# else
  integer, parameter :: N_InPar  = 10
  integer, parameter :: N_OutPar = 10
# endif
  integer, parameter :: GRIB_NUM(N_InPar) = (/   &
! Indicator Of 
! Parameter              
     2       &  ! Mean sea level pressure
    ,33      &  ! 10 metre U wind component
    ,34      &  ! 10 metre V wind component
    ,11      &  ! Air temperature
    ,52      &  ! Relative humidity
    ,73      &  ! Low Cloud Cover
    ,74      &  ! Medium Cloud Cover
    ,75      &  ! High Cloud Cover
    ,71      &  ! Total Cloud Cover
    ,61      &  ! Total precipitation rate  (kg m-2 s-1)
# if defined BULK_FLUX
    ,121     &  ! Latent heat net flux  (W m-2)
    ,122     &  ! Sensible heat net flux  (W m-2)
    ,204     &  ! Downward short-wave radiation flux  (W m-2)
    ,211     &  ! Upward short-wave radiation flux  (W m-2)
    ,205     &  ! Downward long-wave radiation flux  (W m-2)
    ,212     &  ! Upward long-wave radiation flux  (W m-2)
# endif
    /)
#endif
  character(10) :: GRIB_yyyymmddhh = "2002070100"

  integer :: ifile,idx,iret,igrib
  integer :: istart, iend
  integer :: YYYYMMDD(N_InPar), hhmm(N_InPar)
  real(8), allocatable :: values(:)
  integer :: p1,p2,p3
  character(256) :: p4

  character(9) :: FRC_yyyymmdd = "_20020701" 
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: FRC_FILE(2)

  character(256) :: NC_NAME(N_OutPar) = (/    &
     "Pair      "                             &
    ,"Uwind     "                             &
    ,"Vwind     "                             &
    ,"Tair      "                             &
    ,"Qair      "                             &
    ,"lcloud    "                             &
    ,"mcloud    "                             &
    ,"hcloud    "                             &
    ,"cloud     "                             &
    ,"rain      "                             &
#if defined SWRAD
    ,"swrad     "                             &
#elif defined BULK_FLUX
    ,"latent    "                             &
    ,"sensible  "                             &
    ,"swrad     "                             &
    ,"lwrad_down"                             &
    ,"lwrad     "                             &
#endif
    /)
  character(256) :: NC_LNAME(N_OutPar) = (/   &
     "surface air pressure               "    &
    ,"surface u-wind component           "    &
    ,"surface v-wind component           "    &
    ,"surface air temperature            "    &
    ,"surface air relative humidity      "    &
    ,"low cloud fraction                 "    &
    ,"medium cloud fraction              "    &
    ,"high cloud fraction                "    &
    ,"cloud fraction                     "    &
    ,"rain fall rate                     "    &
#if defined SWRAD
    ,"solar shortwave radiation flux     "    &
#elif defined BULK_FLUX
    ,"net latent heat flux               "    &
    ,"net sensible heat flux             "    &
    ,"solar shortwave radiation flux     "    &
    ,"downwelling longwave radiation flux"    &
    ,"net longwave radiation flux        "    &
#endif
    /)
  character(256) :: NC_UNIT(N_OutPar) = (/    &
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
#elif defined BULK_FLUX
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
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
  
  real(8), allocatable :: in_data(:,:, :), in_data2(:,:, :)
  real(8), allocatable :: out_data(:,:,:) ! output forcing data
       
  integer :: Im, Jm
#if defined DSJRA55
  real(8), allocatable :: lat(:, :), lon(:, :)
  real(8), allocatable :: cosAx(:, :), sinAx(:, :)
  real(8), allocatable :: cosAy(:, :), sinAy(:, :)
#else
  real(8), allocatable :: lat(:), lon(:)
#endif
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
  integer :: idays,ihours,ijdate,Sjdate
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
  real(8) :: dpT, sat_VP, VP
  logical :: file_exists
!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
#if defined JMA_MSM
  namelist/frc_jmamsm/MSM_dir
  namelist/frc_jmamsm/FRC_prefix
#elif defined DSJRA55
  namelist/frc_dsjra55/DSJRA55_dir
  namelist/frc_dsjra55/FRC_prefix
#elif defined JRA55
  namelist/frc_jra55/JRA55_dir
  namelist/frc_jra55/FRC_prefix
#endif
  ! Read parameters in namelist file
  
  read (5, nml=grd)
  rewind(5)
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
  rewind(5)
  read (5, nml=refdate)
#if defined JMA_MSM
  rewind(5)
  read (5, nml=frc_jmamsm)
#elif defined DSJRA55
  rewind(5)
  read (5, nml=frc_dsjra55)
#elif defined JRA55
  rewind(5)
  read (5, nml=frc_jra55)
#endif
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
  L = Nxr-1
  M = Nyr-1
  allocate( latr(0:L, 0:M) )
  allocate( lonr(0:L, 0:M) )
  allocate( cosAu(0:L, 0:M) )
  allocate( sinAu(0:L, 0:M) )
  allocate( cosAv(0:L, 0:M) )
  allocate( sinAv(0:L, 0:M) )
  allocate(out_data(N_OutPar, 0:L, 0:M))
!#if defined DSJRA55
!  allocate(ID_cont(8, Nxr*Nyr))
!  allocate(w_cont(3, Nxr*Nyr))
!#else
  allocate(ID_cont(4, Nxr*Nyr))
  allocate(w_cont(4, Nxr*Nyr))
!#endif
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
  GRIB_yyyymmddhh = YYYY//MM//DD//"00"
# if defined JMA_MSM
!---- Set JMA-MSM GRIB2 file --------------------------------
  IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
              GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix

# elif defined DSJRA55
!---- Set DSJRA55 GRIB2 file --------------------------------
  IN_FILE(1) = trim(DSJRA55_dir)//"Hist/Daily/fcst_surf/"//YYYY//MM// &
                 "/fcst_surf."//GRIB_yyyymmddhh
  IN_FILE(2) = trim(DSJRA55_dir)//"Hist/Daily/fcst_phy2m/"//YYYY//MM// &
                 "/fcst_phy2m."//GRIB_yyyymmddhh

# elif defined JRA55
!---- Read JRA55 GRIB2 file --------------------------------
  IN_FILE(1) = trim(JRA55_dir)//"Hist/Daily/fcst_surf125/"//YYYY//MM// &
                 "/fcst_surf125."//GRIB_yyyymmddhh
  IN_FILE(2) = trim(JRA55_dir)//"Hist/Daily/fcst_phy2m125/"//YYYY//MM// &
                 "/fcst_phy2m125."//GRIB_yyyymmddhh

# endif

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
! Allocate variable
  allocate(values(Im*Jm))
# if defined DSJRA55
  allocate(lat(Im,Jm))
  allocate(lon(Im,Jm))
  allocate(cosAx(Im,Jm))
  allocate(sinAx(Im,Jm))
  allocate(cosAy(Im,Jm))
  allocate(sinAy(Im,Jm))
# else
  call codes_get(igrib,'distinctLatitudes',lat)
  call codes_get(igrib,'distinctLongitudes',lon)
# endif      
  call codes_release(igrib)
  call codes_close_file(ifile)
#endif

  allocate(in_data(N_InPar, Im, Jm))
  allocate(in_data2(N_OutPar, Im, Jm))

!---- Calculate weight parameters for interpolation --------------------------------
#if defined DSJRA55
!  ---- Get Lat Lon coordinates from Binary file --------------
!
  LL_FILE = trim(DSJRA55_dir)//"Consts/Lambert5km_latlon.dat"
  open(unit=20, file=LL_FILE, action='read',               &
       form='unformatted', access='direct', recl=4,        &
       CONVERT='BIG_ENDIAN', status='old')
  do j=1, Jm
    do i=1,Im
      read(20, rec=i+Im*(j-1)) bd
      lat(i,j) =dble(bd)
    end do
  end do
  do j=1, Jm
    do i=1,Im
      read(20, rec=Im*Jm+i+Im*(j-1)) bd
      lon(i,j) = dble(bd)
    end do
  end do
  close(20)

  write(*,*) 'NW corner:', lat(1,1),  lon(1,1)
  write(*,*) 'SW corner:', lat(1,Jm), lon(1,Jm)
  write(*,*) 'SE corner:', lat(Im,1), lon(Im,1)
  write(*,*) 'NE corner:', lat(Im,Jm),lon(Im,Jm)
   
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

!  ---- Calc. weight parameters for interpolation --------------  
  write(*,*) "CALC.: weight parameters for interpolation"
  call weight2D_grid2(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)
#else
!  ---- Calc. weight parameters for interpolation --------------  
  write(*,*) "CALC.: weight parameters for interpolation"
  call weight2D_grid(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)  
#endif      
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

  DO iparam=10,N_OutPar
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
  call jd(iyear, imonth, iday, Sjdate)
#if defined JMA_MSM
  call jd(2006, 3, 1, jd_msmnew)
#endif
!---- LOOP1 START --------------------------------
  DO

    ihour = mod(ihours,24)
    ijdate = Sjdate + int(ihours/24)
    call cdate( ijdate, iyear, imonth, iday )
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
    GRIB_yyyymmddhh = YYYY//MM//DD//hh

# if defined JMA_MSM
    ihours = ihours + 3  !!! Files exist 3 hourly interval

    IN_FILE2(1) = IN_FILE(1)  
    IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
                  GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix

# elif defined DSJRA55
!---- Read DSJRA55 GRIB2 file --------------------------------
    ihours = ihours + 1  !!! Files exist 1 hourly interval

    IN_FILE(1) = trim(DSJRA55_dir)//"Hist/Daily/fcst_surf/"//YYYY//MM// &
                   "/fcst_surf."//GRIB_yyyymmddhh
    IN_FILE(2) = trim(DSJRA55_dir)//"Hist/Daily/fcst_phy2m/"//YYYY//MM// &
                   "/fcst_phy2m."//GRIB_yyyymmddhh

# elif defined JRA55
!---- Read JRA55 GRIB2 file --------------------------------
    ihours = ihours + 3  !!! Files exist 3 hourly interval

    IN_FILE(1) = trim(JRA55_dir)//"Hist/Daily/fcst_surf125/"//YYYY//MM// &
                   "/fcst_surf125."//GRIB_yyyymmddhh
    IN_FILE(2) = trim(JRA55_dir)//"Hist/Daily/fcst_phy2m125/"//YYYY//MM// &
                   "/fcst_phy2m125."//GRIB_yyyymmddhh

# endif
#endif

    ! Check GRIB/nc file
    write(*,*) "CHECK: ", trim( IN_FILE(1) )
    inquire(FILE=trim( IN_FILE(1) ), EXIST=file_exists)
#if defined JMA_MSM && !defined NETCDF_INPUT
    if( file_exists ) then
      isp=1
      write(*,*) "PASSED"
    else
      IN_FILE(1) = IN_FILE2(1)
      isp=isp+3
      write(*,*) "Not found..."
      if(isp>13) cycle
    endif
#else
    if( file_exists ) then
      write(*,*) "PASSED"
    else
      write(*,*) "Not found..."
      cycle
    endif
#endif
!  ---- LOOP2 START --------------------------------
#if defined JMA_MSM       
# if defined NETCDF_INPUT
    DO ifc=1,24
# else
    DO ifc=isp,isp+2
# endif
#else
    DO ifc=1,1
#endif
!  ---- LOOP3.1 START --------------------------------
      DO iparam=1,N_InPar
#if defined JMA_MSM
        if(iparam==9 .and. ijdate<jd_msmnew) cycle
#endif
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
        write(*,*) "READ: ", trim( NCIN_NAME(iparam) )
        call check( nf90_inq_varid(ncid, trim( NCIN_NAME(iparam) ), var_id) ) 
        if(ijdate>=jd_msmnew .or.iparam<=5 .or.iparam==10 ) then
          call check( nf90_get_var(ncid, var_id, in_data(iparam,:,:), start=start3D, count=count3D) )
          call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
          call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
          in_data(iparam,:,:)=in_data(iparam,:,:)*sf+off
        else
          call check( nf90_get_var(ncid, var_id, in_data(iparam,:,:), start=start3D, count=count3D) )
        endif
        call check( nf90_close(ncid) )     
   
#else
!      ---- Seek message --------------------------------
#if defined JMA_MSM
        inc=1
#else
        if(iparam>=10) then !!! for rain (rain fall rate)
          inc=2
        else
          inc=1
        end if
#endif
        write(*,*) "OPEN: ", trim( IN_FILE(inc) )
        call codes_open_file(ifile, trim( IN_FILE(inc) ),'r', iret)
        call codes_grib_new_from_file(ifile,igrib, iret)

        DO WHILE (iret /= GRIB_END_OF_FILE)
#if defined JMA_MSM
          call codes_get(igrib,'forecastTime',p1)
          call codes_get(igrib,'parameterName', p4)
          if (p1==GRIB_STEP(ifc) .and.             &
              trim(p4)==trim(GRIB_NAME(iparam))  ) exit
#elif defined DSJRA55
          call codes_get(igrib,'parameterCategory',p1)
          call codes_get(igrib,'parameterNumber', p2)
          call codes_get(igrib,'typeOfFirstFixedSurface', p3)
          call codes_get(igrib,'parameterName', p4)
          if (p1==GRIB_NUM(1,iparam) .and. p2==GRIB_NUM(2,iparam) .and. &
              p3==GRIB_NUM(3,iparam)   ) exit
#elif defined JRA55
          call codes_get(igrib,'indicatorOfParameter',p1)
          call codes_get(igrib,'parameterName', p4)
          if (p1==GRIB_NUM(iparam) ) exit
#endif
          call codes_release(igrib)
          call codes_grib_new_from_file(ifile,igrib, iret)
        END DO
   
        write(*,*) "READ GRIB DATA: ", trim(p4)
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
#if defined JMA_MSM
    ! for Pair (Pressure)
      if(ijdate<jd_msmnew) then
        in_data2(1,:,:) = in_data(1,:,:)  ! jPa -> millibar (= hPa) 
      else   
        in_data2(1,:,:) = in_data(1,:,:)*0.01  ! Pa -> millibar (= hPa) 
      endif
    ! for U V
      in_data2(2,:,:) = in_data(2,:,:)
      in_data2(3,:,:) = in_data(3,:,:)
    ! for Tair
      if(ijdate<jd_msmnew) then
        in_data2(4,:,:) = in_data(4,:,:)   ! degC -> degC 
      else   
        in_data2(4,:,:) = in_data(4,:,:) - 273.15d0  ! K -> degC
      endif
    ! for Qait (Relative humidity)
      in_data2(5,:,:) = in_data(5,:,:)
    ! for cloud (cloud fraction) 
      if(ijdate<jd_msmnew) then
      ! no total cloud cover data
        in_data2(6,:,:) = in_data(6,:,:)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(7,:,:) = in_data(7,:,:)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(8,:,:) = in_data(8,:,:)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(9,:,:) = 1.0d0 - (1.0-in_data2(6,:,:))*(1.0-in_data2(7,:,:)) &
                                 *(1.0-in_data2(8,:,:)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else   
        in_data2(6,:,:) = in_data(6,:,:)*0.01  ! percent -> ratio(0 to 1)
        in_data2(7,:,:) = in_data(7,:,:)*0.01  ! percent -> ratio(0 to 1)
        in_data2(8,:,:) = in_data(8,:,:)*0.01  ! percent -> ratio(0 to 1)
        in_data2(9,:,:) = in_data(9,:,:)*0.01  ! percent -> ratio(0 to 1)
      endif
    ! for rain (Total precipitation rate)
      in_data2(10,:,:) = in_data(10,:,:)/3600.0d0  ! kg m-2 h-1 -> kg m-2 s-1 
# if defined SWRAD
    ! for short-wave radiation
      in_data2(11,:,:) = in_data(11,:,:)  ! W m-2 -> W m-2
# endif 
! --------------------------------------------
#elif defined DSJRA55
    ! for Pair (Pressure)
      in_data2(1,:,:) = in_data(1,:,:)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V: change DSJRA55 Lambert conformal to regular Lat Lon coordinat vectors 
      do j=1, Jm
        do i=1,Im
          in_data2(2,i,j) = in_data(2,i,j)*cosAx(i,j)-in_data(3,i,j)*sinAy(i,j)
          in_data2(3,i,j) = in_data(2,i,j)*sinAx(i,j)+in_data(3,i,j)*cosAy(i,j)
        end do
      end do
    ! for Tair
      in_data2(4,:,:) = in_data(4,:,:) - 273.15d0  ! K -> degC
    ! for Qair: convert Dewpoint depression (K) to Relative humidity (%)
      do j=1, Jm
        do i=1,Im
          ! Dewpoint (oC)
          dpT = in_data(4,i,j)-in_data(5,i,j) ! (oC)
          ! Saturation vapor pressure (hPa) !!! Check equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          sat_VP=6.1078d0*10.0d0**(7.5d0*in_data2(4,i,j)/in_data(4,i,j))
          ! Vapor pressure (hPa)
          VP    =6.1078d0*10.0d0**(7.5d0*dpT/(237.3d0+dpT))
          ! Relative humidity (%)
          in_data2(5,i,j) = VP/sat_VP*100.0d0 ! (%)
        end do
      end do
    ! for cloud (cloud fraction)
      in_data2(6,:,:) = in_data(6,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(7,:,:) = in_data(7,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(8,:,:) = in_data(8,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(9,:,:) = in_data(9,:,:)*0.01  ! percent -> ratio(0 to 1)
    ! for rain (Total precipitation rate)
      in_data2(10,:,:) = in_data(10,:,:)  ! kg m-2 s-1 -> kg m-2 s-1     
# if defined BULK_FLUX
    ! for bulk flux parameters, heating -> positive?, cooling -> negattive?
      in_data2(11,:,:) = -in_data(11,:,:)  ! W m-2 -> W m-2
      in_data2(12,:,:) = -in_data(12,:,:)  ! W m-2 -> W m-2
      in_data2(13,:,:) = in_data(13,:,:) !-in_data(14,:,:) ! solar radiation = downward (?)!!!!!!!!!!!!!!!!
      in_data2(14,:,:) = in_data(15,:,:) ! downward long-wave (?)!!!!!!!!!!!!!!!!
      in_data2(15,:,:) = in_data(15,:,:)-in_data(16,:,:) ! net long-wave  = downward - upward (?)!!!!!!!!!!!!!!!!
# endif
! --------------------------------------------
#elif defined JRA55
    ! for Pair (Pressure)
      in_data2(1,:,:) = in_data(1,:,:)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V
      in_data2(2,:,:) = in_data(2,:,:)
      in_data2(3,:,:) = in_data(3,:,:)
    ! for Tair
      in_data2(4,:,:) = in_data(4,:,:) - 273.15d0  ! K -> degC
    ! for Qait (Relative humidity)
      in_data2(5,:,:) = in_data(5,:,:)
    ! for cloud (cloud fraction)
      in_data2(6,:,:) = in_data(6,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(7,:,:) = in_data(7,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(8,:,:) = in_data(8,:,:)*0.01  ! percent -> ratio(0 to 1)
      in_data2(9,:,:) = in_data(9,:,:)*0.01  ! percent -> ratio(0 to 1)
    ! for rain (Total precipitation rate)
      in_data2(10,:,:) = in_data(10,:,:)*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1  
# if defined BULK_FLUX
    ! for bulk flux parameters, heating -> positive?, cooling -> negattive?
      in_data2(11,:,:) = -in_data(11,:,:)  ! W m-2 -> W m-2
      in_data2(12,:,:) = -in_data(12,:,:)  ! W m-2 -> W m-2
      in_data2(13,:,:) = in_data(13,:,:) !-in_data(14,:,:) ! solar radiation = downward (?)!!!!!!!!!!!!!!!!
      in_data2(14,:,:) = in_data(15,:,:) ! downward long-wave (?)!!!!!!!!!!!!!!!!
      in_data2(15,:,:) = in_data(16,:,:)-in_data(15,:,:) ! net long-wave  = upward - downward (?)!!!!!!!!!!!!!!!!
# endif
#endif         
!  ---- LOOP3.2 START --------------------------------
      DO iparam=1,N_OutPar
        write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )

!#if defined DSJRA55
!        call interp2D_grid2(Im, Jm, in_data2(iparam,:,:)       &   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          , Nxr, Nyr, out_data(iparam,:,:)     &
!                          , Id_cont, w_cont)
!#else      
        call interp2D_grid (Im, Jm, in_data2(iparam,:,:)       &  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          , Nxr, Nyr, out_data(iparam,:,:)     &
                          , Id_cont, w_cont)
!#endif
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
#if defined JRA55
      time(1) = t+1.5d0/24.0d0
#else
      time(1) = t+0.5d0/24.0d0
#endif
      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d( 'time', trim( FRC_FILE(2) )                  &
            , 1, time, start1D, count1D )

      DO iparam=1,N_OutPar
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
      
END PROGRAM frcATM2ROMS
      
