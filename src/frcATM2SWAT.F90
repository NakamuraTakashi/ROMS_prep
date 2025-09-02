
!!!=== Copyright (c) 2014-2025 Takashi NAKAMURA  =====

#if defined JMA_MSM
# undef BULK_FLUX
# if defined NETCDF_INPUT
#  undef SWRAD
# endif
#elif defined DSJRA55 || defined JRA55
# undef SWRAD
# undef NETCDF_INPUT
# define BULK_FLUX
#elif defined ERA5
# define NETCDF_INPUT
#endif

PROGRAM frcATM2SWAT
  use netcdf
  use eccodes
  use mod_roms_netcdf
  use mod_interpolation
  use mod_calendar
 
  implicit none
  
!-------------------------------------------------------------------------------
  real(8) :: Tlat, Blat, Llon, Rlon
  character(256) :: BATH_FILE
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
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
  integer, parameter :: N_OutPar = 11  ! 10 -> 11
# endif
  integer, parameter :: N_OutPar1 = 9

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
    character(len=*), parameter :: NC_LAT_NAME  = 'lat'
    character(len=*), parameter :: NC_LON_NAME  = 'lon'
    character(len=*), parameter :: NC_TIME_NAME = 'time'
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
  integer, parameter :: N_OutPar1 = 9

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
  integer, parameter :: N_OutPar1 = 9

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
#elif defined ERA5
!--- REA5 parameter setting -----------------
  character(256) :: FRC_prefix
!  integer :: jd_msmnew

  integer, parameter :: N_InPar  = 8
  integer, parameter :: N_OutPar = 8
  integer, parameter :: N_OutPar1 = 5

  character(256) :: NCIN_NAME(N_InPar) = (/   &
     "msl  "  &  ! Mean sea level pressure (Pa)
    ,"u10  "  &  ! 10 metre U wind component (m s-1)
    ,"v10  "  &  ! 10 metre V wind component (m s-1)
    ,"t2m  "  &  ! 2 metre temperature (K)
    ,"d2m  "  &  ! 2 metre dewpoint temperature (K)
    ,"tp   "  &  ! Total precipitation (m h-1)
    ,"ssrd "  &  ! Surface solar radiation downward (J h-1 m-2)
    ,"strd "  &  ! Surface thermal radiation downward (J h-1 m-2)
    /)
  character(len=*), parameter :: NC_LAT_NAME  = 'latitude'
  character(len=*), parameter :: NC_LON_NAME  = 'longitude'
  character(len=*), parameter :: NC_TIME_NAME = 'time'
  integer :: d_ref_days
  real(8) :: sf, off

#endif

#if defined ERA5
  TYPE T_NC
    real(8), pointer :: time_all(:)
    integer :: Nt
    integer :: ItS, ItE
  END TYPE T_NC
  TYPE (T_NC), allocatable :: NC(:)
  real(8), allocatable :: atm_time(:)
  real(8), allocatable :: time2(:)
  integer, allocatable :: iNCt(:)
  integer, allocatable :: idt(:)
  integer :: NCnum
  integer :: iNCs, iNCe
  character(256), allocatable :: ATM_FILE(:)
  real(8) :: d_jdate
  real(8) :: d_jdate_Start, d_jdate_Ref, d_jdate_Ref_atm
  integer :: jdate_Start, jdate_End, jdate_Ref_atm
  integer :: N_days
  integer :: iNC, iNCm
  integer :: Nt

#endif

  character(10) :: GRIB_yyyymmddhh = "2002070100"

  integer :: ifile,idx,iret,igrib
  integer :: istart, iend
  integer :: jdate_Ref, jd_now
  integer :: YYYYMMDD(N_InPar), hhmm(N_InPar)
  real(8), allocatable :: values(:)
  integer :: p1,p2,p3
  character(256) :: p4

  character(9) :: FRC_yyyymmdd = "_20020701" 
!  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: FRC_FILE(2)

  real(8), parameter :: PI = 3.141592653589793d0
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: cosAu(:, :)
  real(8), allocatable :: sinAu(:, :)
  real(8), allocatable :: cosAv(:, :)
  real(8), allocatable :: sinAv(:, :)
  
  real(8), allocatable :: in_data(:,:, :), in_data2(:,:, :)
  real(8), allocatable :: out_data(:,:,:) ! output forcing data
  logical, allocatable :: in_range(:,:)
       
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
  integer :: start1D(1), count1D(1)
  integer :: start3D(3), count3D(3)
  
  character(256) :: IN_FILE(2), IN_FILE2(2)
  integer :: iyear, imonth, iday, iyear2
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
  
  integer :: iparam,ifc,iin,isp
  real(8) :: d_lat, d_lon
  real(8) :: u, v
  real(8) :: dpT, sat_VP, VP
  logical :: file_exists

  integer :: Nout
  integer :: iout
  character(256), allocatable :: OUT_FILE(:,:)
  character(4) :: CNUM
  real(8) :: pcp, tmp, hmd, slr, wnd
  real(8) :: tot_pcp, max_tmp, min_tmp, ave_hmd, tot_slr, ave_wnd
  integer :: nbyr
  integer :: count_subday, count_day 
  integer :: iost
  integer :: itstep

  real(8), allocatable :: lat_all(:), lon_all(:)
  real(8), allocatable :: elev(:,:)
  integer :: N_lat_all, N_lon_all
  real(8) :: plat, plon, pelev
  real(8) :: cff, cff1, cff2,zenith,LatRad,Dangle,Hangle
  real(8) :: dhour

  real(8), parameter :: deg2rad = PI / 180.0d0
  real(8), parameter :: Csolar = 1353.0d0          ! 1360-1380 W/m2

!-------------------------------------------------------------------------------
  namelist/range/Tlat, Blat, Llon, Rlon
  namelist/gebco/BATH_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
#if defined JMA_MSM
  namelist/frc_jmamsm/MSM_dir
  namelist/frc_jmamsm/FRC_prefix
#elif defined DSJRA55
  namelist/frc_dsjra55/DSJRA55_dir
  namelist/frc_dsjra55/FRC_prefix
#elif defined JRA55
  namelist/frc_jra55/JRA55_dir
  namelist/frc_jra55/FRC_prefix
#elif defined ERA5
  namelist/frc_era5_1/NCnum
  namelist/frc_era5_2/ATM_FILE
  namelist/frc_era5_2/FRC_prefix
#endif
  ! Read parameters in namelist file

  read (5, nml=range)
  rewind(5)
  read (5, nml=gebco)
  rewind(5) 
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
#if defined JMA_MSM
  rewind(5)
  read (5, nml=frc_jmamsm)
#elif defined DSJRA55
  rewind(5)
  read (5, nml=frc_dsjra55)
#elif defined JRA55
  rewind(5)
  read (5, nml=frc_jra55)
#elif defined ERA5
  rewind(5)
  read (5, nml=frc_era5_1)
  allocate( ATM_FILE(NCnum) )
  allocate( NC(NCnum) )
  rewind(5)
  read (5, nml=frc_era5_2)
#endif

!---- Modify time-unit description ---------------------------------

!  TIME_ATT(12:15)=YYYY
!  TIME_ATT(17:18)=MM
!  TIME_ATT(20:21)=DD
      
  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  write (DD, "(I2.2)") Sday

!==== Read ATM file coordinates ==================================

#if defined NETCDF_INPUT
!---- NetCDF input file -------------------------------------------
# if defined JMA_MSM
!---- Read JMA-MSM NetCDF file --------------------------------
  IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//DD//".nc"
  IN_FILE(2) = trim(MSM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

# elif defined ERA5
!---- Read ERA5 NetCDF file --------------------------------
  IN_FILE(1) = trim(ATM_FILE(1))
# endif

  write(*,*) "OPEN: ", trim( IN_FILE(1) )
  !Open NetCDF file
  call check( nf90_open(trim( IN_FILE(1) ), nf90_nowrite, ncid) )

  ! Get dimension data
  call get_dimension(ncid, NC_LAT_NAME, Jm)
  call get_dimension(ncid, NC_LON_NAME, Im)
  write(*,*) Im,Jm
  
  ! Allocate variable
  allocate(lat(Jm))
  allocate(lon(Im))
  
  call check( nf90_inq_varid(ncid, NC_LAT_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lat) )
  call check( nf90_inq_varid(ncid, NC_LON_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lon) )
  call check( nf90_close(ncid) )

#else

!---- GRIB2 input file -------------------------------------------
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
!  ---- Get Lat Lon coordinates from Binary file --------------
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

# else

  call codes_get(igrib,'distinctLatitudes',lat)
  call codes_get(igrib,'distinctLongitudes',lon)

# endif 

  call codes_release(igrib)
  call codes_close_file(ifile)
#endif

  allocate(in_range(Im, Jm))
  allocate(in_data(Im, Jm, N_InPar))
  allocate(in_data2(Im, Jm, N_OutPar))

!---- Read GEBCO grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( BATH_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( BATH_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'lat', N_lat_all)
  call get_dimension(ncid, 'lon', N_lon_all)
  ! Allocate variable
  allocate(lat_all(N_lat_all))
  allocate(lon_all(N_lon_all))
  allocate(elev(N_lon_all,N_lat_all))
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat', var_id) )
  call check( nf90_get_var(ncid, var_id, lat_all) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) )
  call check( nf90_get_var(ncid, var_id, lon_all) )
  call check( nf90_inq_varid(ncid, 'elevation', var_id) )
  call check( nf90_get_var(ncid, var_id, elev) )

!==== Pickup indices inside the range. ===============================================
  Nout = 0
  do j=1, Jm
    do i=1,Im
#if defined DSJRA55
      if( Blat<=lat(i,j) .and. lat(i,j)<=Tlat .and.   &
          Llon<=lon(i,j) .and. lon(i,j)<=Rlon ) then
#else
      if( Blat<=lat(j) .and. lat(j)<=Tlat .and.   &
          Llon<=lon(i) .and. lon(i)<=Rlon ) then
#endif
        in_range(i,j) = .true.
        Nout = Nout + 1
      else
        in_range(i,j) = .false.
      endif  
    end do
  end do

!==== Create output files ===================================

  allocate(OUT_FILE(Nout, N_OutPar))

  FRC_yyyymmdd(2:5)=YYYY
  FRC_yyyymmdd(6:7)=MM
  FRC_yyyymmdd(8:9)=DD

  iout = 0
  do j=1, Jm
    do i=1,Im
      if(in_range(i,j)) then

        iout = iout + 1
        write (CNUM, "(I4.4)") iout
        OUT_FILE(iout,1) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'_tmp.txt'
    
        write(*,*) "CREATE: ", trim( OUT_FILE(iout,1) )
        open(10+iout, file=OUT_FILE(iout,1), status='replace')
#if defined DSJRA55
        plat = lat(i,j)
        plon = lon(i,j)
#else
        plat = lat(j)
        plon = lon(i)
#endif
        call LinearInterpolation2D_point(                     &
                N_lon_all, N_lat_all, lon_all, lat_all, elev  & 
              , plon, plat, pelev )

        pelev = max(pelev, 0.0d0)      
        write(10+iout,*) "LAT     LONG     ELEV"
        write(10+iout,*) plat, plon, pelev
        write(10+iout,*) "Year     Month     Day     Hour     pcp     tmp     hmd     slr     wnd"

      endif
    end do
  end do

!==== LOOP set up ==========================================
  itime = 1
  ihours = 0
  iyear = Syear
  imonth = Smonth
  iday = Sday
  call jd(iyear, imonth, iday, Sjdate)
#if defined JMA_MSM
  call jd(2006, 3, 1, jd_msmnew)
#endif

!  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  call jd(2000, 1, 1, jdate_Ref)

#if defined ERA5
!======= Set starting time and Ending time ========================

  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)

  jdate_End = jdate_Start + N_days
  
  write(*,*) jdate_Start, jdate_End, N_days

  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref

  call jd(1900, 1, 1, jdate_Ref_atm)! ERA5 time: hours since 1900-01-01 00:00:00
  d_jdate_Ref_atm = dble(jdate_Ref_atm)
  write(*,*) d_jdate_Ref_atm

! Check time
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    ! Open NetCDF file
    write(*,*) "OPEN: ", trim( ATM_FILE(iNC) )
    call check( nf90_open(trim( ATM_FILE(iNC) ), nf90_nowrite, ncid) )
    call get_dimension(ncid, NC_TIME_NAME, NC(iNC)%Nt)
    write(*,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    allocate( time2(NC(iNC)%Nt) )
    call check( nf90_inq_varid(ncid, NC_TIME_NAME, var_id) )
    call check( nf90_get_var(ncid, var_id, time2) )
    call check( nf90_close(ncid) )

    NC(iNC)%time_all = time2/24.0d0 + d_jdate_Ref_atm ! hour -> day

    deallocate(time2)
  end do

  write(*,*) NC(:)%Nt

  do iNC=1, NCnum-1
    do i=2,NC(iNC)%Nt
      if( NC(iNC+1)%time_all(1) <= NC(iNC)%time_all(i) ) then
        NC(iNC)%Nt = i-1
        exit
      end if
    end do
  end do
  
  write(*,*) NC(:)%Nt

  write(*,*) "******************************************************************"

  NC(:)%ItE = -1
  NC(:)%ItS = -1
  
  do iNC=NCnum,1,-1
    do i=NC(iNC)%Nt,1,-1
      d_jdate = NC(iNC)%time_all(i)
      if(d_jdate < dble(jdate_End)) then
        write(*,*) '*** FOUND: Ending point @ ATM_FILE',iNC
        NC(iNC)%ItE=i
        exit
      endif
    end do
  end do
  write(*,*) NC(:)%ItE 
  
  do iNC=1,NCnum
    do i=NC(iNC)%ItE,1,-1
      d_jdate = NC(iNC)%time_all(i)
      if(d_jdate < dble(jdate_Start)) then
!        write(*,*) '*** FOUND: Starting point @ ATM_FILE',iNC
        exit
      endif
      NC(iNC)%ItS=i
    end do
  end do
  write(*,*) NC(:)%ItS 

  Nt = 0
  iNCs = NCnum
  iNCe = 1

  do iNC=1,NCnum
    if(NC(iNC)%ItS==-1) then
      cycle
    end if
    Nt = Nt + NC(iNC)%ItE - NC(iNC)%ItS + 1
    iNCs = min(iNCs,iNC)
    iNCe = max(iNCe,iNC)
  enddo

  allocate(atm_time(Nt))
  allocate(idt(Nt))
  allocate(iNCt(Nt))

  i=1
  do iNC=iNCs,iNCe
    do k=0,NC(iNC)%ItE-NC(iNC)%ItS
      iNCt(i+k) = iNC
      idt(i+k)  = NC(iNC)%ItS + k
    enddo
    j=i+NC(iNC)%ItE-NC(iNC)%ItS
    atm_time(i:j) = NC(iNC)%time_all( NC(iNC)%ItS : NC(iNC)%ItE ) &
                    - d_jdate_Ref
    i=j+1
  enddo

  write(*,*) "*************************************"
  
  iNC = 0
#endif

!===== LOOP1 START ================================================
  DO
#if defined ERA5
    ! Check end date
    if(itime>Nt) then
      write(*,*) "Completed!!!"
      EXIT
    endif
#else
    ihour = mod(ihours,24)
    ijdate = Sjdate + int(ihours/24)
    call cdate( ijdate, iyear, imonth, iday )
    ! Check end date
    if(iyear==Eyear .and. imonth==Emonth .and. iday==Eday) then
      write(*,*) "Completed!!!"
      EXIT
    endif

    write (YYYY, "(I4.4)") iyear
    write (MM, "(I2.2)") imonth
    write (DD, "(I2.2)") iday ! 1+int((itime-1)*1/24)
    write (hh, "(I2.2)") ihour ! mod((itime-1)*1,24)
#endif

#if defined NETCDF_INPUT
# if defined JMA_MSM
    ihours = ihours + 24  !!! Files exist daily interval

    IN_FILE2(1) = IN_FILE(1)  
    IN_FILE2(2) = IN_FILE(2)  
    IN_FILE(1) = trim(MSM_dir)//YYYY//"/"//MM//DD//".nc"
    IN_FILE(2) = trim(MSM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

# elif defined ERA5
 
    iNCm = iNC
    iNC  = iNCt(itime)
    if(iNCm < iNC) then
      IN_FILE(1) = trim(ATM_FILE(iNC))
    endif

# endif
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
#elif defined ERA5
    DO ifc=NC(iNC)%ItS,NC(iNC)%ItE
#else
    DO ifc=1,1
#endif
!  ---- LOOP3.1 START --------------------------------
      DO iparam=1,N_InPar
#if defined JMA_MSM
        if(iparam==N_OutPar1 .and. ijdate<jd_msmnew) cycle
#endif
#if defined NETCDF_INPUT
# if defined JMA_MSM
        if(iparam>N_OutPar1) then !!! for rain (rain fall rate)
          iin=2
        else
          iin=1
        end if
# elif defined ERA5
        iin=1
# endif

        write(*,*) "OPEN: ", trim( IN_FILE(iin) )
        !Open NetCDF file
        call check( nf90_open(trim( IN_FILE(iin) ), nf90_nowrite, ncid) )

        start3D = (/ 1,  1,  ifc /)
        count3D = (/ Im, Jm, 1   /)
        ! Get variable id
        write(*,*) "READ: ", trim( NCIN_NAME(iparam) )
        call check( nf90_inq_varid(ncid, trim( NCIN_NAME(iparam) ), var_id) ) 
# if defined JMA_MSM
        if(ijdate>=jd_msmnew .or.iparam<=5 .or.iparam==10 ) then
          call check( nf90_get_var(ncid, var_id, in_data(:,:,iparam), start=start3D, count=count3D) )
          call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
          call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
          in_data(:,:,iparam)=in_data(:,:,iparam)*sf+off
        else
          call check( nf90_get_var(ncid, var_id, in_data(:,:,iparam), start=start3D, count=count3D) )
        endif
# elif defined ERA5
        call check( nf90_get_var(ncid, var_id, in_data(:,:,iparam), start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        in_data(:,:,iparam)=in_data(:,:,iparam)*sf+off
# endif
        call check( nf90_close(ncid) )     
   
#else
!      ---- Seek message --------------------------------
# if defined JMA_MSM || defined ERA5
        iin=1
# else
        if(iparam>=10) then !!! for rain (rain fall rate)
          iin=2
        else
          iin=1
        end if
# endif
        write(*,*) "OPEN: ", trim( IN_FILE(iin) )
        call codes_open_file(ifile, trim( IN_FILE(iin) ),'r', iret)
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
#elif defined ERA5
          call codes_get(igrib,'dataDate' , p1)
          call codes_get(igrib,'dataTime' , p2)
          call codes_get(igrib,'endStep'  , p3)
          call codes_get(igrib,'shortName', p4)
          if (p1==GRIB_STEP(ifc) .and.             &
              trim(p4)==trim(GRIB_NAME(iparam))  ) exit
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
          in_data(:,i,iparam) = values(istart:iend)
        end do
#endif
      END DO

#if defined NETCDF_INPUT
    ! Set date & time
      call check( nf90_open(trim( IN_FILE(1) ), nf90_nowrite, ncid) )
      start1D = (/ ifc /)
      count1D = (/ 1 /)
      call check( nf90_inq_varid(ncid, NC_TIME_NAME, var_id) )  !!!  not Japan time (00:00:00+09:00)
      call check( nf90_get_var(ncid, var_id, time, start=start1D, count=count1D) )
      call check( nf90_close(ncid) )
# if defined JMA_MSM
    ! JMA_MSM time: hours since 2000-01-01 00:00:00
      call ndays(imonth, iday, iyear, 1, 1, 2000, d_ref_days)
      t = time(1)/24.0d0 + dble(d_ref_days)

# elif defined ERA5
    ! ERA5 time:    hours since 1900-01-01 00:00:00
      t = atm_time(itime)
# endif
      jd_now=int(t)+jdate_Ref
      call cdate(jd_now, iyear, imonth, iday)
      call ndays(imonth, iday, iyear, 1, 1, iyear, idays)
!      dhour = (t-dble(floor(t)))*24.0d0
      ihour = nint(mod(t,1.0d0)*24.0d0) 

#else
    ! Set date & time
      iyear  = YYYYMMDD(1)/10000
      imonth = (YYYYMMDD(1)-iyear*10000)/100
      iday   = YYYYMMDD(1)-iyear*10000-imonth*100
      ihour  = hhmm(1)/100
      imin   = hhmm(1)-100*ihour

    !  ihour  = ihour-1 ! since time for precipitation is set +1 hour

      call ndays(imonth, iday, iyear, 1, 1, 2000, idays)
      
      t = dble(idays)+dble(ihour)/24.0d0+dble(imin)/1440.0d0
#endif
!  ---- LOOP3.1 END --------------------------------
#if defined JMA_MSM
    ! for Pair (Pressure)
      if(ijdate<jd_msmnew) then
        in_data2(:,:,1) = in_data(:,:,1)  ! hPa -> millibar (= hPa) 
      else   
        in_data2(:,:,1) = in_data(:,:,1)*0.01  ! Pa -> millibar (= hPa) 
      endif
    ! for U V
      in_data2(:,:,2) = in_data(:,:,2)
      in_data2(:,:,3) = in_data(:,:,3)
    ! for Tair
      if(ijdate<jd_msmnew) then
        in_data2(:,:,4) = in_data(:,:,4)   ! degC -> degC 
      else   
        in_data2(:,:,4) = in_data(:,:,4) - 273.15d0  ! K -> degC
      endif
    ! for Qait (Relative humidity)
      in_data2(:,:,5) = in_data(:,:,5)
    ! for cloud (cloud fraction) 
      if(ijdate<jd_msmnew) then
      ! no total cloud cover data
        in_data2(:,:,6) = in_data(:,:,6)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(:,:,7) = in_data(:,:,7)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(:,:,8) = in_data(:,:,8)/9.0d0  ! 0-9 level -> ratio(0 to 1)
        in_data2(:,:,9) = 1.0d0 - (1.0-in_data2(:,:,6))*(1.0-in_data2(:,:,7)) &
                                 *(1.0-in_data2(:,:,8)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else   
        in_data2(:,:,6) = in_data(:,:,6)*0.01  ! percent -> ratio(0 to 1)
        in_data2(:,:,7) = in_data(:,:,7)*0.01  ! percent -> ratio(0 to 1)
        in_data2(:,:,8) = in_data(:,:,8)*0.01  ! percent -> ratio(0 to 1)
        in_data2(:,:,9) = in_data(:,:,9)*0.01  ! percent -> ratio(0 to 1)
      endif
    ! for rain (Total precipitation rate)
      in_data2(:,:,10) = in_data(:,:,10)/3600.0d0  ! kg m-2 h-1 -> kg m-2 s-1 (= mm s-1) 
# if defined SWRAD
    ! for short-wave radiation
      in_data2(:,:,11) = in_data(:,:,11)  ! W m-2 -> W m-2
#else

!      jd_now=int(t)+jdate_Ref
!      call cdate(jd_now, iyear, imonth, iday)
!      call ndays(imonth, iday, iyear, 1, 1, iyear, idays)
!      dhour = (t-dble(floor(t)))*24.0d0

      dhour = dble(ihour)

      Dangle=23.44d0*COS((172.0d0-dble(idays))*2.0d0*pi/365.25d0)
      Dangle=Dangle*deg2rad
      Hangle=(12.0d0-dhour)*pi/12.0d0

      do j=1, Jm
        do i=1,Im
!  Local daylight is a function of the declination (Dangle) and hour
!  angle adjusted for the local meridian (Hangle-lonr(i,j)/15.0).
!  The 15.0 factor is because the sun moves 15 degrees every hour.
!
          LatRad=lat(j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
!
!  Estimate variation in optical thickness of the atmosphere over
!  the course of a day under cloudless skies (Zillman, 1972). To
!  obtain net incoming shortwave radiation multiply by (1.0-0.6*c**3),
!  where c is the fractional cloud cover.
!
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606).
!
          in_data2(i,j,11) = 0.0d0 ! for short-wave radiation

          zenith=cff1+cff2*COS(Hangle-lon(i)*deg2rad)
          IF (zenith.gt.0.0d0) THEN
            sat_VP=6.1078d0*exp(5423.0d0*(1.0d0/273.15d0-1.0d0/(in_data2(i,j,4)+273.15d0)))   ! saturation vapor pressure (hPa=mbar)
            VP    =sat_VP*in_data2(i,j,5)/100.0d0 ! water vapor pressure (hPa=mbar)  ! in_data2(i,j,5) -> humidity (%)
            in_data2(i,j,11) =Csolar*zenith*zenith*                           & ! for short-wave radiation
                (1.0d0-1.5d0*in_data2(i,j,9)**2+0.98d0*in_data2(i,j,9)**3)/   & ! (1) for JMAMSM data   ! in_data2(i,j,9) -> cloud fraction (0-1)
!               (1.0d0-0.6d0*cloud(i,j)**3)/                                  & ! (2) default
                ((zenith+2.7d0)*VP*1.0d-3+1.085d0*zenith+0.1d0)
          END IF
        end do
      end do
# endif 
! --------------------------------------------
#elif defined DSJRA55
    ! for Pair (Pressure)
      in_data2(:,:,1) = in_data(:,:,1)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V: change DSJRA55 Lambert conformal to regular Lat Lon coordinat vectors 
      do j=1, Jm
        do i=1,Im
          in_data2(i,j,2) = in_data(i,j,2)*cosAx(i,j)-in_data(i,j,3)*sinAy(i,j)
          in_data2(i,j,3) = in_data(i,j,2)*sinAx(i,j)+in_data(i,j,3)*cosAy(i,j)
        end do
      end do
    ! for Tair
      in_data2(:,:,4) = in_data(:,:,4) - 273.15d0  ! K -> degC
    ! for Qair: convert Dewpoint depression (K) to Relative humidity (%)
      do j=1, Jm
        do i=1,Im
          ! Dewpoint (K)
          dpT = in_data(i,j,4)-in_data(i,j,5) ! (K)
          ! Saturation vapor pressure (hPa)
          sat_VP=6.1078d0*exp(5423.0d0*(1.0d0/273.15d0-1.0d0/in_data(i,j,4)))
          ! Vapor pressure (hPa)
          VP    =6.1078d0*exp(5423.0d0*(1.0d0/273.15d0-1.0d0/dpT))
          ! Relative humidity (%)
          in_data2(i,j,5) = VP/sat_VP*100.0d0 ! (%)
        end do
      end do
    ! for cloud (cloud fraction)
      in_data2(:,:,6) = in_data(:,:,6)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,7) = in_data(:,:,7)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,8) = in_data(:,:,8)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,9) = in_data(:,:,9)*0.01  ! percent -> ratio(0 to 1)
    ! for rain (Total precipitation rate)
      in_data2(:,:,10) = in_data(:,:,10)  ! kg m-2 s-1 -> kg m-2 s-1  (= mm s-1) 
# if defined BULK_FLUX
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,11) = -in_data(:,:,11)  ! W m-2 -> W m-2
      in_data2(:,:,12) = -in_data(:,:,12)  ! W m-2 -> W m-2
      in_data2(:,:,13) =  in_data(:,:,13) !-in_data(14,:,:) ! solar radiation = downward
      in_data2(:,:,14) =  in_data(:,:,14) ! downward long-wave
      in_data2(:,:,15) =  in_data(:,:,15)-in_data(:,:,16) ! net long-wave  = downward - upward
# endif
! --------------------------------------------
#elif defined JRA55
    ! for Pair (Pressure)
      in_data2(:,:,1) = in_data(:,:,1)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V
      in_data2(:,:,2) = in_data(:,:,2)
      in_data2(:,:,3) = in_data(:,:,3)
    ! for Tair
      in_data2(:,:,4) = in_data(:,:,4) - 273.15d0  ! K -> degC
    ! for Qait (Relative humidity)
      in_data2(:,:,5) = in_data(:,:,5)
    ! for cloud (cloud fraction)
      in_data2(:,:,6) = in_data(:,:,6)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,7) = in_data(:,:,7)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,8) = in_data(:,:,8)*0.01  ! percent -> ratio(0 to 1)
      in_data2(:,:,9) = in_data(:,:,9)*0.01  ! percent -> ratio(0 to 1)
    ! for rain (Total precipitation rate)
      in_data2(:,:,10) = in_data(:,:,10)*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1 (= mm s-1) 
# if defined BULK_FLUX
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,11) = -in_data(:,:,11)  ! W m-2 -> W m-2
      in_data2(:,:,12) = -in_data(:,:,12)  ! W m-2 -> W m-2
      in_data2(:,:,13) =  in_data(:,:,13) !-in_data(14,:,:) ! solar radiation = downward
      in_data2(:,:,14) =  in_data(:,:,14) ! downward long-wave
      in_data2(:,:,15) =  in_data(:,:,15)-in_data(:,:,16) ! net long-wave  = downward - upward
# endif
! --------------------------------------------
#elif defined ERA5
    ! for Pair (Pressure)
      in_data2(:,:,1) = in_data(:,:,1)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V
      in_data2(:,:,2) = in_data(:,:,2)
      in_data2(:,:,3) = in_data(:,:,3)
    ! for Tair
      in_data2(:,:,4) = in_data(:,:,4) - 273.15d0  ! K -> degC
    ! for Qair: convert Dewpoint temperature (K) to Relative humidity (%)
      do j=1, Jm
        do i=1,Im
          ! Saturation vapor pressure (hPa)
          sat_VP=6.1078d0*exp(5423.0d0*(1.0d0/273.15d0-1.0d0/in_data(i,j,4)))
          ! Vapor pressure (hPa)
          VP    =6.1078d0*exp(5423.0d0*(1.0d0/273.15d0-1.0d0/in_data(i,j,5)))
          ! Relative humidity (%)
          in_data2(i,j,5) = VP/sat_VP*100.0d0 ! (%)
        end do
      end do
    ! for rain (Total precipitation rate)
      in_data2(:,:,6) = in_data(:,:,6)*1000.0d0/3600.0d0  ! m h-1 -> kg m-2 s-1 (= mm s-1)   
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,7) = in_data(:,:,7)/3600.0d0  ! J h-1 m-2 -> W m-2
      in_data2(:,:,8) = in_data(:,:,8)/3600.0d0  ! J h-1 m-2 -> W m-2
#endif
    
!  ---- Write data to temporary files --------------------------------  
      
      iout = 0
      do j=1, Jm
        do i=1,Im
          if(in_range(i,j)) then
            iout = iout + 1
#if defined JMA_MSM
            pcp = max( in_data2(i,j,10)*3600.0d0, 0.0d0 ) ! kg m-2 s-1 (= mm s-1)  -> mm h-1 
            tmp = in_data2(i,j,4)  ! degC
            hmd = in_data2(i,j,5)  ! %
            wnd = sqrt( in_data2(i,j,2)*in_data(i,j,2)     &
                       +in_data2(i,j,3)*in_data2(i,j,3) )  ! m s-1
            slr = in_data2(i,j,11)*3.6d-3  ! W m-2 -> MJ m-2 h-1

#elif defined DSJRA55 || defined JRA55
            pcp = max( in_data2(i,j,10)*3600.0d0, 0.0d0 ) ! kg m-2 s-1 (= mm s-1)  -> mm h-1 
            tmp = in_data2(i,j,4)  ! degC
            hmd = in_data2(i,j,5)  ! %
            wnd = sqrt( in_data2(i,j,2)*in_data(i,j,2)     &
                       +in_data2(i,j,3)*in_data2(i,j,3) )  ! m s-1
# if defined BULK_FLUX
            slr = in_data2(i,j,13)*3.6d-3  ! W m-2 -> MJ m-2 h-1 ! W m-2
# else
            slr = To be developed  *3.6d-3  ! W m-2 -> MJ m-2 h-1! W m-2
# endif
#elif defined ERA5
            pcp = max( in_data2(i,j,6)*3600.0d0, 0.0d0 )  ! kg m-2 s-1 (= mm s-1) -> mm h-1 
            tmp = in_data2(i,j,4)  ! degC
            hmd = in_data2(i,j,5)  ! %
            wnd = sqrt( in_data2(i,j,2)*in_data(i,j,2)     &
                       +in_data2(i,j,3)*in_data2(i,j,3) )  ! m s-1
            slr = in_data2(i,j,7)*3.6d-3  ! W m-2 -> MJ m-2 h-1
#endif
!            write(10+iout,*) t, pcp, tmp, hmd, slr, wnd
            write(10+iout,*) iyear, imonth, iday, ihour, pcp, tmp, hmd, slr, wnd

          endif
        end do
      end do
!  ---- LOOP3.3 END --------------------------------         
      itime = itime + 1
    END DO
! ---- LOOP2 END --------------------------------
  END DO
!---- LOOP1 END --------------------------------

  DO iout=1,Nout
    close(10+iout)
  END DO

!==== CREATE SWAT+ weather input files ===========================

  open ( 9, file='pcp_sd.cli', status='replace')
  write( 9,'("Hourly precipitation file names")') 
  write( 9,'("FILENAME")') 

  open (10, file='pcp.cli', status='replace')
  write(10,'("Daily precipitation file names")') 
  write(10,'("FILENAME")') 

  open (11, file='tmp.cli', status='replace')
  write(11,'("Daily temperature file names")') 
  write(11,'("FILENAME")') 

  open (12, file='hmd.cli', status='replace')
  write(12,'("Daily humidity file names")') 
  write(12,'("FILENAME")') 

  open (13, file='slr.cli', status='replace')
  write(13,'("Daily solar radiation file names")') 
  write(13,'("FILENAME")') 

  open (14, file='wnd.cli', status='replace')
  write(14,'("Daily wind speed file names")') 
  write(14,'("FILENAME")') 

  nbyr = iyear - Syear + 1

  DO iout=1,Nout
    open(15, file=OUT_FILE(iout,1), status='old')
!   Read header information
    read (15,*)
    read (15,*) plat, plon, pelev
    read (15,*)

!   Create SWAT+ input files
    write (CNUM, "(I4.4)") iout
    OUT_FILE(iout,2) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'_sd.pcp'
    OUT_FILE(iout,3) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'.pcp'
    OUT_FILE(iout,4) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'.tmp'
    OUT_FILE(iout,5) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'.hmd'
    OUT_FILE(iout,6) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'.slr'
    OUT_FILE(iout,7) = trim(FRC_prefix)//FRC_yyyymmdd//'_'//CNUM//'.wnd'

    open (16, file=OUT_FILE(iout,2), status='replace')
    write(16,'("Hourly precipitation (mm h-1)")') 
    write(16,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(16,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 24, plat, plon, pelev

    open (17, file=OUT_FILE(iout,3), status='replace')
    write(17,'("Daily precipitation (mm day-1)")') 
    write(17,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(17,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 0, plat, plon, pelev

    open (18, file=OUT_FILE(iout,4), status='replace')
    write(18,'("Daily maximum & minimum temperature (degC)")') 
    write(18,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(18,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 0, plat, plon, pelev

    open (19, file=OUT_FILE(iout,5), status='replace')
    write(19,'("Daily mean humidity (%)")') 
    write(19,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(19,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 0, plat, plon, pelev

    open (20, file=OUT_FILE(iout,6), status='replace')
    write(20,'("Daily solar radiation (MJ m-2)")') 
    write(20,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(20,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 0, plat, plon, pelev

    open (21, file=OUT_FILE(iout,7), status='replace')
    write(21,'("Daily mean wind speed (m s-1)")') 
    write(21,'("NBYR  TSTEP        LAT       LONG     ELEV")') 
    write(21,'(i4,i7,f11.5,f11.5,f9.2)') nbyr, 0, plat, plon, pelev

! -- loop start ----
    count_day = 0
    iyear2 = Syear

    OUTER: DO
      count_day = count_day + 1

      max_tmp = -100.0d0
      min_tmp = 100.0d0
      tot_pcp = 0.0d0
      ave_hmd = 0.0d0
      tot_slr = 0.0d0
      ave_wnd = 0.0d0

#if defined JRA55
      itstep=8   ! 3-hourly
#else
      itstep=24  ! 1-hourly 
#endif
      count_subday = 0

      DO i=1,itstep
!        read (15,*,iostat=iost) t, pcp, tmp, hmd, slr, wnd
        read (15,*,iostat=iost) iyear, imonth, iday, ihour, pcp, tmp, hmd, slr, wnd
        if(iost/=0) exit OUTER
        if(i==1) then
!          call cdate( nint(t)+jdate_Ref, iyear, imonth, iday )
          if(iyear/=iyear2) then
            iyear2 = iyear
            count_subday = 1
            count_day = 1
          endif
        endif
        ! write *_sd.pcp file 
        write(16,'(i4,i5,i5,f9.2)') iyear, count_day, count_subday, pcp

        max_tmp = max( max_tmp, tmp )
        min_tmp = min( min_tmp, tmp )
        tot_pcp = tot_pcp + pcp
        ave_hmd = ave_hmd + hmd
        tot_slr = tot_slr + slr
        ave_wnd = ave_wnd + wnd

        count_subday = count_subday + 1
  
      END DO

      ave_hmd = ave_hmd/dble(itstep)
      ave_wnd = ave_wnd/dble(itstep)
      ! write *.pcp file 
      write(17,'(i4,i5,f9.2)') iyear, count_day, tot_pcp
      ! write *.tmp file 
      write(18,'(i4,i5,f9.2,f9.2)') iyear, count_day, max_tmp, min_tmp
      ! write *.hmd file 
      write(19,'(i4,i5,f9.2)') iyear, count_day, ave_hmd
      ! write *.slr file 
      write(20,'(i4,i5,f9.2)') iyear, count_day, tot_slr
      ! write *.wnd file 
      write(21,'(i4,i5,f9.2)') iyear, count_day, ave_wnd

    END DO OUTER

    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)

    write( 9,'(a)') trim( OUT_FILE(iout,2) )
    write(10,'(a)') trim( OUT_FILE(iout,3) )
    write(11,'(a)') trim( OUT_FILE(iout,4) )
    write(12,'(a)') trim( OUT_FILE(iout,5) )
    write(13,'(a)') trim( OUT_FILE(iout,6) )
    write(14,'(a)') trim( OUT_FILE(iout,7) )
  
  END DO

  close( 9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
 
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM frcATM2SWAT
      
