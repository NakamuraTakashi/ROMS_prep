
!!!=== Copyright (c) 2014-2026 Takashi NAKAMURA  =====

#if defined JMA_MSM
# undef BULK_FLUX
# if defined NETCDF_INPUT
#  undef SWRAD
# endif
# if defined JMA_MSM_CLOUD_ONLY
#  undef SWRAD
# endif

#elif defined DSJRA55 || defined JRA55
# undef NETCDF_INPUT
# undef SWRAD
# undef JMA_MSM_CLOUD_ONLY

#elif defined JMA_LSM
# undef NETCDF_INPUT
# undef BULK_FLUX
# undef SWRAD
# undef JMA_MSM_CLOUD_ONLY

#elif defined ERA5
# undef JMA_MSM_CLOUD_ONLY

#elif defined FORP_ATM
# define NETCDF_INPUT

#endif

! Every donor uses the shared mod_infile time list + the single "do itime=1,Nt"
! engine, so the input-file bookkeeping below is unconditional.

PROGRAM frcATM2ROMS
  use netcdf
  use eccodes
  use mod_roms_netcdf
  use mod_interpolation
  use mod_calendar
  use mod_infile   ! shared T_INFILE / infile_check_time (ERA5 + FORP_ATM)

  implicit none
  
!-------------------------------------------------------------------------------
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  character(256) :: GRID_FILE
  integer :: Ryear, Rmonth, Rday
  character(256) :: ATM_dir
  character(256) :: FRC_prefix
!----------------------------------------------------------------------
#if defined JMA_MSM
!--- JMA_MSM parameter setting -----------------
  integer :: jd_msmnew

# if defined SWRAD
  integer, parameter :: N_InPar  = 11
  integer, parameter :: N_OutPar = 11
# else
  integer, parameter :: N_InPar  = 10
  integer, parameter :: N_OutPar = 10
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

#elif defined JMA_LSM
!--- JMA_MSM parameter setting -----------------
  integer, parameter :: N_InPar  = 5
  integer, parameter :: N_OutPar = 5
  integer, parameter :: N_OutPar1 = 5

  character(len=*), parameter :: GRIB_prefix  = "LANAL_"
  character(len=*), parameter :: GRIB_suffix  = ".grb2"

  character(256) :: GRIB_NAME(N_InPar) = (/   &
     "10u  "     &
    ,"10v  "     &
    ,"t    "     &
    ,"r    "     &
    ,"prmsl"     &
    /)
  integer :: GRIB_LEVEL(5) = (/ 10,10,2,2,0 /)


#elif defined JRA55
!--- JRA55 parameter setting -----------------
! (ATM_dir / FRC_prefix are already declared globally near the top of the program;
!  the earlier local re-declaration here was a compile error.)

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
  integer :: jd_msmnew

  integer, parameter :: N_InPar  = 8
  integer, parameter :: N_OutPar = 8
  integer, parameter :: N_OutPar1 = 5

# if defined NETCDF_INPUT
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
  real(8) :: sf, off

# else
  character(256) :: GRIB_NAME(N_InPar) = (/   &
     "msl "  &  ! Mean sea level pressure (Pa)
    ,"10u "  &  ! 10 metre U wind component (m s-1)
    ,"10v "  &  ! 10 metre V wind component (m s-1)
    ,"2t  "  &  ! 2 metre temperature (K)
    ,"2d  "  &  ! 2 metre dewpoint temperature (K)
    ,"tp  "  &  ! Total precipitation (m h-1)
    ,"ssrd"  &  ! Surface solar radiation downward (J h-1 m-2)
    ,"strd"  &  ! Surface thermal radiation downward (J h-1 m-2)
    /)
# endif

#elif defined FORP_ATM
!--- FORP parameter setting -----------------
  integer, parameter :: N_InPar  = 8
  integer, parameter :: N_OutPar = 8
  integer, parameter :: N_OutPar1 = 5

  character(256) :: NCIN_NAME(N_InPar) = (/   &
     "pso   "  &  ! Sea Water Pressure at Sea Water Surface (hPa)
    ,"uas   "  &  ! Eastward Near-Surface Wind (cm s-1)
    ,"vas   "  &  ! Northward Near-Surface Wind (cm s-1)
    ,"tas   "  &  ! Near-Surface Air Temperature (degC)
    ,"evs   "  &  ! Water Evaporation Flux Where Ice Free Ocean Over Sea (kg m-2 s-1)
    ,"pr    "  &  ! Precipitation (kg m-2 s-1)
    ,"rsntds"  &  ! Net Downward Shortwave Radiation at Sea Water Surface (W m-2)
    ,"rlntds"  &  ! Surface Net Downward Longwave Radiation (W m-2)
    /)
  character(len=*), parameter :: NC_LAT_NAME  = 'lat'
  character(len=*), parameter :: NC_LON_NAME  = 'lon'
  character(len=*), parameter :: NC_TIME_NAME = 'time'
  real(8) :: sf, off

#endif

  real(8), allocatable :: atm_time(:)
  real(8), allocatable :: time2(:)
  integer, allocatable :: iNCt(:)
  integer, allocatable :: idt(:)
  integer :: NCnum
  character(256), allocatable :: ATM_FILE(:)
  real(8) :: d_jdate
  real(8) :: d_jdate_Ref, d_jdate_Ref_atm
  integer :: jdate_Start, jdate_End, jdate_Ref, jdate_Ref_atm
  integer :: N_days
  integer :: iNC, iNCm
  integer :: Nt
  real(8), allocatable :: amask(:, :), amaskv(:, :)

  character(10) :: GRIB_yyyymmddhh = "2002070100"

  integer :: ifile(2),iret(2),igrib(2)
  integer :: istart, iend
  real(8), allocatable :: values(:)
  integer :: p1,p2,p3
  character(256) :: p4
  integer :: param_count
  integer :: ips, ipe

  character(9) :: FRC_yyyymmdd = "_20020701" 
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: FRC_FILE(2)

  character(256) :: NC_NAME(N_OutPar) = (/    &
     "Pair      "                             &
    ,"Uwind     "                             &
    ,"Vwind     "                             &
    ,"Tair      "                             &
    ,"Qair      "                             &
#if defined JMA_MSM || defined DSJRA55 || defined JRA55
    ,"lcloud    "                             &
    ,"mcloud    "                             &
    ,"hcloud    "                             &
    ,"cloud     "                             &
    ,"rain      "                             &
# if defined SWRAD
    ,"swrad     "                             &
# elif defined BULK_FLUX
    ,"latent    "                             &
    ,"sensible  "                             &
    ,"swrad     "                             &
    ,"lwrad_down"                             &
    ,"lwrad     "                             &
# endif
#elif defined ERA5
    ,"rain      "                             &
    ,"swrad     "                             &
    ,"lwrad_down"                             &
#elif defined FORP_ATM
    ,"rain      "                             &
    ,"swrad     "                             &
    ,"lwrad     "                             &
#endif
    /)
  character(256) :: NC_LNAME(N_OutPar) = (/   &
     "surface air pressure               "    &
    ,"surface u-wind component           "    &
    ,"surface v-wind component           "    &
    ,"surface air temperature            "    &
    ,"surface air relative humidity      "    &
#if defined JMA_MSM || defined DSJRA55 || defined JRA55
    ,"low cloud fraction                 "    &
    ,"medium cloud fraction              "    &
    ,"high cloud fraction                "    &
    ,"cloud fraction                     "    &
    ,"rain fall rate                     "    &
# if defined SWRAD
    ,"solar shortwave radiation flux     "    &
# elif defined BULK_FLUX
    ,"net latent heat flux               "    &
    ,"net sensible heat flux             "    &
    ,"solar shortwave radiation flux     "    &
    ,"downwelling longwave radiation flux"    &
    ,"net longwave radiation flux        "    &
# endif
#elif defined ERA5
    ,"rain fall rate                     "    &
    ,"solar shortwave radiation flux     "    &
    ,"downwelling longwave radiation flux"    &
#elif defined FORP_ATM
    ,"rain fall rate                     "    &
    ,"solar shortwave radiation flux     "    &
    ,"net longwave radiation flux        "    &
#endif
    /)
  character(256) :: NC_UNIT(N_OutPar) = (/    &
     "millibar                 "              &
    ,"meter second-1           "              &
    ,"meter second-1           "              &
    ,"Celsius                  "              &
    ,"percentage               "              &
#if defined JMA_MSM || defined DSJRA55 || defined JRA55
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"0 to 1                   "              &
    ,"kilogram meter-2 second-1"              &
# if defined SWRAD
    ,"watt meter-2             "              &
# elif defined BULK_FLUX
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
# endif
#elif defined ERA5
    ,"kilogram meter-2 second-1"              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
#elif defined FORP_ATM
    ,"kilogram meter-2 second-1"              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
#endif
    /)

  real(8), parameter :: PI = 3.141592653589793d0
  real(8), parameter :: rad2deg = 180.0d0/PI
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: cosAu(:, :)
  real(8), allocatable :: sinAu(:, :)
  real(8), allocatable :: cosAv(:, :)
  real(8), allocatable :: sinAv(:, :)
  
  real(8), allocatable :: in_data(:,:, :), in_data2(:,:, :)
  real(8), allocatable :: out_data(:,:,:) ! output forcing data
       
  integer :: Im, Jm
  integer :: Irdg_min, Irdg_max, Jrdg_min, Jrdg_max  ! Minimum and maximum ID of donor grid index.

#if defined DSJRA55 || defined JMA_LSM
  real(8), allocatable :: lat(:, :), lon(:, :)
  real(8), allocatable :: cosAx(:, :), sinAx(:, :)
  real(8), allocatable :: cosAy(:, :), sinAy(:, :)
#else
  real(8), allocatable :: lat(:), lon(:)
#endif
  real(8), allocatable :: alat(:, :), alon(:, :)
  real(8), allocatable :: alatv(:, :), alonv(:, :)

  real(8) :: time(1)
  integer, allocatable :: ID_cont(:,:)
  real(8), allocatable :: w_cont(:,:)
  integer, allocatable :: ID_contv(:,:)
  real(8), allocatable :: w_contv(:,:)
  integer :: start1D(1), count1D(1)
  integer :: start3D(3), count3D(3)
  
  character(256), allocatable :: IN_FILE(:)
  integer :: iyear, imonth, iday
  integer :: i,j,k
  integer :: ihours,ijdate
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
  
  integer :: iparam,ifc,iin
  real(8) :: d_lat, d_lon
  real(8) :: s_lat, s_lon
  real(8) :: u, v
  real(8) :: dpT, sat_VP, VP
  logical :: file_exists
  real(8) :: xp,yp

!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
#if defined ERA5
  namelist/frc_era5_1/NCnum
  namelist/frc_era5_2/ATM_FILE
  namelist/frc_era5_2/FRC_prefix
#else
  namelist/frc_atm/ATM_dir
  namelist/frc_atm/FRC_prefix
#endif
  ! Read parameters in namelist file

  read (5, nml=grd)
  rewind(5)
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
  rewind(5)
  read (5, nml=refdate)
  rewind(5)
#if defined ERA5
  read (5, nml=frc_era5_1)
  allocate( ATM_FILE(NCnum) )
  rewind(5)
  read (5, nml=frc_era5_2)
#else
  read (5, nml=frc_atm)
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
  allocate(out_data(0:L, 0:M, N_OutPar))
  allocate(ID_cont(4, Nxr*Nyr))
  allocate(w_cont(4, Nxr*Nyr))
#if defined FORP_ATM
  allocate(ID_contv(4, Nxr*Nyr))
  allocate(w_contv(4, Nxr*Nyr))
#endif

  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr) )
  
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
  do j=0,M
    do i=0,L-1
      d_lat=latr(i+1,j)-latr(i,j)
      d_lon=lonr(i+1,j)-lonr(i,j)
      d_lon=d_lon*cos(latr(i,j)/180.0d0*PI)
      cosAu(i,j) = d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAu(i,j) = d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo
  enddo
  cosAu(L,:) = cosAu(L-1,:)
  sinAu(L,:) = sinAu(L-1,:)

  do j=0,M-1
    do i=0,L
      d_lat=latr(i,j+1)-latr(i,j)
      d_lon=lonr(i,j)-lonr(i,j+1)
      d_lon=d_lon*cos(latr(i,j)/180.0d0*PI)
      cosAv(i,j) = d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAv(i,j) = d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo         
  enddo
  cosAv(:,M) = cosAv(:,M-1)
  sinAv(:,M) = sinAv(:,M-1)

  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  write (DD, "(I2.2)") Sday

!==== Read ATM file coordinates ==================================

#if defined NETCDF_INPUT
!---- NetCDF input file -------------------------------------------
# if defined JMA_MSM
!---- Read JMA-MSM NetCDF file --------------------------------
  allocate( IN_FILE(2) )
  IN_FILE(1) = trim(ATM_dir)//YYYY//"/"//MM//DD//".nc"
  IN_FILE(2) = trim(ATM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

# elif defined ERA5
!---- Read ERA5 NetCDF file --------------------------------
  allocate( IN_FILE(1) )
  IN_FILE(1) = trim(ATM_FILE(1))

# elif defined FORP_ATM
!---- Read ERA5 NetCDF file --------------------------------
  allocate( IN_FILE(N_InPar) )
  do iin=1,N_InPar
    IN_FILE(iin) = trim(ATM_dir)//"forp-jpn-v4_"//trim(NCIN_NAME(iin))// &
                   "_dy_"//YYYY//MM//".nc"
  end do

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
  allocate(amask(Im,Jm))
  allocate(alat(Im,Jm))
  allocate(alon(Im,Jm))
  
  call check( nf90_inq_varid(ncid, NC_LAT_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lat) )
  call check( nf90_inq_varid(ncid, NC_LON_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lon) )
  call check( nf90_close(ncid) )

  do j=1,Jm
    do i=1,Im
      alat(i,j) = lat(j)
      alon(i,j) = lon(i)
    enddo
  enddo

# if defined FORP_ATM
  allocate(amaskv(Im,Jm))
  allocate(alatv(Im,Jm))
  allocate(alonv(Im,Jm))
  write(*,*) "OPEN: ", trim( IN_FILE(2) )
  !Open NetCDF file
  call check( nf90_open(trim( IN_FILE(2) ), nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, NC_LAT_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lat) )
  call check( nf90_inq_varid(ncid, NC_LON_NAME, var_id) )
  call check( nf90_get_var(ncid, var_id, lon) )
  call check( nf90_close(ncid) )

  do j=1,Jm
    do i=1,Im
      alatv(i,j) = lat(j)
      alonv(i,j) = lon(i)
    enddo
  enddo
# endif

#else

!---- GRIB2 input file -------------------------------------------
  GRIB_yyyymmddhh = YYYY//MM//DD//"00"
# if defined JMA_MSM
!---- Set JMA-MSM GRIB2 file --------------------------------
  allocate( IN_FILE(2) )
  IN_FILE(1) = trim(ATM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
              GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix
  IN_FILE(2) = trim(ATM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
              GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix

# elif defined DSJRA55
!---- Set DSJRA55 GRIB2 file --------------------------------
  allocate( IN_FILE(2) )
  IN_FILE(1) = trim(ATM_dir)//"Hist/Daily/fcst_surf/"//YYYY//MM// &
                 "/fcst_surf."//GRIB_yyyymmddhh
  IN_FILE(2) = trim(ATM_dir)//"Hist/Daily/fcst_phy2m/"//YYYY//MM// &
                 "/fcst_phy2m."//GRIB_yyyymmddhh

# elif defined JMA_LSM
!---- Set JMA-LSM GRIB2 file --------------------------------
  allocate( IN_FILE(1) )
  IN_FILE(1) = trim(ATM_dir)//YYYY//MM//"/LA"//YYYY//"_"//MM//"/"// &
              GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix

# elif defined JRA55
!---- Read JRA55 GRIB2 file --------------------------------
  allocate( IN_FILE(2) )
  IN_FILE(1) = trim(ATM_dir)//"Hist/Daily/fcst_surf125/"//YYYY//MM// &
                 "/fcst_surf125."//GRIB_yyyymmddhh
  IN_FILE(2) = trim(ATM_dir)//"Hist/Daily/fcst_phy2m125/"//YYYY//MM// &
                 "/fcst_phy2m125."//GRIB_yyyymmddhh

# elif defined ERA5
!---- Set ERA5 GRIB2 file --------------------------------
  allocate( IN_FILE(1) )
  IN_FILE(1) = trim(ATM_FILE(1))
# endif

  !Open GRIB file
  call codes_grib_multi_support_on	(	iret(1)	)	
  write(*,*) "OPEN: ", trim( IN_FILE(1) )
  call codes_open_file(ifile(1), trim( IN_FILE(1) ),'r', iret(1))
  call codes_grib_new_from_file(ifile(1),igrib(1), iret(1))
!
# if defined JMA_MSM || defined DSJRA55 || defined JMA_LSM || defined JRA55
! ! Get dimension data
  call codes_get(igrib(1),'Ny', Jm)
  call codes_get(igrib(1),'Nx', Im)
# elif defined ERA5
! ! Get dimension data
  call codes_get(igrib(1),'Nj', Jm)
  call codes_get(igrib(1),'Ni', Im)
# endif
  write(*,*) Im,Jm
! Allocate variable
  allocate(values(Im*Jm))

# if defined DSJRA55 || defined JMA_LSM
  allocate(lat(Im,Jm))
  allocate(lon(Im,Jm))
  allocate(cosAx(Im,Jm))
  allocate(sinAx(Im,Jm))
  allocate(cosAy(Im,Jm))
  allocate(sinAy(Im,Jm))

#  if defined DSJRA55
!  ---- Get Lat Lon coordinates from Binary file --------------
  LL_FILE = trim(ATM_dir)//"Consts/Lambert5km_latlon.dat"
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

#  elif defined JMA_LSM
!  ---- Calculate Lambert Lat Lon coordinates --------------
!        Please see "users_guide_local_analisys.pdf" (p.2) 
  do j=1, Jm
    do i=1,Im
      xp = (dble(i)-449.0d0)*5000.0d0
      yp = 7.71061d6+(dble(j)-361.0d0)*5000.0d0
      lat(i,j) = 90.0d0-2.0d0 * atan( 1.37003d-10*(xp*xp+yp*yp)**0.69875d0 )*rad2deg
      lon(i,j) = 1.39750d0 * atan(xp/yp)*rad2deg + 140
    end do
  end do

#  endif

  write(*,*) 'NW corner:', lat(1,1),  lon(1,1)
  write(*,*) 'SW corner:', lat(1,Jm), lon(1,Jm)
  write(*,*) 'NE corner:', lat(Im,1), lon(Im,1)
  write(*,*) 'SE corner:', lat(Im,Jm),lon(Im,Jm)
   
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

# elif defined JMA_MSM || defined JRA55
!  ---- Get Lat Lon coordinates from GRIB file --------------
  call codes_get(igrib(1),'distinctLatitudes',lat)
  call codes_get(igrib(1),'distinctLongitudes',lon)

# elif defined ERA5
!  ---- Get Lat Lon coordinates from GRIB file --------------
  allocate(lat(Jm))
  allocate(lon(Im))
  call codes_get(igrib(1),'latitudeOfFirstGridPointInDegrees',  s_lat)
  call codes_get(igrib(1),'longitudeOfFirstGridPointInDegrees', s_lon)
  call codes_get(igrib(1),'jDirectionIncrementInDegrees', d_lat)
  call codes_get(igrib(1),'iDirectionIncrementInDegrees', d_lon)
  do j=1,Jm
    lat(j) = s_lat - d_lat*dble(j-1)
  end do
  do i=1,Im
    lon(i) = s_lon + d_lon*dble(i-1)
  end do

  write(*,*) 'NW corner:', lat(1),  lon(1)
  write(*,*) 'SW corner:', lat(Jm), lon(1)
  write(*,*) 'NE corner:', lat(1), lon(Im)
  write(*,*) 'SE corner:', lat(Jm),lon(Im)

# endif 

  call codes_release(igrib(1))
  call codes_close_file(ifile(1))
  write(*,*) "CLOSE: ", trim( IN_FILE(1) )
#endif

  allocate(in_data(Im, Jm, N_InPar))
  allocate(in_data2(Im, Jm, N_OutPar))

# if defined FORP_ATM
! create atm land mask grid
  amask = 1.0d0
  amaskv = 1.0d0

  write(*,*) "OPEN: ", trim( IN_FILE(1) )
  !Open NetCDF file
  call check( nf90_open(trim( IN_FILE(1) ), nf90_nowrite, ncid) )
  start3D = (/ 1,  1,  1 /)
  count3D = (/ Im, Jm, 1 /)
  call check( nf90_inq_varid(ncid, trim( NCIN_NAME(1) ), var_id) ) 
  call check( nf90_get_var(ncid, var_id, in_data(:,:,1), start=start3D, count=count3D) )
  call check( nf90_close(ncid) )     
  
  write(*,*) "CREATE: FORP atm T points land mask"
  do j=1, Jm
    do i=1,Im
      if(in_data(i,j,1)<-1.0d10) then
        amask(i,j)=0.0d0
      endif
    end do
  end do

  write(*,*) "OPEN: ", trim( IN_FILE(2) )
  !Open NetCDF file
  call check( nf90_open(trim( IN_FILE(2) ), nf90_nowrite, ncid) )
  start3D = (/ 1,  1,  1 /)
  count3D = (/ Im, Jm, 1 /)
  call check( nf90_inq_varid(ncid, trim( NCIN_NAME(2) ), var_id) ) 
  call check( nf90_get_var(ncid, var_id, in_data(:,:,2), start=start3D, count=count3D) )
  call check( nf90_close(ncid) )     
  
  write(*,*) "CREATE: FORP atm V points land mask"
  do j=1, Jm
    do i=1,Im
      if(in_data(i,j,2)<-1.0d10) then
        amaskv(i,j)=0.0d0
      endif
    end do
  end do
  Irdg_min = 1
  Irdg_max = Im
  Jrdg_min = 1
  Jrdg_max = Jm

# endif

!==== Calculate weight parameters for interpolation ====================
!  ---- Calc. weight parameters for interpolation --------------  
  write(*,*) "CALC.: weight parameters for interpolation"
#if defined DSJRA55 || defined JMA_LSM
  call weight2D_grid2(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)
#elif defined FORP_ATM
!  call weight2D_grid2(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)
  write(*,*) "CALC.: weight parameters for T points"
  call weight2D_grid3_2(                              &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
          , alon(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
          , alat(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
          , amask(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
          , 0, L, 0, M                        &
          , lonr(0:L,0:M)                     &
          , latr(0:L,0:M)                     &
          , ID_cont,w_cont )
  write(*,*) "CALC.: weight parameters for V points"
  call weight2D_grid3_2(                              &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
          , alonv(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
          , alatv(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
          , amaskv(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
          , 0, L, 0, M                        &
          , lonr(0:L,0:M)                     &
          , latr(0:L,0:M)                     &
          , ID_contv,w_contv )
#else
  call weight2D_grid(Im,Jm,lon,lat,Nxr,Nyr,lonr,latr,ID_cont,w_cont)  
#endif

!==== Create the forcing netCDF file ===================================

  FRC_yyyymmdd(2:5)=YYYY
  FRC_yyyymmdd(6:7)=MM
  FRC_yyyymmdd(8:9)=DD

#if defined JMA_MSM_CLOUD_ONLY
  FRC_FILE(1) = trim(FRC_prefix)//FRC_yyyymmdd//'_cloud.nc'
#else
  FRC_FILE(1) = trim(FRC_prefix)//FRC_yyyymmdd//'_1.nc'
#endif  
  write(*,*) "CREATE: ", trim( FRC_FILE(1) )
  call check( nf90_create(trim( FRC_FILE(1) ), nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xi_rho_dimid) )
  call check( nf90_def_dim(ncid, 'eta_rho',Nyr, eta_rho_dimid) )
  call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid) )
  dimids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)
  call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, var_id) )
  call check( nf90_put_att(ncid, var_id, 'long_name', 'atmospheric forcing time') )
  call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

  DO iparam=1,N_OutPar1
#if defined JMA_MSM_CLOUD_ONLY
    if(iparam<6.or.iparam>9) cycle
#endif
    call check( nf90_def_var(ncid, trim( NC_NAME(iparam) ), NF90_REAL, dimids, var_id) )
    call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(iparam) )) )
    call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(iparam) ) ) )
    call check( nf90_put_att(ncid, var_id, 'time',      'time') )
  END DO
! End define mode.
  call check( nf90_enddef(ncid) )
  call check( nf90_close(ncid) )

#if !defined JMA_LSM  && !defined JMA_MSM_CLOUD_ONLY
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

  DO iparam=N_OutPar1+1,N_OutPar
    call check( nf90_def_var(ncid, trim( NC_NAME(iparam) ), NF90_REAL, dimids, var_id) )
    call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(iparam) )) )
    call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(iparam) ) ) )
    call check( nf90_put_att(ncid, var_id, 'time',      'time') )
  END DO
! End define mode.
  call check( nf90_enddef(ncid) )
  call check( nf90_close(ncid) )
#endif

!==== LOOP set up ==========================================
#if defined JMA_MSM
  call jd(2006, 3, 1, jd_msmnew)   ! MSM input-format regime switch (2006-03-01)
#endif

#if defined ERA5
!======= Set starting time and Ending time (ERA5) =================
!  Build INFILE(:) (file names + per-file step times) then let the shared
!  mod_infile::infile_check_time produce Nt / atm_time / iNCt / idt. This
!  replaces the old hand-rolled NC(:) time-list build.
  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)

  jdate_End = jdate_Start + N_days

  write(*,*) jdate_Start, jdate_End, N_days

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref

  call jd(1900, 1, 1, jdate_Ref_atm)   ! ERA5 time: hours since 1900-01-01 00:00:00
  d_jdate_Ref_atm = dble(jdate_Ref_atm)
  write(*,*) d_jdate_Ref_atm

  allocate( INFILE(NCnum) )

# if defined NETCDF_INPUT
!== ERA5 NetCDF format ========================
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    write(*,*) "OPEN: ", trim( ATM_FILE(iNC) )
    call check( nf90_open(trim( ATM_FILE(iNC) ), nf90_nowrite, ncid) )
    call get_dimension(ncid, NC_TIME_NAME, INFILE(iNC)%Nt)
    write(*,*) INFILE(iNC)%Nt
    allocate( INFILE(iNC)%NAME(1) )
    INFILE(iNC)%NAME(1) = ATM_FILE(iNC)
    allocate( INFILE(iNC)%time_all(INFILE(iNC)%Nt) )
    allocate( time2(INFILE(iNC)%Nt) )
    call check( nf90_inq_varid(ncid, NC_TIME_NAME, var_id) )
    call check( nf90_get_var(ncid, var_id, time2) )
    call check( nf90_close(ncid) )

    INFILE(iNC)%time_all = time2/24.0d0 + d_jdate_Ref_atm ! hour -> day

    deallocate(time2)
  end do

# else
!== ERA5 GRIB format ========================
  allocate( time2(100000) ) ! temporary array
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    write(*,*) "OPEN: ", trim( ATM_FILE(iNC) )
    call codes_open_file(ifile(1), trim( ATM_FILE(iNC) ),'r', iret(1))
    call codes_grib_new_from_file(ifile(1),igrib(1), iret(1))
    allocate( INFILE(iNC)%NAME(1) )
    INFILE(iNC)%NAME(1) = ATM_FILE(iNC)
    INFILE(iNC)%Nt=0
    do while (iret(1) /= CODES_END_OF_FILE)
      call codes_get(igrib(1),'dataDate' , p1)
      call codes_get(igrib(1),'dataTime' , p2)
      call codes_get(igrib(1),'shortName', p4)
      if (trim(p4)==trim(GRIB_NAME(1))  ) then
        write(*,*) 'CHECK: Time', p1, p2
        INFILE(iNC)%Nt = INFILE(iNC)%Nt + 1
        call dble_jd2(p1, p2, time2( INFILE(iNC)%Nt ))
      endif
      call codes_release(igrib(1))
      call codes_grib_new_from_file(ifile(1),igrib(1), iret(1))
    enddo

    write(*,*) INFILE(iNC)%Nt
    allocate( INFILE(iNC)%time_all(INFILE(iNC)%Nt) )

    INFILE(iNC)%time_all = time2( 1:INFILE(iNC)%Nt ) ! julian date

    call codes_close_file(ifile(1))
  end do
  deallocate(time2)
# endif

  call infile_check_time( NCnum, dble(jdate_Start), dble(jdate_End), Nt, atm_time, iNCt, idt )

  write(*,*) "*************************************"

  iNC = 0

#elif defined FORP_ATM
!======= Set starting time and Ending time (FORP) =================
!  Build INFILE(:) the same way as ERA5: each monthly file group holds one file
!  per input variable (NAME(1)=pso is the time source); infile_check_time then
!  yields Nt / atm_time / iNCt / idt.
  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)

  jdate_End = jdate_Start + N_days

  write(*,*) jdate_Start, jdate_End, N_days

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref

  call jd(1, 1, 1, jdate_Ref_atm)! FORP time: hours since 1-01-01 00:00:00
  d_jdate_Ref_atm = dble(jdate_Ref_atm)-2.0d0
  write(*,*) d_jdate_Ref_atm

  NCnum = 12*(Eyear-Syear) + Emonth - Smonth + 1
  NCnum = min( NCnum, 12 )  ! FORP data limit to 1 year
  write(*,*) 'NCnum = ', NCnum
  allocate( ATM_FILE(NCnum) )
  allocate( INFILE(NCnum) )
  do iNC=1, NCnum
    write (YYYY, "(I4.4)") Syear
    write (MM, "(I2.2)") Smonth + iNC -1
    allocate( INFILE(iNC)%NAME(N_InPar) )
    do iin=1, N_InPar
      INFILE(iNC)%NAME(iin) = trim(ATM_dir)//"forp-jpn-v4_"//trim(NCIN_NAME(iin))// &
                     "_dy_"//YYYY//MM//".nc"
    end do
    ATM_FILE(iNC) = INFILE(iNC)%NAME(1)
  end do

! Check time
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    ! Open NetCDF file
    write(*,*) "OPEN: ", trim( INFILE(iNC)%NAME(1) )
    call check( nf90_open(trim( INFILE(iNC)%NAME(1) ), nf90_nowrite, ncid) )
    call get_dimension(ncid, NC_TIME_NAME, INFILE(iNC)%Nt)
    write(*,*) INFILE(iNC)%Nt
    allocate( INFILE(iNC)%time_all(INFILE(iNC)%Nt) )
    allocate( time2(INFILE(iNC)%Nt) )
    call check( nf90_inq_varid(ncid, NC_TIME_NAME, var_id) )
    call check( nf90_get_var(ncid, var_id, time2) )
    call check( nf90_close(ncid) )

    INFILE(iNC)%time_all = time2/24.0d0 + d_jdate_Ref_atm ! hour -> day

    deallocate(time2)
  end do

  call infile_check_time( NCnum, dble(jdate_Start), dble(jdate_End), Nt, atm_time, iNCt, idt )

  write(*,*) "*************************************"

  iNC = 0

#elif defined JMA_MSM && defined NETCDF_INPUT
!======= Set starting time and Ending time (JMA-MSM NetCDF) =================
!  One daily file (24 steps) per day; rain comes from a parallel r1h file.
!  Enumerate candidate days [Start, End-1], skip missing days, and build
!  INFILE(:) (NAME(1)=main, NAME(2)=r1h; time_all = julian date of each step).
!  The per-file time variable is "hours since <that day> 00:00", so
!  time_all = (julian day of the file) + hours/24.
  call jd(Syear, Smonth, Sday, jdate_Start)
  call jd(Eyear, Emonth, Eday, jdate_End)
  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) jdate_Start, jdate_End

  allocate( INFILE(jdate_End - jdate_Start + 1) )   ! upper bound on day count
  NCnum = 0
  do ijdate = jdate_Start, jdate_End-1
    call cdate( ijdate, iyear, imonth, iday )
    write (YYYY,"(I4.4)") iyear
    write (MM,"(I2.2)") imonth
    write (DD,"(I2.2)") iday
    IN_FILE(1) = trim(ATM_dir)//YYYY//"/"//MM//DD//".nc"
    inquire(FILE=trim(IN_FILE(1)), EXIST=file_exists)
    if(.not. file_exists) then
      write(*,*) "Not found (skip): ", trim(IN_FILE(1))
      cycle
    endif
    NCnum = NCnum + 1
    allocate( INFILE(NCnum)%NAME(2) )
    INFILE(NCnum)%NAME(1) = trim(ATM_dir)//YYYY//"/"//MM//DD//".nc"
    INFILE(NCnum)%NAME(2) = trim(ATM_dir)//"r1h/"//YYYY//"/"//MM//DD//".nc"

    write(*,*) "OPEN: ", trim( INFILE(NCnum)%NAME(1) )
    call check( nf90_open(trim( INFILE(NCnum)%NAME(1) ), nf90_nowrite, ncid) )
    call get_dimension(ncid, NC_TIME_NAME, INFILE(NCnum)%Nt)
    allocate( INFILE(NCnum)%time_all(INFILE(NCnum)%Nt) )
    allocate( time2(INFILE(NCnum)%Nt) )
    call check( nf90_inq_varid(ncid, NC_TIME_NAME, var_id) )
    call check( nf90_get_var(ncid, var_id, time2) )
    call check( nf90_close(ncid) )

    INFILE(NCnum)%time_all = dble(ijdate) + time2/24.0d0   ! hour-of-day -> julian date
    deallocate(time2)
  end do
  write(*,*) 'NCnum = ', NCnum

  call infile_check_time( NCnum, dble(jdate_Start), dble(jdate_End), Nt, atm_time, iNCt, idt )

  write(*,*) "*************************************"

  iNC = 0

#elif defined JMA_MSM
!======= Set starting time and Ending time (JMA-MSM GRIB) =================
!  One 3-hourly GPV file (FH00-15) per base time; each base contributes 3 hourly
!  outputs (forecast hours 0,1,2). Enumerate base times Start..(End, exclusive)
!  every 3 h, skip missing files (unified skip; no FH-lookback gap fill), and
!  build INFILE(:) with NAME(1)=NAME(2)=file, Nt=3, time_all = base + {0,1,2}/24.
!  The validity time is base+FH, so times come straight from the filename --
!  no need to open the GRIB here.
  call jd(Syear, Smonth, Sday, jdate_Start)
  call jd(Eyear, Emonth, Eday, jdate_End)
  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) jdate_Start, jdate_End

  allocate( INFILE( (jdate_End - jdate_Start)*8 + 8 ) )   ! upper bound (8 files/day)
  NCnum = 0
  ihours = 0
  do while ( dble(jdate_Start) + dble(ihours)/24.0d0 < dble(jdate_End) )
    ijdate = jdate_Start + ihours/24
    call cdate( ijdate, iyear, imonth, iday )
    write (YYYY,"(I4.4)") iyear
    write (MM,"(I2.2)") imonth
    write (DD,"(I2.2)") iday
    write (hh,"(I2.2)") mod(ihours,24)
    GRIB_yyyymmddhh = YYYY//MM//DD//hh
    IN_FILE(1) = trim(ATM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
                 GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix
    inquire(FILE=trim(IN_FILE(1)), EXIST=file_exists)
    if(file_exists) then
      NCnum = NCnum + 1
      allocate( INFILE(NCnum)%NAME(2) )
      INFILE(NCnum)%NAME(1) = IN_FILE(1)
      INFILE(NCnum)%NAME(2) = IN_FILE(1)   ! MSM reads both halves from one file
      INFILE(NCnum)%Nt = 3                 ! forecast hours 0,1,2
      allocate( INFILE(NCnum)%time_all(3) )
      do i=1,3
        INFILE(NCnum)%time_all(i) = dble(ijdate) + dble(mod(ihours,24) + (i-1))/24.0d0
      end do
    else
      write(*,*) "Not found (skip): ", trim(IN_FILE(1))
    endif
    ihours = ihours + 3   !!! files exist at 3-hourly interval
  end do
  write(*,*) 'NCnum = ', NCnum

  call infile_check_time( NCnum, dble(jdate_Start), dble(jdate_End), Nt, atm_time, iNCt, idt )

  write(*,*) "*************************************"

  iNC = 0

#elif defined DSJRA55 || defined JRA55 || defined JMA_LSM
!======= Set starting time and Ending time (DSJRA55 / JRA55 / JMA_LSM GRIB) =====
!  One GRIB time-slice per file. Enumerate nominal times Start..(End, exclusive)
!  at the donor interval, skip missing files, and build INFILE(:) with Nt=1 and
!  time_all(1) = julian date of the file's nominal (= validity) time. The field
!  reads below select records by parameter (not by time), so the field data is
!  independent of this list; only the output time axis uses it.
  call jd(Syear, Smonth, Sday, jdate_Start)
  call jd(Eyear, Emonth, Eday, jdate_End)
  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) jdate_Start, jdate_End

# if defined JRA55
  i = 3      ! JRA55 files are 3-hourly
# else
  i = 1      ! DSJRA55 / JMA_LSM files are hourly
# endif
  allocate( INFILE( (jdate_End - jdate_Start)*24/i + 24 ) )   ! upper bound
  NCnum = 0
  ihours = 0
  do while ( dble(jdate_Start) + dble(ihours)/24.0d0 < dble(jdate_End) )
    ijdate = jdate_Start + ihours/24
    call cdate( ijdate, iyear, imonth, iday )
    write (YYYY,"(I4.4)") iyear
    write (MM,"(I2.2)") imonth
    write (DD,"(I2.2)") iday
    write (hh,"(I2.2)") mod(ihours,24)
    GRIB_yyyymmddhh = YYYY//MM//DD//hh
# if defined DSJRA55
    IN_FILE(1) = trim(ATM_dir)//"Hist/Daily/fcst_surf/"//YYYY//MM//"/fcst_surf."//GRIB_yyyymmddhh
    IN_FILE(2) = trim(ATM_dir)//"Hist/Daily/fcst_phy2m/"//YYYY//MM//"/fcst_phy2m."//GRIB_yyyymmddhh
# elif defined JRA55
    IN_FILE(1) = trim(ATM_dir)//"Hist/Daily/fcst_surf125/"//YYYY//MM//"/fcst_surf125."//GRIB_yyyymmddhh
    IN_FILE(2) = trim(ATM_dir)//"Hist/Daily/fcst_phy2m125/"//YYYY//MM//"/fcst_phy2m125."//GRIB_yyyymmddhh
# elif defined JMA_LSM
    IN_FILE(1) = trim(ATM_dir)//YYYY//MM//"/LA"//YYYY//"_"//MM//"/"// &
                 GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix
# endif
    inquire(FILE=trim(IN_FILE(1)), EXIST=file_exists)
    if(file_exists) then
      NCnum = NCnum + 1
# if defined JMA_LSM
      allocate( INFILE(NCnum)%NAME(1) )
      INFILE(NCnum)%NAME(1) = IN_FILE(1)
# else
      allocate( INFILE(NCnum)%NAME(2) )
      INFILE(NCnum)%NAME(1) = IN_FILE(1)
      INFILE(NCnum)%NAME(2) = IN_FILE(2)
# endif
      INFILE(NCnum)%Nt = 1
      allocate( INFILE(NCnum)%time_all(1) )
      INFILE(NCnum)%time_all(1) = dble(ijdate) + dble(mod(ihours,24))/24.0d0
    else
      write(*,*) "Not found (skip): ", trim(IN_FILE(1))
    endif
    ihours = ihours + i
  end do
  write(*,*) 'NCnum = ', NCnum

  call infile_check_time( NCnum, dble(jdate_Start), dble(jdate_End), Nt, atm_time, iNCt, idt )

  write(*,*) "*************************************"

  iNC = 0
#endif

!===== LOOP1 START ================================================
  do itime=1,Nt
    iNCm = iNC
    iNC  = iNCt(itime)
    ifc  = idt(itime)
#if defined JMA_MSM
    ijdate = int( INFILE(iNC)%time_all(ifc) )   ! julian day of this step (jd_msmnew test)
#endif

#if defined JMA_MSM && defined NETCDF_INPUT
    ! iNCm/iNC/ifc/ijdate already set at the top of the do itime loop.
    if(iNCm < iNC) then        ! donor-file change: repoint main + r1h files
      IN_FILE(1) = trim( INFILE(iNC)%NAME(1) )
      IN_FILE(2) = trim( INFILE(iNC)%NAME(2) )
    endif

#elif defined ERA5
    ! iNCm/iNC/ifc already set at the top of the do itime loop.
    if(iNCm < iNC) then        ! donor-file change: (re)point / (re)open
      IN_FILE(1) = trim(ATM_FILE(iNC))
# if !defined NETCDF_INPUT
      if(iNCm >= 1) call codes_close_file(ifile(1))   ! close previous GRIB file
      call codes_grib_multi_support_on( iret(1) )
      write(*,*) "OPEN: ", trim( IN_FILE(1) )
      call codes_open_file(ifile(1), trim( IN_FILE(1) ),'r', iret(1))
      call codes_grib_new_from_file(ifile(1),igrib(1), iret(1))
# endif
    endif

#elif defined FORP_ATM
    ! iNCm/iNC/ifc already set at the top of the do itime loop.
    if(iNCm < iNC) then        ! donor-file change: repoint the per-variable files
      do iin=1, N_InPar
        IN_FILE(iin) = trim( INFILE(iNC)%NAME(iin) )
      end do
    endif

#elif defined JMA_MSM
    ! GRIB: iNCm/iNC/ifc/ijdate already set at the top of the do itime loop.
    ! Each 3-hourly file is read through two handles (params 1-9 / 10-N), so open
    ! both on a file change and keep them across the file's 3 forecast steps.
    if(iNCm < iNC) then
      IN_FILE(1) = trim( INFILE(iNC)%NAME(1) )
      IN_FILE(2) = trim( INFILE(iNC)%NAME(2) )
      if(iNCm >= 1) then
        call codes_close_file(ifile(1))
        call codes_close_file(ifile(2))
      endif
      DO iin=1,2
        call codes_grib_multi_support_on( iret(iin) )
        write(*,*) "OPEN: ", trim( IN_FILE(iin) )
        call codes_open_file(ifile(iin), trim( IN_FILE(iin) ),'r', iret(iin))
        call codes_grib_new_from_file(ifile(iin),igrib(iin), iret(iin))
      END DO
    endif

#elif defined DSJRA55 || defined JRA55 || defined JMA_LSM
    ! GRIB: iNCm/iNC/ifc set at the top. One time-slice per file (Nt=1), so a new
    ! file is opened every step; DSJRA55/JRA55 read two files (surf / phy2m),
    ! JMA_LSM one.
    if(iNCm < iNC) then
# if defined JMA_LSM
      IN_FILE(1) = trim( INFILE(iNC)%NAME(1) )
      if(iNCm >= 1) call codes_close_file(ifile(1))
      DO iin=1,1
# else
      IN_FILE(1) = trim( INFILE(iNC)%NAME(1) )
      IN_FILE(2) = trim( INFILE(iNC)%NAME(2) )
      if(iNCm >= 1) then
        call codes_close_file(ifile(1))
        call codes_close_file(ifile(2))
      endif
      DO iin=1,2
# endif
        call codes_grib_multi_support_on( iret(iin) )
        write(*,*) "OPEN: ", trim( IN_FILE(iin) )
        call codes_open_file(ifile(iin), trim( IN_FILE(iin) ),'r', iret(iin))
        call codes_grib_new_from_file(ifile(iin),igrib(iin), iret(iin))
      END DO
    endif
#endif

! ======= NetCDF file input ================================================
#if defined NETCDF_INPUT
!  ---- LOOP3.1 START --------------------------------
      DO iparam=1,N_InPar
# if defined JMA_MSM
        if(iparam==N_OutPar1 .and. ijdate<jd_msmnew) cycle
#  if defined JMA_MSM_CLOUD_ONLY
        if(iparam<6.or.iparam>9) cycle
#  endif
        if(iparam>N_OutPar1) then !!! for rain (rain fall rate)
          iin=2
        else
          iin=1
        end if
# elif defined ERA5
        iin=1
# elif defined FORP_ATM
        iin = iparam
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
# elif defined FORP_ATM
        call check( nf90_get_var(ncid, var_id, in_data(:,:,iparam), start=start3D, count=count3D) )
# endif
        call check( nf90_close(ncid) )     
      END DO
!  ---- LOOP3.1 END --------------------------------   
#else

! ======= GRIB2 file input ================================================
      param_count = 0  ! This countor is only used for GRIB input.
# if defined JMA_MSM || defined DSJRA55 || defined JRA55
      DO iin=1,2
        if(iin==1) then
          ips = 1
          ipe = 9
        else
          ips = 10
          ipe = N_InPar
        end if
# else
      DO iin=1,1
        ipe = N_InPar
# endif
        DO WHILE (param_count/=ipe)

          SEARCH_LOOP: DO WHILE (iret(iin) /= CODES_END_OF_FILE)
# if defined JMA_MSM
            call codes_get(igrib(iin),'forecastTime',p1)
            call codes_get(igrib(iin),'parameterName', p4)
!            write(*,*) 'DEBUG1 :', iin, p1, trim(p4)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=ips,ipe
              if (p1==GRIB_STEP(ifc) .and.             &
                  trim(p4)==trim(GRIB_NAME(i))  ) then
                iparam = i
                param_count = param_count + 1
!                write(*,*) 'DEBUG2 :', iin, iparam, param_count, ifc  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                exit SEARCH_LOOP
              endif
            end do
# elif defined DSJRA55
            call codes_get(igrib(iin),'parameterCategory',p1)
            call codes_get(igrib(iin),'parameterNumber', p2)
            call codes_get(igrib(iin),'typeOfFirstFixedSurface', p3)
            call codes_get(igrib(iin),'parameterName', p4)
            do i=ips,ipe
              if (p1==GRIB_NUM(1,i) .and. p2==GRIB_NUM(2,i) .and. &
                  p3==GRIB_NUM(3,i)   ) then
                iparam = i
                param_count = param_count + 1
                exit SEARCH_LOOP
              endif
            end do
# elif defined JMA_LSM
            call codes_get(igrib(iin),'level',p1)
            call codes_get(igrib(iin),'shortName', p4)
            do i=1,N_InPar
              if (p1==GRIB_LEVEL(i) .and.             &
                trim(p4)==trim(GRIB_NAME(i))  ) then
                iparam = i
                param_count = param_count + 1
                exit SEARCH_LOOP
              endif
            end do
# elif defined JRA55
            call codes_get(igrib(iin),'indicatorOfParameter',p1)
            call codes_get(igrib(iin),'parameterName', p4)
            if(iin==1) then
              ips = 1
              ipe = 9
            else
              ips = 10
              ipe = N_InPar
            end if
            do i=1,N_InPar
              if (p1==GRIB_NUM(i) ) then
                iparam = i
                param_count = param_count + 1
                exit SEARCH_LOOP
              endif
            end do
# elif defined ERA5
            call codes_get(igrib(iin),'dataDate' , p1)
            call codes_get(igrib(iin),'dataTime' , p2)
            call codes_get(igrib(iin),'endStep'  , p3)
            call codes_get(igrib(iin),'shortName', p4)
            call dble_jd2( p1, p2, d_jdate)
            d_jdate = d_jdate + dble(p3)/24.0d0
            if (d_jdate == INFILE(iNC)%time_all(ifc) ) then
              do i=1,N_InPar
                if (trim(p4)==trim(GRIB_NAME(i))  ) then
                  iparam = i
                  param_count = param_count + 1
                  exit SEARCH_LOOP
                endif
              end do
            endif
# endif
            call codes_release(igrib(iin))
            call codes_grib_new_from_file(ifile(iin),igrib(iin), iret(iin))
          END DO SEARCH_LOOP
   
          write(*,*) "READ GRIB DATA: ", trim(p4)

          call codes_get(igrib(iin),'values', values)
          call codes_release(igrib(iin))
          call codes_grib_new_from_file(ifile(iin),igrib(iin), iret(iin))

          do i=1, Jm
            istart = 1 + Im*(i-1)
            iend   = Im*i
            in_data(:,i,iparam) = values(istart:iend)
          end do
        END DO

      END DO
#endif

!  ---- LOOP3.1 END --------------------------------
#if defined JMA_MSM
# if !defined JMA_MSM_CLOUD_ONLY
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
# endif
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
# if !defined JMA_MSM_CLOUD_ONLY
    ! for rain (Total precipitation rate)
      in_data2(:,:,10) = in_data(:,:,10)/3600.0d0  ! kg m-2 h-1 -> kg m-2 s-1 
#  if defined SWRAD
    ! for short-wave radiation
      in_data2(:,:,11) = in_data(:,:,11)  ! W m-2 -> W m-2
#  endif
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
      in_data2(:,:,10) = in_data(:,:,10)  ! kg m-2 s-1 -> kg m-2 s-1     
# if defined BULK_FLUX
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,11) = -in_data(:,:,11)  ! W m-2 -> W m-2
      in_data2(:,:,12) = -in_data(:,:,12)  ! W m-2 -> W m-2
      in_data2(:,:,13) =  in_data(:,:,13) !-in_data(14,:,:) ! solar radiation = downward
      in_data2(:,:,14) =  in_data(:,:,14) ! downward long-wave
      in_data2(:,:,15) =  in_data(:,:,15)-in_data(:,:,16) ! net long-wave  = downward - upward
# endif
#elif defined JMA_LSM
    ! for Pair (Pressure) index 5->1
      in_data2(:,:,1) = in_data(:,:,5)*0.01  ! Pa -> millibar (= hPa) 
    ! for U V: change DSJRA55 Lambert conformal to regular Lat Lon coordinat vectors
      ! indecis: 1->2, 2-> 3 
      do j=1, Jm
        do i=1,Im
          in_data2(i,j,2) = in_data(i,j,1)*cosAx(i,j)-in_data(i,j,2)*sinAy(i,j)
          in_data2(i,j,3) = in_data(i,j,1)*sinAx(i,j)+in_data(i,j,2)*cosAy(i,j)
        end do
      end do
    ! for Tair index 3->4
      in_data2(:,:,4) = in_data(:,:,3) - 273.15d0  ! K -> degC
    ! for Qait (Relative humidity) index 4-> 5
      in_data2(:,:,5) = in_data(:,:,4)
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
      in_data2(:,:,10) = in_data(:,:,10)*1.0d0/86400.0d0  ! mm day-1 -> kg m-2 s-1  
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
      in_data2(:,:,6) = in_data(:,:,6)*1000.0d0/3600.0d0  ! m h-1 -> kg m-2 s-1  
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,7) = in_data(:,:,7)/3600.0d0  ! J h-1 m-2 -> W m-2
      in_data2(:,:,8) = in_data(:,:,8)/3600.0d0  ! J h-1 m-2 -> W m-2
! --------------------------------------------
#elif defined FORP_ATM
    ! for Pair (Pressure)
      in_data2(:,:,1) = in_data(:,:,1)  ! hPa -> millibar (= hPa) 
    ! for U V
      in_data2(:,:,2) = in_data(:,:,2)*0.01d0  ! cm s-1 -> m s-1 
      in_data2(:,:,3) = in_data(:,:,3)*0.01d0  ! cm s-1 -> m s-1 
    ! for Tair
      in_data2(:,:,4) = in_data(:,:,4)  ! degC -> degC
    ! for Qair: Water Evaporation Flux (kg m-2 s-1) to Relative humidity (%)
      do j=1, Jm
        do i=1,Im
          !!! TO BE UPDATED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          in_data2(i,j,5) = in_data(i,j,5)
        end do
      end do
    ! for rain (Total precipitation rate)
      in_data2(:,:,6) = in_data(:,:,6)  ! kg m-2 s-1 -> kg m-2 s-1  
    ! for bulk flux parameters, heating -> positive, cooling -> negattive
      in_data2(:,:,7) = in_data(:,:,7)  ! W m-2 -> W m-2
      in_data2(:,:,8) = in_data(:,:,8)  ! W m-2 -> W m-2
#endif         
!  ---- LOOP3.2 START --------------------------------
      DO iparam=1,N_OutPar
#if defined JMA_MSM
# if defined JMA_MSM_CLOUD_ONLY
        if(iparam<6.or.iparam>9) cycle
# endif
#endif
        write(*,*) 'Linear Interporation: ',trim( NC_NAME(iparam) )

#if defined FORP_ATM
        if(iparam==2 .or. iparam==3) then
          call interp2D_grid3_2(                                              &
              Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, in_data2(:,:,iparam)  &
            , 0, L, 0, M                                      &
            , ID_contv, w_contv                               &
            , out_data(0:L,0:M,iparam)  )
        else
          call interp2D_grid3_2(                                              &
              Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, in_data2(:,:,iparam)  &
            , 0, L, 0, M                                      &
            , Id_cont, w_cont                                 &
            , out_data(0:L,0:M,iparam)  )
        endif
#else 
        call interp2D_grid (Im, Jm, in_data2(:,:,iparam)       &  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          , Nxr, Nyr, out_data(:,:,iparam)     &
                          , Id_cont, w_cont)
#endif
      END DO
!  ---- LOOP3.2 END --------------------------------
#if !defined JMA_MSM_CLOUD_ONLY
  !!! for U V: change regular Lat Lon to ROMS grid coordinat vectors 
      do j=0,M
        do i=0,L
          u = out_data(i,j,2)*cosAu(i,j)+out_data(i,j,3)*sinAu(i,j)
          v =-out_data(i,j,2)*sinAv(i,j)+out_data(i,j,3)*cosAv(i,j)
          out_data(i,j,2) = u
          out_data(i,j,3) = v
        enddo
      enddo
#endif
!  ---- LOOP3.3 START --------------------------------
    ! atm_time(itime) is a julian date; subtract d_jdate_Ref -> days since &refdate
      time(1) = atm_time(itime)-d_jdate_Ref
      write(*,*) time(1),TIME_ATT
      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d( 'time', trim( FRC_FILE(1) )                  &
            , 1, time, start1D, count1D )

#if !defined JMA_LSM && !defined JMA_MSM_CLOUD_ONLY
# if defined JRA55
      time(1) = atm_time(itime)-d_jdate_Ref+1.5d0/24.0d0
# elif defined ERA5
      time(1) = atm_time(itime)-d_jdate_Ref-0.5d0/24.0d0
# elif defined FORP_ATM
      time(1) = atm_time(itime)-d_jdate_Ref
# else
      time(1) = atm_time(itime)-d_jdate_Ref+0.5d0/24.0d0
# endif

      start1D = (/ itime /)
      count1D = (/ 1 /)
      call writeNetCDF_1d( 'time', trim( FRC_FILE(2) )                  &
            , 1, time, start1D, count1D )
#endif

      DO iparam=1,N_OutPar
#if defined JMA_MSM
# if defined JMA_MSM_CLOUD_ONLY
        if(iparam<6.or.iparam>9) cycle
# endif
#endif
        if(iparam>N_OutPar1) then !!! for rain (rain fall rate)
          iin=2
        else
          iin=1
        end if
  
        start3D = (/ 1,  1,  itime /)
        count3D = (/ Nxr, Nyr, 1 /)     
        call writeNetCDF_3d(trim( NC_NAME(iparam) ), trim( FRC_FILE(iin))   &
            , Nxr, Nyr, 1, out_data(:,:,iparam), start3D, count3D )
        
      END DO
!  ---- LOOP3.3 END --------------------------------
  end do
!---- LOOP1 END (do itime=1,Nt) --------------------------------
#if !defined NETCDF_INPUT
  write(*,*) "CLOSE: ", trim( IN_FILE(1) )
  if(iNC >= 1) call codes_close_file(ifile(1))   ! close last GRIB file
# if defined JMA_MSM || defined DSJRA55 || defined JRA55
  if(iNC >= 1) call codes_close_file(ifile(2))   ! second handle (MSM same file / surf+phy2m)
# endif
#endif

  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM frcATM2ROMS
      
