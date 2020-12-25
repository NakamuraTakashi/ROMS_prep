
!!!=== Copyright (c) 2014-2021 Takashi NAKAMURA  =====

PROGRAM getHYCOM
  use netcdf
  use mod_roms_netcdf
  use mod_calendar
 
  implicit none
! -------------------------------------------------------------------------
  real(8) :: Tlat, Blat, Llon, Rlon
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  integer :: Ryear, Rmonth, Rday
  character(256) :: HYCOM_prefix
  integer :: mode
! -------------------------------------------------------------------------
  character(31) :: TIME_ATT   = "hours since 2000-01-01 00:00:00"
  character(10) :: HYCOM_suffix = "_200001.nc"
  character(256) :: HYCOM_FILE
  
  TYPE T_NC
    real(8), pointer :: time_all(:)
    integer :: Nt
    integer :: ItS
    integer :: ItE
  END TYPE T_NC
  TYPE (T_NC) :: NC(NCnum)
  real(8), allocatable :: lat_all(:), lon_all(:)
  real(8), allocatable :: lat(:), lon(:), depth(:),tau(:)
  real(8), allocatable :: time2(:)
  real(8), allocatable :: water_u(:, :, :, :), water_v(:, :, :, :)
  real(8), allocatable :: salinity(:, :, :, :),water_temp(:, :, :, :)
  real(8), allocatable :: surf_el(:, :, :)
  integer :: start1D(1), count1D(1)
  integer :: start3D(3), count3D(3)
  integer :: start4D(4), count4D(4)
  
  integer :: N_days
  integer :: jdate_Start, jdate_End, jdate_Ref
  integer :: jdate_20000101
  integer :: iyear, imonth, iday
  integer :: i,j,k
  integer :: idays
!      integer :: jdate
  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(11) :: YYYYMMDDpHH
  real(8) :: d_jdate_20000101, d_jdate     
  integer :: ncid,var_id
  integer :: ncid2,var_id2
  integer :: N_xi_rho, N_eta_rho
  real(8) :: sf, off
  integer :: Im, Jm, Nz
  integer :: Im_all, Jm_all
  integer :: IL, IR, JB, JT
  integer :: iNC, itime
  integer :: end_flag
  integer :: create_flag = 0
  integer :: st = 1
  real(8) :: time(1)
  integer :: itry, status
  namelist/range/Tlat, Blat, Llon, Rlon
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/intpmode/mode
  namelist/hycom/HYCOM_prefix
      
  ! Read parameters in namelist file
  
  read (*, nml=range)
  read (*, nml=sdate)
  read (*, nml=edate)
  read (*, nml=refdate)
  read (*, nml=intpmode)
  read (*, nml=hycom)
  
  call jd(2000, 1, 1, jdate_20000101)
  d_jdate_20000101 = dble(jdate_20000101)
  
  write(*,*) d_jdate_20000101
  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)
  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  jdate_End = jdate_Start + N_days
  
  write(*,*) jdate_Start, jdate_End, N_days

!---- Modify time-unit description ---------------------------------
  
  write (YYYY, "(I4.4)") Ryear
  write (MM, "(I2.2)") Rmonth
  write (DD, "(I2.2)") Rday
  
  TIME_ATT(13:16)=YYYY
  TIME_ATT(18:19)=MM
  TIME_ATT(21:22)=DD
  
!---- Create the Output file --------------------------------
  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  HYCOM_suffix(2:5)=YYYY
  HYCOM_suffix(6:7)=MM
  HYCOM_OUT_FILE = trim( HYCOM_prefix )//HYCOM_suffix

!---- Allocate lat, lon dimension data ---------------------------------
  
  write(*,*) "OPEN: ", HYCOM_FILE(1)
  call try_nf_open(HYCOM_FILE(1), nf90_nowrite, ncid)
  
  ! Get dimension data
  call get_dimension(ncid, 'lat', Jm_all)
  call get_dimension(ncid, 'lon', Im_all)
  call get_dimension(ncid, 'depth', Nz)
  write(*,*) Jm_all, Im_all, Nz
      
!      ! Allocate variable
  allocate(lat_all(Jm_all))
  allocate(lon_all(Im_all))
  allocate(depth(Nz))
      
  call readNetCDF_1d(ncid, 'lat', Jm_all, lat_all)
  call readNetCDF_1d(ncid, 'lon', Im_all, lon_all)
  call readNetCDF_1d(ncid, 'depth', Nz, depth)
  call check( nf90_close(ncid) )         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST code
  write(*,*) "CLOSE: ", HYCOM_FILE(1)     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST code
  write(*,*) Im, lon_all(1), lon_all(Im_all)
  JB=1
  do j=1,Jm_all
    if(lat_all(j)>=Blat) exit
    JB=j
  end do
  JT=JB
  do j=JB,Jm_all
    if(lat_all(j)>Tlat) exit
    JT=j+1
  end do
  IL=1
  do i=IL,Im_all
    if(lon_all(i)>=Llon) exit
    IL=i
  end do
  IR=IL
  do i=IL,Im_all
    if(lon_all(i)>Rlon) exit
    IR=i+1
  end do
  write(*,*) JB, lat_all(JB), JT, lat_all(JT)
  write(*,*) IL, lon_all(IL), IR, lon_all(IR)
  
  Jm=JT-JB+1
  Im=IR-IL+1
  
  allocate(lat(Jm))
  allocate(lon(Im))
  
  lat(:) = lat_all(JB:JT)
  lon(:) = lon_all(IL:IR)
  write(*,*) Im,Jm,Nz
  
  ! Initialize count and start for reading data
  
  allocate(surf_el(Im,Jm,1))
  allocate(water_temp(Im,Jm,Nz,1))
  allocate(salinity(Im,Jm,Nz,1))
  allocate(water_u(Im,Jm,Nz,1))
  allocate(water_v(Im,Jm,Nz,1))
          
!---- Create the extracted HYCOM netCDF file --------------------------------
      
  call createNetCDF_HYCOM(  &
!      input parameters
      trim( HYCOM_OUT_FILE ) &
    , TIME_ATT               &  
    , Im, Jm, Nz, N_days+1   &   
    )
  
  call check( nf90_open( trim( HYCOM_OUT_FILE ), NF90_WRITE, ncid2 ) )

  start1D = (/ 1  /)
  count1D = (/ Im /) 
  call check( nf90_inq_varid(ncid2, 'lon', var_id2) )
  call check( nf90_put_var(ncid2, var_id2, lon, start = start1D, count = count1D) )

  start1D = (/ 1  /)
  count1D = (/ Jm /)
  call check( nf90_inq_varid(ncid2, 'lat', var_id2) )
  call check( nf90_put_var(ncid2, var_id2, lat, start = start1D, count = count1D) )

  start1D = (/ 1  /)
  count1D = (/ Nz /)
  call check( nf90_inq_varid(ncid2, 'depth', var_id2) )
  call check( nf90_put_var(ncid2, var_id2, depth, start = start1D, count = count1D) ) 

  call check( nf90_close(ncid2) )

!---- Read HYCOM netCDF file --------------------------------
  write(*,*) "******************************************************************"

#if defined GOFS_31
# if defined ANALYSIS_Y
  open(50, file='time_HYCOM_GOF_31_analysisY.dat')
# elif defined ANALYSIS
  open(50, file='time_HYCOM_GOF_31_analysis.dat')
# elif defined REANALYSIS
  open(50, file='time_HYCOM_GOF_31_reanalysis.dat')
# endif
#elif defined GOFS_30
# if defined ANALYSIS
  open(50, file='time_HYCOM_GOF_30_analysis.dat')
# elif defined REANALYSIS
  open(50, file='time_HYCOM_GOF_30_reanalysis.dat')
# endif
#endif
#if defined SKIP_CHECK_TIME
  write(*,*) 'READ: Time'
  do iNC=1, NCnum
    read(50,*) NC(iNC)%Nt
    write(*,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    read(50,*) NC(iNC)%time_all
  end do
#else
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    ! Open NetCDF file
    write(*,*) "OPEN: ", HYCOM_FILE(iNC)
    call try_nf_open(HYCOM_FILE(iNC), nf90_nowrite, ncid)
    call get_dimension(ncid, 'time', NC(iNC)%Nt)
    write(*,*) NC(iNC)%Nt
    write(50,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    allocate( time2(NC(iNC)%Nt) )
    call readNetCDF_1d(ncid, 'time', NC(iNC)%Nt, time2)
    call check( nf90_close(ncid) )
    write(*,*) "CLOSE: ", HYCOM_FILE(iNC)
    NC(iNC)%time_all = time2
    write(50,*) NC(iNC)%time_all
    deallocate(time2)
  end do
#endif
  close(50)   

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
    do i=NC(iNC)%Nt-1,1,-1
      d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24.0d0
      if(d_jdate < dble(jdate_End)) then
        write(*,*) '*** FOUND: Ending point @ HYCOM_FILE',iNC
        NC(iNC)%ItE=i+1
        exit
      endif
    end do
  end do
  write(*,*) NC(:)%ItE 
  
  do iNC=1,NCnum
    do i=2,NC(iNC)%ItE
      d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24.0d0
      if(d_jdate>dble(jdate_Start)) then
        write(*,*) '*** FOUND: Starting point @ HYCOM_FILE',iNC
        NC(iNC)%ItS=i-1
        exit
      endif
    end do
  end do
  write(*,*) NC(:)%ItS 
    
!==== LOOP start!! =====================================================

  do iNC=1,NCnum
    if(NC(iNC)%ItS==-1) then
      cycle
    end if
    ! Seek lon index IL
    write(*,*) "OPEN: ", HYCOM_FILE(iNC)
    call try_nf_open(HYCOM_FILE(iNC), nf90_nowrite, ncid)
    call readNetCDF_1d(ncid, 'lon', Im_all, lon_all)
    write(*,*) Im_all, lon_all(1), lon_all(Im_all)
    IL=1
    do i=IL,Im_all
      if(lon_all(i)>=Llon) exit
      IL=i
    end do
    write(*,*) IL, lon_all(IL)
    
    do itime=NC(iNC)%ItS, NC(iNC)%ItE
    
      time(1) = NC(iNC)%time_all(itime) + dble(jdate_20000101-jdate_Ref)*24.0d0
    
      write(*,*) "******************************************************************"

      CALL oceantime2cdate(time(1)*3600,Ryear,Rmonth,Rday,YYYYMMDDpHH)
      write(*,*) 'time = ', YYYYMMDDpHH
          
! --- Read NetCDF file ------------------------

      start3D = (/ IL, JB, itime /)
      count3D = (/ Im, Jm, 1  /)
      call readNetCDF_3d(  ncid, 'surf_el'      &
        , Im, Jm, 1, start3D, count3D           &
        , surf_el )
  
      start4D = (/ IL, JB, 1,  itime /)
      count4D = (/ Im, Jm, Nz, 1  /)
      call readNetCDF_4d_2( ncid, 'water_temp'  &
        , Im, Jm, Nz, 1, start4D, count4D       &
        , water_temp )
      call readNetCDF_4d_2( ncid, 'salinity'    &
        , Im, Jm, Nz, 1, start4D, count4D       &
        , salinity )
      call readNetCDF_4d_2( ncid, 'water_u'     &
        , Im, Jm, Nz, 1, start4D, count4D       &
        , water_u )
      call readNetCDF_4d_2( ncid, 'water_v'     &
        , Im, Jm, Nz, 1, start4D, count4D       &
        , water_v )
     
! --- Write NetCDF file ------------------------
      write(*,*) "Write: HYCOM data"

      call check( nf90_open( trim( HYCOM_OUT_FILE ), NF90_WRITE, ncid2 ) )
            
      start1D = (/ st  /)
      count1D = (/ 1 /)
      call check( nf90_inq_varid(ncid2, 'time', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, time, start = start1D, count = count1D) )   
          
      start3D = (/ 1,  1,  st /)
      count3D = (/ Im, Jm, 1 /)
      call check( nf90_inq_varid(ncid2, 'surf_el', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, surf_el, start = start3D, count = count3D) )   
           
      start4D = (/ 1,  1,  1,  st /)
      count4D = (/ Im, Jm, Nz, 1 /)
      call check( nf90_inq_varid(ncid2, 'water_temp', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, water_temp, start = start4D, count = count4D) )   
      call check( nf90_inq_varid(ncid2, 'salinity', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, salinity, start = start4D, count = count4D) )   
      call check( nf90_inq_varid(ncid2, 'water_u', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, water_u, start = start4D, count = count4D) )   
      call check( nf90_inq_varid(ncid2, 'water_v', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, water_v, start = start4D, count = count4D) )   

      call check( nf90_close(ncid2) ) 
      
      st = st+1
    
    end do
    call check( nf90_close(ncid) )
    write(*,*) "CLOSE: ", HYCOM_FILE(iNC)
    
    write(*,*) "******************************************************************"
  end do
!==== LOOP END ========================================================
  

  write(*,*) 'FINISH!!'
   
END PROGRAM getHYCOM
      
