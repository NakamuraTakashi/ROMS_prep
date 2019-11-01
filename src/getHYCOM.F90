
!!!=== Copyright (c) 2014-2018 Takashi NAKAMURA  =====

    PROGRAM getHYCOM
      use netcdf
      use mod_roms_netcdf
      use mod_calendar
     
      implicit none
      
#if defined GOFS_31
# if defined ANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2014-07-01 to present) -----
      integer, parameter :: NCnum   = 6
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3"  &  ! lon(4500): -180 - 179.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.2"  &  ! lon(4500): -180 - 179.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.8"  &  ! lon(4500): 0 - 359.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7"  &  ! lon(4500): -180 - 179.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9"  &  ! lon(4500): 0 - 359.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0"  &  ! lon(4500): 0 - 359.92
     &  /)
# elif defined REANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Renalysis (since 1994-01-01 to 2015-12-31) -----
      integer, parameter :: NCnum   = 22
      character(63) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1994"  &  ! lon(4500): -180 - 179.92
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1995"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1996"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1997"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1998"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1999"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2000"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2001"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2002"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2003"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2004"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2005"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2006"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2007"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2008"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2009"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2010"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2011"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2012"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2013"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2014"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015"  &  ! lon(4500): -180 - 179.92
     &  /)
# endif
#elif defined GOFS_30
# if defined ANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to present) -----
      integer, parameter :: NCnum   = 4
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2"  &
     &  /)
# elif defined REANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31) -----
      integer, parameter :: NCnum   = 2
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"  &
     &  /)
# endif
#endif
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
      
!      call cdate(julian_date, iyear, imonth, iday)

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
      HYCOM_FILE = trim( HYCOM_prefix )//HYCOM_suffix

!---- Allocate lat, lon dimension data ---------------------------------
      
      write(*,*) "OPEN: ", NC_FILE(1)
      call try_nf_open(NC_FILE(1), nf90_nowrite, ncid)
      
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
      write(*,*) "CLOSE: ", NC_FILE(1)     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST code
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
      
      call createNetCDF(         &
!      input parameters
     &    trim( HYCOM_FILE )     &
     &  , TIME_ATT               &  
     &  , Im, Jm, Nz, N_days+1   &   
     &  )
      start1D = (/ 1  /)
      count1D = (/ Im /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lon'                  &
     &  , trim( HYCOM_FILE )     &
     &  , Im                     &
     &  , lon                    &
     &  , start1D, count1D       &
     &  )
      start1D = (/ 1  /)
      count1D = (/ Jm /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lat'                  &
     &  , trim( HYCOM_FILE )     &
     &  , Jm                     &
     &  , lat                    &
     &  , start1D, count1D       &
     &  )

      start1D = (/ 1  /)
      count1D = (/ Nz /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'depth'                &
     &  , trim( HYCOM_FILE )     &
     &  , Nz                     &
     &  , depth                  &
     &  , start1D, count1D       &
     &  )
     
      
!---- Read HYCOM netCDF file --------------------------------
      write(*,*) "******************************************************************"

#if defined GOFS_31
# if defined ANALYSIS
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
        write(*,*) "OPEN: ", NC_FILE(iNC)
        call try_nf_open(NC_FILE(iNC), nf90_nowrite, ncid)
        call get_dimension(ncid, 'time', NC(iNC)%Nt)
        write(*,*) NC(iNC)%Nt
        write(50,*) NC(iNC)%Nt
        allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
        allocate( time2(NC(iNC)%Nt) )
        call readNetCDF_1d(ncid, 'time', NC(iNC)%Nt, time2)
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE(iNC)
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
            write(*,*) '*** FOUND: Ending point @ NC_FILE',iNC
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
            write(*,*) '*** FOUND: Starting point @ NC_FILE',iNC
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
        write(*,*) "OPEN: ", NC_FILE(iNC)
        call try_nf_open(NC_FILE(iNC), nf90_nowrite, ncid)
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
          
          call readNetCDF_3d(        &
!          input parameters
     &        ncid                   &
     &      , 'surf_el'              &
     &      , Im, Jm, 1              &
     &      , start3D, count3D       &
!          output parameters
     &      , surf_el                &
     &)
     
          start4D = (/ IL, JB, 1,  itime /)
          count4D = (/ Im, Jm, Nz, 1  /)

          call readNetCDF_4d_2(      &
!          input parameters
     &        ncid                   &
     &      , 'water_temp'           &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , water_temp             &
     &)
          call readNetCDF_4d_2(      &
!          input parameters
     &        ncid                   &
     &      , 'salinity'             &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , salinity               &
     &)
          call readNetCDF_4d_2(      &
!          input parameters
     &        ncid                   &
     &      , 'water_u'              &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , water_u                &
     &)
          call readNetCDF_4d_2(      &
!          input parameters
     &        ncid                   &
     &      , 'water_v'              &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , water_v                &
     &)
     
! --- Write NetCDF file ------------------------
            
          start1D = (/ st  /)
          count1D = (/ 1 /)
          
          call writeNetCDF_1d(       &
!          input parameters
     &        'time'                 &
     &      , trim( HYCOM_FILE )     &
     &      , 1                      &
     &      , time(1)                &
     &      , start1D, count1D       &
     &      )
          start3D = (/ 1,  1,  st /)
          count3D = (/ Im, Jm, 1 /)
          
          call writeNetCDF_3d(       &
!          input parameters
     &        'surf_el'              &
     &      , trim( HYCOM_FILE )     &
     &      , Im, Jm, 1              &
     &      , surf_el                &
     &      , start3D, count3D       &
     &      )
          
          start4D = (/ 1,  1,  1,  st /)
          count4D = (/ Im, Jm, Nz, 1 /)
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_temp'           &
     &      , trim( HYCOM_FILE )     &
     &      , Im, Jm, Nz, 1          &
     &      , water_temp             &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'salinity'             &
     &      , trim( HYCOM_FILE )     &
     &      , Im, Jm, Nz, 1          &
     &      , salinity               &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_u'              &
     &      , trim( HYCOM_FILE )     &
     &      , Im, Jm, Nz, 1          &
     &      , water_u                &
     &      , start4D, count4D       &
     &      )
          
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_v'              &
     &      , trim( HYCOM_FILE )     &
     &      , Im, Jm, Nz, 1          &
     &      , water_v                &
     &      , start4D, count4D       &
     &      )
     
          st = st+1
        
        end do

        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE(iNC)
        
        write(*,*) "******************************************************************"

      end do
      
!==== LOOP END ========================================================

      write(*,*) 'FINISH!!'

!**** End of Main program **********************************************
        
    CONTAINS
    
!**** create HYCOM NetCDF **********************************************

      SUBROUTINE createNetCDF(   &
!        input parameters
     &      OUT_FILE             &
     &    , TIME_ATT             &  
     &    , Im, Jm, Nz, Nt       &   
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Im, Jm, Nz, Nt
      
      integer :: ncid2,var_id2
      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: dim3Dids(3), dim4Dids(4)
      
!---- Create the extracted HYCOM netCDF file --------------------------------

      write(*,*) "CREATE: ", trim( OUT_FILE )

      call check( nf90_create(trim( OUT_FILE ), nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'lat', Jm, lat_dimid) )
      call check( nf90_def_dim(ncid2, 'lon', Im, lon_dimid) )
      call check( nf90_def_dim(ncid2, 'depth',Nz, depth_dimid) )
      call check( nf90_def_dim(ncid2, 'time', NF90_UNLIMITED, time_dimid) )

      dim3Dids = (/ lon_dimid, lat_dimid, time_dimid /)
      dim4Dids = (/ lon_dimid, lat_dimid, depth_dimid, time_dimid /)
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid2, 'time', NF90_DOUBLE, time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Valid Time') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      call check( nf90_def_var(ncid2, 'depth', NF90_DOUBLE, depth_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Depth') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter' ) )

      call check( nf90_def_var(ncid2, 'lon', NF90_DOUBLE, lon_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Longitude') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'degrees_east') )

      call check( nf90_def_var(ncid2, 'lat', NF90_DOUBLE, lat_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Latitude') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'degrees_north' ) )

      call check( nf90_def_var(ncid2, 'surf_el', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Water Surface Elevation') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'water_temp', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Water Temperature') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'degC') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'salinity', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Salinity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'psu') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'water_u', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Eastward Water Velocity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'm/s') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

      call check( nf90_def_var(ncid2, 'water_v', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Northward Water Velocity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'm/s') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )

  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDF
      
! -------------------------------------------------------------------------
      
    END PROGRAM getHYCOM
      
