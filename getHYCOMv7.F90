
!!!=== Copyright (c) 2014-2018 Takashi NAKAMURA  =====

#define GOFS_31
/*#define GOFS_30_3H*/
/*#define GOFS_30*/

#define ANALYSIS
/*#define REANALYSIS*/

    PROGRAM getHYCOM
      use netcdf
      use mod_roms_netcdf
      use mod_calendar
     
      implicit none
      
! SETTINGS of HYCOM Data reading area  ------------------------------------
! -------------------------------------------------------------------------
! NetCDF file for the box corners bounded by (Llon,Blat) and (Rlon,Tlat).
!
!                 ______ (Rlon,Tlat)
!                |      |
!                |      |
!                |______|
!     (Llon,Blat)                     
!
! On Input:
!
!    Llon         Box left-edge   longitude (degrees, 0 - 360)
!    Rlon         Box right-edge  longitude (degrees, 0 - 360)
!    Blat         Box bottom-edge latitude  (degress, -90 - 90 )
!    Tlat         Box top-edge    latitude  (degress, -90 - 90 )
!
! Geographical and Grid parameters --------

      real(8), parameter :: Tlat = 28.0d0     ! Latitude  (degrees) of the top-right corner of the grid.
      real(8), parameter :: Blat = 22.0d0     ! Latitude  (degrees) of the bottom-left corner of the grid.
      real(8), parameter :: Llon = 120.0d0    ! Longitude (degrees) of the bottom-left corner of the grid. 
      real(8), parameter :: Rlon = 128.0d0    ! Longitude (degrees) of the top-right corner of the grid. 
!      real(8), parameter :: Tlat =  32.0d0    ! Latitude  (degrees) of the top-right corner of the grid.
!      real(8), parameter :: Blat = -16.0d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
!      real(8), parameter :: Llon =  89.0d0    ! Longitude (degrees) of the bottom-left corner of the grid. 
!      real(8), parameter :: Rlon = 158.0d0    ! Longitude (degrees) of the top-right corner of the grid. 

      integer, parameter :: Syear  = 2018   ! Starting year
      integer, parameter :: Smonth = 4      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
      
      integer, parameter :: Eyear  = 2018   ! Ending year   !!!4/30, 5/31, 6/30, 7/31, 8/31
      integer, parameter :: Emonth = 5      ! Ending month
      integer, parameter :: Eday   = 1      ! Ending day
      
      integer, parameter :: Ryear  = 2000   ! Reference year
      integer, parameter :: Rmonth = 1      ! Reference month
      integer, parameter :: Rday   = 1      ! Reference day
      
     ! NetCDF file     
      character(len=*), parameter :: OUT_prefix  = &
     &  "D:/ROMS/Data/Yaeyama/HYCOM_extracted_data/Yaeyama1_HYCOM_extracted"
!      character(len=*), parameter :: OUT_prefix  = &
!     &  "D:/ROMS/Data/Coral_Triangle/HYCOM_extracted_data/CT_HYCOM_extracted_1604.nc"

#if defined GOFS_31
# if defined ANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2014-July to present) -----
      integer, parameter :: NCnum   = 6
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.2"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.8"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0"  &
     &  /)
# elif defined REANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Renalysis (coming soon) -----
      integer, parameter :: NCnum   = 0
      character(53) :: NC_FILE(NCnum) = (/                          &
     &  /)
# endif
#elif defined GOFS_30_3H
# if defined ANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to present) -----
      integer, parameter :: NCnum   = 0
      character(53) :: NC_FILE(NCnum) = (/                          &
     &  /)
# elif defined REANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31) -----
      integer, parameter :: NCnum   = 2
      character(59) :: NC_FILE(NCnum) = (/                          &  !!! Reanalysis�Ƃ�Lon�̔z�u���Ⴄ�悤���B
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0/3hrly"  &  !!! �v����
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/3hrly"  &
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
      character(53) :: NC_FILE(NCnum) = (/                          &  !!! Reanalysis�Ƃ�Lon�̔z�u���Ⴄ�悤���B
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0"  &  !!! �v����
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"  &
     &  /)
# endif
#endif
! -------------------------------------------------------------------------

      character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
      character(10) :: OUT_suffix    = "_200001.nc"
      character(256) :: OUT_FILE
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation
      
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
      real(8) :: d_jdate_20000101, d_jdate     
      integer :: ncid,var_id
      integer :: ncid2,var_id2
      integer :: N_xi_rho, N_eta_rho
      real(8) :: sf, off
      integer :: Im, Jm, Nz
      integer :: IL, IR, JB, JT
      integer :: iNC, itime
      integer :: end_flag
      integer :: create_flag = 0
      integer :: st = 1
      real(8) :: time(1)
      
!      integer :: time_varid, Uwind_varid
      
      
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
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      
!---- Create the Output file --------------------------------

      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth

      OUT_suffix(2:5)=YYYY
      OUT_suffix(6:7)=MM
      OUT_FILE = OUT_prefix//OUT_suffix

!---- Allocate lat, lon dimension data ---------------------------------
      
      write(*,*) "OPEN: ", NC_FILE(NCnum)
      call check( nf90_open(NC_FILE(NCnum), nf90_nowrite, ncid) )
      
      ! Get dimension data
      call get_dimension(ncid, 'lat', Jm)
      call get_dimension(ncid, 'lon', Im)
      call get_dimension(ncid, 'depth', Nz)
      
!      ! Allocate variable
      allocate(lat_all(Jm))
      allocate(lon_all(Im))
      allocate(depth(Nz))
          
      call check( nf90_inq_varid(ncid, 'lat', var_id) )
      call check( nf90_get_var(ncid, var_id, lat_all) )
      call check( nf90_inq_varid(ncid, 'lon', var_id) )
      call check( nf90_get_var(ncid, var_id, lon_all) )
      call check( nf90_inq_varid(ncid, 'depth', var_id) )
      call check( nf90_get_var(ncid, var_id, depth) )
      call check( nf90_close(ncid) )         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST code
      write(*,*) "CLOSE: ", NC_FILE(NCnum)     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST code

      JB=1
      do j=1,Jm
        if(lat_all(j)>=Blat) exit
        JB=j
      end do
      JT=JB
      do j=JB,Jm
        if(lat_all(j)>Tlat) exit
        JT=j+1
      end do
      IL=1
      do i=IL,Im
        if(lon_all(i)>=Llon) exit
        IL=i
      end do
      IR=IL
      do i=IL,Im
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
     &    trim( OUT_FILE )       &
     &  , TIME_ATT               &  
     &  , Im, Jm, Nz, N_days+1   &   
     &  )
      start1D = (/ 1  /)
      count1D = (/ Im /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lon'                  &
     &  , trim( OUT_FILE )       &
     &  , Im                     &
     &  , lon                    &
     &  , start1D, count1D       &
     &  )
      start1D = (/ 1  /)
      count1D = (/ Jm /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lat'                  &
     &  , trim( OUT_FILE )       &
     &  , Jm                     &
     &  , lat                    &
     &  , start1D, count1D       &
     &  )

      start1D = (/ 1  /)
      count1D = (/ Nz /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'depth'                &
     &  , trim( OUT_FILE )       &
     &  , Nz                     &
     &  , depth                  &
     &  , start1D, count1D       &
     &  )
     
      
!---- Read HYCOM netCDF file --------------------------------
      write(*,*) "******************************************************************"
      do iNC=1, NCnum
        ! Open NetCDF file
        write(*,*) "OPEN: ", NC_FILE(iNC)
      write(*,*) 'CHECK: Time'
        call check( nf90_open(NC_FILE(iNC), nf90_nowrite, ncid) )
        call get_dimension(ncid, 'time', NC(iNC)%Nt)
        allocate(NC(iNC)%time_all(NC(iNC)%Nt))
        call check( nf90_inq_varid(ncid, 'time', var_id) )
        call check( nf90_get_var(ncid, var_id, NC(iNC)%time_all) )
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE(iNC)
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
        do i=NC(iNC)%Nt-1,1,-1
          d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24
          if(d_jdate < dble(jdate_End)) then
            write(*,*) '*** FOUND: Ending point @ NC_FILE',iNC
            NC(iNC)%ItE=i+1
            exit
          endif
        end do
      end do
      write(*,*) NC(:)%ItE 
      
      do iNC=1,NCnum
        do i=1,NC(iNC)%ItE
          d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24
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
        
        do itime=NC(iNC)%ItS, NC(iNC)%ItE
        
          time(1) = NC(iNC)%time_all(itime)/24.0d0 + dble(jdate_20000101-jdate_Ref)
        
          write(*,*) "******************************************************************"
          write(*,*) 'time = ', time(1), TIME_ATT
          
! --- Read NetCDF file ------------------------
          
          start3D = (/ IL, JB, itime /)
          count3D = (/ Im, Jm, 1  /)
          
          call readNetCDF_3d(        &
!          input parameters
     &        NC_FILE(iNC)           &
     &      , 'surf_el'              &
     &      , Im, Jm, 1              &
     &      , start3D, count3D       &
!          output parameters
     &      , surf_el                &
     &)
     
          start4D = (/ IL, JB, 1,  itime /)
          count4D = (/ Im, Jm, Nz, 1  /)

          call readNetCDF_4d(        &
!          input parameters
     &        NC_FILE(iNC)           &
     &      , 'water_temp'           &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , water_temp             &
     &)
          call readNetCDF_4d(        &
!          input parameters
     &        NC_FILE(iNC)           &
     &      , 'salinity'             &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , salinity               &
     &)
          call readNetCDF_4d(        &
!          input parameters
     &        NC_FILE(iNC)           &
     &      , 'water_u'              &
     &      , Im, Jm, Nz, 1          &
     &      , start4D, count4D       &
!          output parameters
     &      , water_u                &
     &)
          call readNetCDF_4d(        &
!          input parameters
     &        NC_FILE(iNC)           &
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
     &      , trim( OUT_FILE )       &
     &      , 1                      &
     &      , time(1)                &
     &      , start1D, count1D       &
     &      )
          start3D = (/ 1,  1,  st /)
          count3D = (/ Im, Jm, 1 /)
          
          call writeNetCDF_3d(       &
!          input parameters
     &        'surf_el'              &
     &      , trim( OUT_FILE )       &
     &      , Im, Jm, 1              &
     &      , surf_el                &
     &      , start3D, count3D       &
     &      )
          
          start4D = (/ 1,  1,  1,  st /)
          count4D = (/ Im, Jm, Nz, 1 /)
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_temp'           &
     &      , trim( OUT_FILE )       &
     &      , Im, Jm, Nz, 1          &
     &      , water_temp             &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'salinity'             &
     &      , trim( OUT_FILE )       &
     &      , Im, Jm, Nz, 1          &
     &      , salinity               &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_u'              &
     &      , trim( OUT_FILE )       &
     &      , Im, Jm, Nz, 1          &
     &      , water_u                &
     &      , start4D, count4D       &
     &      )
          
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_v'              &
     &      , trim( OUT_FILE )       &
     &      , Im, Jm, Nz, 1          &
     &      , water_v                &
     &      , start4D, count4D       &
     &      )
     
          st = st+1
        
        end do
        
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
      
