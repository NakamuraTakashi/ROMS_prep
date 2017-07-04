
!!!=== Copyright (c) 2014-2017 Takashi NAKAMURA  =====

#define EXPT_9X

    PROGRAM getHYCOM
      use netcdf
      use mod_calendar
     
      implicit none
      
! SETTINGS of HYCOM Data reading area  ------------------------------------
! -------------------------------------------------------------------------
! NetCDF file for the box corners bounded by (Llon,Blat) and (Rlon,Tlat).
!
!  　　　　　　　 ______ (Rlon,Tlat)
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
!      real(8), parameter :: Tlat =  31.0d0    ! Latitude  (degrees) of the top-right corner of the grid.
!      real(8), parameter :: Blat = -18.0d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
!      real(8), parameter :: Llon =  98.0d0    ! Longitude (degrees) of the bottom-left corner of the grid. 
!      real(8), parameter :: Rlon = 140.0d0    ! Longitude (degrees) of the top-right corner of the grid. 

      integer, parameter :: Syear  = 2016   ! Starting year
      integer, parameter :: Smonth = 4      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
      
      integer, parameter :: Eyear  = 2016   ! Ending year   !!!4/30, 5/31, 6/30, 7/31, 8/31
      integer, parameter :: Emonth = 5      ! Ending month
      integer, parameter :: Eday   = 1      ! Ending day
      
      integer, parameter :: Ryear  = 2000   ! Reference year
      integer, parameter :: Rmonth = 1      ! Reference month
      integer, parameter :: Rday   = 1      ! Reference day
      
     ! NetCDF file     
      character(len=*), parameter :: OCEAN_FILE  = &
     &  "D:/ROMS/Data/Yaeyama/HYCOM_extracted_data/Yaeyama1_HYCOM_extracted_1604.nc"
!      character(len=*), parameter :: OCEAN_FILE  = &
!     &  "D:/ROMS/Data/Coral_Triangle/HYCOM_extracted_data/CT_HYCOM_extracted_1604.nc"

#ifdef EXPT_9X
! ----- HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to present) -----
      integer, parameter :: NCnum   = 4
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1"  &
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2"  &
     &  /)
#else
! ----- HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31) -----
      integer, parameter :: NCnum   = 2
      character(53) :: NC_FILE(NCnum) = (/                          &  !!! ReanalysisとはLonの配置が違うようだ。
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0"  &  !!! 要検討
     &    ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"  &
     &  /)
#endif
! -------------------------------------------------------------------------

      character(30) :: TIME_ATT  = "days since 1992-01-01 00:00:00"
     
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
      integer :: jdate
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      
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

      write(*,*) jdate_20000101

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
     &    OCEAN_FILE             &
     &  , TIME_ATT               &  
     &  , Im, Jm, Nz, N_days+1   &   
     &  )
      start1D = (/ 1  /)
      count1D = (/ Im /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lon'                  &
     &  , OCEAN_FILE             &
     &  , Im                     &
     &  , lon                    &
     &  , start1D, count1D       &
     &  )
      start1D = (/ 1  /)
      count1D = (/ Jm /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'lat'                  &
     &  , OCEAN_FILE             &
     &  , Jm                     &
     &  , lat                    &
     &  , start1D, count1D       &
     &  )

      start1D = (/ 1  /)
      count1D = (/ Nz /)
      call writeNetCDF_1d(       &
!      input parameters
     &    'depth'                &
     &  , OCEAN_FILE             &
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
          jdate=jdate_20000101+nint(NC(iNC)%time_all(i))/24
          if(jdate < jdate_End) then
            write(*,*) '*** FOUND: Ending point @ NC_FILE',iNC
            NC(iNC)%ItE=i+1
            exit
          endif
        end do
      end do
      
      do iNC=1,NCnum
        do i=2,NC(iNC)%ItE
          jdate=jdate_20000101+nint(NC(iNC)%time_all(i))/24
          if(jdate>jdate_Start) then
            write(*,*) '*** FOUND: Starting point @ NC_FILE',iNC
            NC(iNC)%ItS=i-1
            exit
          endif
        end do
      end do
      
      write(*,*) NC(:)%ItS, NC(:)%ItE 
        
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
     &      , OCEAN_FILE             &
     &      , 1                      &
     &      , time(1)                &
     &      , start1D, count1D       &
     &      )
          start3D = (/ 1,  1,  st /)
          count3D = (/ Im, Jm, 1 /)
          
          call writeNetCDF_3d(       &
!          input parameters
     &        'surf_el'              &
     &      , OCEAN_FILE             &
     &      , Im, Jm, 1              &
     &      , surf_el                &
     &      , start3D, count3D       &
     &      )
          
          start4D = (/ 1,  1,  1,  st /)
          count4D = (/ Im, Jm, Nz, 1 /)
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_temp'           &
     &      , OCEAN_FILE             &
     &      , Im, Jm, Nz, 1          &
     &      , water_temp             &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'salinity'             &
     &      , OCEAN_FILE             &
     &      , Im, Jm, Nz, 1          &
     &      , salinity               &
     &      , start4D, count4D       &
     &      )
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_u'              &
     &      , OCEAN_FILE             &
     &      , Im, Jm, Nz, 1          &
     &      , water_u                &
     &      , start4D, count4D       &
     &      )
          
          call writeNetCDF_4d(       &
!          input parameters
     &        'water_v'              &
     &      , OCEAN_FILE             &
     &      , Im, Jm, Nz, 1         &
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
    
!**** createNetCDF **********************************************

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

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

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
      
!**** writeNetCDF_1d **********************************************
      
      SUBROUTINE writeNetCDF_1d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im                     &
     &    , data                   &
     &    , start1D, count1D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im 
      real(8), intent( in) :: data(Im )
      integer, intent( in) :: start1D(1), count1D(1)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start1D, count = count1D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_1d
      
!**** writeNetCDF_3d **********************************************
      
      SUBROUTINE writeNetCDF_3d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nt             &
     &    , data                   &
     &    , start3D, count3D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nt 
      real(8), intent( in) :: data(Im, Jm, Nt )
      integer, intent( in) :: start3D(3), count3D(3)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_3d
      
!**** writeNetCDF_4d **********************************************
      
      SUBROUTINE writeNetCDF_4d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nz, Nt         &
     &    , data                   &
     &    , start4D, count4D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nz, Nt
      real(8), intent( in) :: data(Im, Jm, Nz, Nt)
      integer, intent( in) :: start4D(4), count4D(4)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start4D, count = count4D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_4d
      
!**** readNetCDF_3d **********************************************
      
      SUBROUTINE readNetCDF_3d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nt             &
     &    , start3D, count3D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nt
      integer, intent( in) :: start3D(3), count3D(3)
      real(8), intent(out) :: data(Im, Jm, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start3D, count=count3D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do

      data(:,:,:)=data(:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_3d
      
!**** readNetCDF_4d **********************************************
      
      SUBROUTINE readNetCDF_4d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start4D, count4D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start4D(4), count4D(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start4D, count=count4D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do
      
      data(:,:,:,:)=data(:,:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_4d
      
!**** NetCDF utility **********************************************
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension

! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          err_flag = 1
!          stop "Stopped"
      end if
      
      END SUBROUTINE check2
      
! -------------------------------------------------------------------------
      
    END PROGRAM getHYCOM
      
