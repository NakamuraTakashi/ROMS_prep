
!!!=== ver 2014/12/19   Copyright (c) 2014 Takashi NAKAMURA  =====

    PROGRAM getHYCOM
      use netcdf
      use mod_interpolation
      use mod_calendar
     
      implicit none
      
! SETTINGS of HYCOM Data reading area  -----------------------------------------------------------
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
!    Blat         Box bottom-edge latitude  (degress,  -90 - 90 )
!    Tlat         Box top-edge    latitude  (degress,  -90 - 90 )
!
! Geographical and Grid parameters --------

      real(8), parameter :: Tlat = 28.0d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
      real(8), parameter :: Blat = 22.0d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
      real(8), parameter :: Llon = 120.0d0   ! Longitude (degrees) of the bottom-left corner of the grid. 
      real(8), parameter :: Rlon = 128.0d0   ! Longitude (degrees) of the bottom-left corner of the grid. 

      integer, parameter :: Syear  = 2014   ! Start year
      integer, parameter :: Smonth = 4      ! Start month
      integer, parameter :: Sday   = 1      ! Start day
      
      integer, parameter :: Eyear  = 2014   ! End year
      integer, parameter :: Emonth = 10     ! End month
      integer, parameter :: Eday   = 1      ! End day
      
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "ROMS_grid_Yaeyama/Yaeyama1_grd.nc"
      character(len=*), parameter :: OCEAN_FILE  = "Yaeyama1_HYCOM140401_141001.nc"
!      character(len=*), parameter :: BRY_FILE  = "Yaeyama1_140401_141001_bry.nc"
!      character(len=*), parameter :: INI_FILE  = "Yaeyama1_140401_ini.nc"
      
      integer, parameter :: NCnum   = 3
      character(53) :: NC_FILE(NCnum) = (/                          &
     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9"  &
     &   , "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0"  &
     &   , "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1"  &
     &  /)
!      character(53) :: NC_FILE(NCnum) = (/                          &
!     &     "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0"  &
!     &   , "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"  &
!     &  /)
     
      character(30) :: TIME_ATT  = "days since 1992-01-01 00:00:00"
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation

      real(8), allocatable :: time_all(:), lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:), depth(:),time(:),tau(:)
      real(8), allocatable :: time2(:)
      real(8), allocatable :: water_u(:, :, :, :), water_v(:, :, :, :)
      real(8), allocatable :: salinity(:, :, :, :),water_temp(:, :, :, :)
      real(8), allocatable :: surf_el(:, :, :)
      integer :: start1D(1), count1D(1)
      integer :: start3D(3), count3D(3)
      integer :: start4D(4), count4D(4)
      
      integer :: N_days
      integer :: jdate_Start, jdate_End
      integer :: jdate_20000101
      integer :: jdate_Start2
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
      integer :: Im, Jm, Nz, Nt
      integer :: IL, IR, JB, JT
      integer :: ItS, ItE
      integer :: iNC
      integer :: end_flag
      integer :: create_flag = 0
      integer :: st = 1
      
!      integer :: time_varid, Uwind_varid
      
      
      call jd(2000, 1, 1, jdate_20000101)

      write(*,*) jdate_20000101

      call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
      call jd(Syear, Smonth, Sday, jdate_Start)
      jdate_Start2 = jdate_Start
      jdate_End = jdate_Start + N_days
      
      write(*,*) jdate_Start, jdate_End, N_days
      
!      call cdate(julian_date, iyear, imonth, iday)
      
      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth
      write (DD, "(I2.2)") Sday
      
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      

      
      do iNC=1, NCnum
      
        write(*,*) "OPEN: ", NC_FILE(iNC)
              
!---- Read HYCOM netCDF file --------------------------------

        !Open NetCDF file
        call check( nf90_open(NC_FILE(iNC), nf90_nowrite, ncid) )
        call get_dimension(ncid, 'time', Nt)
        allocate(time_all(Nt))
        call check( nf90_inq_varid(ncid, 'time', var_id) )
        call check( nf90_get_var(ncid, var_id, time_all) )
        
        write(*,*) 'CHECK: Time interval'

        ItS=0
        do i=1,Nt
          jdate=jdate_20000101+nint(time_all(i))/24
          if(jdate>=jdate_Start) then
            write(*,*) '*** FOUND: Starting point'
            write(*,*) jdate_Start, jdate
            ItS=i
            exit
          endif
        end do
        if(ItS==0) then
          deallocate(time_all)
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE(iNC)
          write(*,*) "*********************************************"
          cycle
        end if
        
        end_flag = 1
        do i=ItS,Nt
          jdate=jdate_20000101+nint(time_all(i))/24
          if(jdate>=jdate_End) then
            write(*,*) '*** FOUND: Ending point'
            write(*,*) jdate_End, jdate
            ItE=i
            end_flag = 0
            exit
          endif
          ItE=i
          jdate_Start = jdate+1
        end do
        write(*,*) jdate_End, jdate
        write(*,*) Nt, ItE, ItS
        
        Nt=ItE-ItS+1
        write(*,*) Nt, ItE, ItS

        allocate(time(Nt))
        
        time(:) = time_all(ItS:ItE)/24.0d0 + dble(jdate_20000101-jdate_Start2)
        
        ! Get dimension data
        call get_dimension(ncid, 'lat', Jm)
        call get_dimension(ncid, 'lon', Im)
        call get_dimension(ncid, 'depth', Nz)
        
!        ! Allocate variable
        allocate(lat_all(Jm))
        allocate(lon_all(Im))
        allocate(depth(Nz))
            
        call check( nf90_inq_varid(ncid, 'lat', var_id) )
        call check( nf90_get_var(ncid, var_id, lat_all) )
        call check( nf90_inq_varid(ncid, 'lon', var_id) )
        call check( nf90_get_var(ncid, var_id, lon_all) )
        call check( nf90_inq_varid(ncid, 'depth', var_id) )
        call check( nf90_get_var(ncid, var_id, depth) )

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

        write(*,*) Im,Jm,Nz,Nt
        

        if(create_flag == 0) then
        
          create_flag = 1
          
     !---- Create the trimed HYCOM netCDF file --------------------------------
          
          call createNetCDF(         &
!          input parameters
     &        OCEAN_FILE             &
     &      , TIME_ATT               &  
     &      , Im, Jm, Nz, N_days+1   &   
     &      )
          start1D = (/ 1  /)
          count1D = (/ Im /)
          call writeNetCDF_1d(       &
!          input parameters
     &        'lon'                  &
     &      , OCEAN_FILE             &
     &      , Im                     &
     &      , lon                    &
     &      , start1D, count1D       &
     &      )
          start1D = (/ 1  /)
          count1D = (/ Jm /)
          call writeNetCDF_1d(       &
!          input parameters
     &        'lat'                  &
     &      , OCEAN_FILE             &
     &      , Jm                     &
     &      , lat                    &
     &      , start1D, count1D       &
     &      )

          start1D = (/ 1  /)
          count1D = (/ Nz /)
          call writeNetCDF_1d(       &
!          input parameters
     &        'depth'                &
     &      , OCEAN_FILE             &
     &      , Nz                     &
     &      , depth                  &
     &      , start1D, count1D       &
     &      )
     
        end if

        ! Initialize count and start for reading data
        start3D = (/ IL, JB, ItS /)
        count3D = (/ Im, Jm, Nt  /)
        start4D = (/ IL, JB, 1,  ItS /)
        count4D = (/ Im, Jm, Nz, Nt  /)
        
        allocate(surf_el(Im,Jm,Nt))
        allocate(water_temp(Im,Jm,Nz,Nt))
        allocate(salinity(Im,Jm,Nz,Nt))
        allocate(water_u(Im,Jm,Nz,Nt))
        allocate(water_v(Im,Jm,Nz,Nt))

        write(*,*) 'DOWNLOAD surf_el'
!     Get variable id
        call check( nf90_inq_varid(ncid, 'surf_el', var_id) ) ! Water Surface Elevation (m)
        call check( nf90_get_var(ncid, var_id, surf_el, start=start3D, count=count3D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        surf_el(:,:,:)=surf_el(:,:,:)*sf+off
        write(*,*) '*** SUCCESS'
        
        write(*,*) 'DOWNLOAD water_temp'
!     Get variable id
        call check( nf90_inq_varid(ncid, 'water_temp', var_id) ) ! Water Temperature (degC)
        call check( nf90_get_var(ncid, var_id, water_temp, start=start4D, count=count4D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        water_temp(:,:,:,:)=water_temp(:,:,:,:)*sf+off
        write(*,*) '*** SUCCESS'
        
        write(*,*) 'DOWNLOAD salinity'
!     Get variable id
        call check( nf90_inq_varid(ncid, 'salinity', var_id) ) ! Salinity (psu)
        call check( nf90_get_var(ncid, var_id, salinity, start=start4D, count=count4D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        salinity(:,:,:,:)=salinity(:,:,:,:)*sf+off
        write(*,*) '*** SUCCESS'
        
        write(*,*) 'DOWNLOAD water_u'
!     Get variable id
        call check( nf90_inq_varid(ncid, 'water_u', var_id) ) ! Eastward Water Velocity (m/s)
        call check( nf90_get_var(ncid, var_id, water_u, start=start4D, count=count4D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        water_u(:,:,:,:)=water_u(:,:,:,:)*sf+off
        write(*,*) '*** SUCCESS'

        write(*,*) 'DOWNLOAD water_v'
!     Get variable id
        call check( nf90_inq_varid(ncid, 'water_v', var_id) ) ! Northward Water Velocity (m/s)
        call check( nf90_get_var(ncid, var_id, water_v, start=start4D, count=count4D) )
        call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
        call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
        water_v(:,:,:,:)=water_v(:,:,:,:)*sf+off
        write(*,*) '*** SUCCESS'
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE(iNC)
          
!        open(50, file='time.dat')
!        do i=1,Nt
!          write(50,*) time(i)
!        end do
!        close(50)
!        write(50,*) water_u
        
        start1D = (/ st  /)
        count1D = (/ Nt /)
        
        call writeNetCDF_1d(       &
!        input parameters
     &      'time'                 &
     &    , OCEAN_FILE             &
     &    , Nt                     &
     &    , time                   &
     &    , start1D, count1D       &
     &    )
        start3D = (/ 1,  1,  st /)
        count3D = (/ Im, Jm, Nt /)
        
        call writeNetCDF_3d(       &
!        input parameters
     &      'surf_el'              &
     &    , OCEAN_FILE             &
     &    , Im, Jm, Nt             &
     &    , surf_el                &
     &    , start3D, count3D       &
     &    )
        
        start4D = (/ 1,  1,  1,  st /)
        count4D = (/ Im, Jm, Nz, Nt /)
        call writeNetCDF_4d(       &
!        input parameters
     &      'water_temp'           &
     &    , OCEAN_FILE             &
     &    , Im, Jm, Nz, Nt         &
     &    , water_temp             &
     &    , start4D, count4D       &
     &    )
        call writeNetCDF_4d(       &
!        input parameters
     &      'salinity'             &
     &    , OCEAN_FILE             &
     &    , Im, Jm, Nz, Nt         &
     &    , salinity               &
     &    , start4D, count4D       &
     &    )
        call writeNetCDF_4d(       &
!        input parameters
     &      'water_u'              &
     &    , OCEAN_FILE             &
     &    , Im, Jm, Nz, Nt         &
     &    , water_u                &
     &    , start4D, count4D       &
     &    )
        
        call writeNetCDF_4d(       &
!        input parameters
     &      'water_v'              &
     &    , OCEAN_FILE             &
     &    , Im, Jm, Nz, Nt         &
     &    , water_v                &
     &    , start4D, count4D       &
     &    )
     
        st = Nt+1
        
        deallocate(time_all)        
        deallocate(lat_all)
        deallocate(lon_all)
        deallocate(depth)
        deallocate(time)        
        deallocate(lat)
        deallocate(lon)
        deallocate(surf_el)
        deallocate(water_temp)
        deallocate(salinity)
        deallocate(water_u)
        deallocate(water_v)
        
        write(*,*) "*********************************************"
          
        if(end_flag == 0) exit
        
      end do

      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
        
    CONTAINS
    
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
      
!---- Create the HYCOM netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'lat', Jm, lat_dimid) )
      call check( nf90_def_dim(ncid2, 'lon', Im, lon_dimid) )
      call check( nf90_def_dim(ncid2, 'depth',Nz, depth_dimid) )
      call check( nf90_def_dim(ncid2, 'time', NF90_UNLIMITED, time_dimid) )

      dim3Dids = (/ lon_dimid, lat_dimid, time_dimid /)
      dim4Dids = (/ lon_dimid, lat_dimid, depth_dimid, time_dimid /)
      
    ! Define the netCDF variables for the pressure and temperature data.
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
      
      
!---- NetCDF utility -------------------------------------------------
     
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension


      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
    END PROGRAM getHYCOM
      
