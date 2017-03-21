
!!!=== ver 2015/02/04   Copyright (c) 2014-2015 Takashi NAKAMURA  =====

    PROGRAM getJMACWM
      use netcdf
      use mod_interpolation
      use mod_calendar
      use grib_mod
     
      implicit none
      
! SETTINGS of JMA MSM Data reading area  -----------------------------------------------------------
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
!    Llon         Box left-edge   longitude (degrees, -180 - 180)
!    Rlon         Box right-edge  longitude (degrees, -180 - 180)
!    Blat         Box bottom-edge latitude  (degress,  -90 - 90 )
!    Tlat         Box top-edge    latitude  (degress,  -90 - 90 )
!
! Geographical and Grid parameters --------

      real(8), parameter :: Tlat = 27.5d0  ! 24.9d0   !27.5d0    ! Latitude  (degrees) of the bottom-left corner of the grid.
      real(8), parameter :: Blat = 22.5d0  ! 23.7d0   !22.5d0    ! Latitude  (degrees) of the top-right corner of the grid.
      real(8), parameter :: Llon = 120.0d0 ! 123.5d0  !120.0d0   ! Longitude (degrees) of the bottom-left corner of the grid. 
      real(8), parameter :: Rlon = 127.0d0 ! 125.0d0  !127.0d0   ! Longitude (degrees) of the top-right corner of the grid. 

      integer, parameter :: Syear  = 2014   ! Starting year
      integer, parameter :: Smonth = 4      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
      
      integer, parameter :: Eyear  = 2014   ! Ending year
      integer, parameter :: Emonth = 10     ! Ending month
      integer, parameter :: Eday   = 1      ! Ending day
      
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "D:/ROMS/Yaeyama/Data/Yaeyama1_grd_ver3.nc"
      character(len=*), parameter :: OUT_FILE  = "Yaeyama1_wave_JMA_140401_141001.dat"

      character(len=*), parameter :: GRIB2_FILE = &
     & "Data/Z__C_RJTD_20140401000000_CWM_GPV_Rjp_Gll0p05deg_FD0000-0300_grib2.bin"
     
      character(30) :: TIME_ATT  = "days since 1992-01-01 00:00:00"
     
      integer, parameter :: mode = 1           ! mode=1, linear intrtporation

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      
      real(8) :: time2(1) ! Ocean time
      real(8), allocatable :: Uwind(:,:) ! surface u-wind component (meter second-1)
      real(8), allocatable :: Vwind(:,:) ! surface v-wind component (meter second-1)
      real(8), allocatable :: Pair (:,:) ! surface air pressure (milibar=hPa)
      real(8), allocatable :: Tair (:,:) ! surface air temperature (Celsius)
      real(8), allocatable :: Qair (:,:) ! surface air relative humidity (percentage)
      real(8), allocatable :: rain (:,:) ! rain fall rate (kilogram meter-2 second-1)
      real(8), allocatable :: cloud(:,:) ! cloud fraction (0 to 1)    
           
      real(8), allocatable :: lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:), time(:)
      real(8), allocatable :: u(:, :, :), v(:, :, :),temp(:, :, :)
      real(8), allocatable :: psea(:, :, :), rh(:, :, :),ncld(:, :, :)
      real(8), allocatable :: r1h(:, :, :)
      integer :: start1D(1), count1D(1)
      integer :: start3D(3), count3D(3)
      real(8), allocatable :: u2(:, :)
      
      integer :: N_days
      integer :: julian_date
      integer :: iyear, imonth, iday
      integer :: i,j,k
      integer :: idays
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      real(8) :: sf, off
      integer :: Im, Jm, Nt
      integer :: IL, IR, JB, JT
      
      integer :: ncid2,var_id2
      integer :: start1Dw(1), count1Dw(1)
      integer :: start3Dw(3), count3Dw(3)
      integer :: dimids(3)
      integer :: xi_rho_dimid, eta_rho_dimid, time_dimid
!      integer :: time_varid, Uwind_varid
! for test      
      integer :: number
      type(gribfield) :: gfld
      integer,dimension(200) :: jids,jpdt,jgdt
      logical :: unpack=.true.
      ifile=10


      call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
      call jd(Syear, Smonth, Sday, julian_date)
      
!      write(*,*) N_days
      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(Uwind(N_xi_rho, N_eta_rho))
      allocate(Vwind(N_xi_rho, N_eta_rho))
      allocate(Pair (N_xi_rho, N_eta_rho))
      allocate(Tair (N_xi_rho, N_eta_rho))
      allocate(Qair (N_xi_rho, N_eta_rho))
      allocate(rain (N_xi_rho, N_eta_rho))
      allocate(cloud(N_xi_rho, N_eta_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      
      call cdate(julian_date, iyear, imonth, iday)
      
      write (YYYY, "(I4.4)") iyear
      write (MM, "(I2.2)") imonth
      write (DD, "(I2.2)") iday
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      
!---- Update file name  -----------------------------------------------
       
!      NC_FILE(77:80)=YYYY
!      NC_FILE(82:83)=MM
!      NC_FILE(84:85)=DD
!      
!      NC_FILE2(81:84)=YYYY
!      NC_FILE2(86:87)=MM
!      NC_FILE2(88:89)=DD
      
      write(*,*) "OPEN: ", GRIB2_FILE
            
!---- Read JMA CWM GRIB2 file --------------------------------

  ! Open GRIB2 file 
      call baopenr(ifile,GRIB2_FILE,iret)
      .
  ! Set GRIB2 field identification values to search for
      jdisc=
      jids(?)=
      jpdtn=
      jpdt(?)=
      jgdtn=
      jgdt(?)=

  ! Get field from file
      call getgb2(ifile,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     &            unpack,j,gfld,iret)

  ! Process field ...
      firstval=gfld%fld(1)
      lastval=gfld%fld(gfld%ndpts)
      fldmax=maxval(gfld%fld)
      fldmin=minval(gfld%fld)
      
  ! Free memory when done with field
      call gf_free(gfld)



      julian_date = julian_date - 2
      
!---- LOOP START --------------------------------
!
!      do idays=1, N_days+1
!
!! --- Update  -----------------------------------------------
!
!        julian_date = julian_date + 1
!        
!        call cdate(julian_date, iyear, imonth, iday)
!      
!        write (YYYY, "(I4.4)") iyear
!        write (MM, "(I2.2)") imonth
!        write (DD, "(I2.2)") iday
!         
!        NC_FILE(77:80)=YYYY
!        NC_FILE(82:83)=MM
!        NC_FILE(84:85)=DD
!        
!        NC_FILE2(81:84)=YYYY
!        NC_FILE2(86:87)=MM
!        NC_FILE2(88:89)=DD
!      
!! --- Read JMA MSM netCDF file --------------------------------
!      
!        write(*,*) "OPEN: ", NC_FILE
!      
!        !Open NetCDF file
!        
!        
!        
!        write(*,*) "*** SUCCESS reading JMA MSM netcdf file!"
!        
!        do k=1,Nt
!          
!     
!     
!          write(*,*) '*************************************************'
!        end do
!        
!
!      end do
      
      write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
        
    CONTAINS
    
      
      
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
      
    END PROGRAM getJMACWM
      
