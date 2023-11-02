
!!!=== Copyright (c) 2023 Takashi NAKAMURA  =====

PROGRAM windROMS2SWAN
  use netcdf
  use mod_roms_netcdf
  use mod_interpolation
  use mod_calendar
 
  implicit none
!--SWAN parameters ------------------------------------------------------ 
  real(8) :: xpc, ypc      ! geographic location of the origin of the computational grid (m)
  real(8) :: alpc          ! direction of the positive x−axis of the computational grid (degree)
  real(8) :: xlenc, ylenc  ! length of the computational grid (m/degree)
  integer :: mxc, myc      ! number of meshes in computational grid
                           !   *this number is one less than the number of grid points in this domain!
                           !   *Maybe because indices start from 0???
  real(8) :: dxinp, dyinp  ! mesh size of the input grid (m/degree)
!-------------------------------------------------------------------------------
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  character(256) :: GRID_FILE
  integer :: Ryear, Rmonth, Rday
  integer :: NCnum
  character(256), allocatable :: ATM_FILE(:)
  character(256) :: OUTPUT_Dir
  real(8) :: grid_size, angle

  TYPE T_NC
    real(8), pointer :: time_all(:)
    integer :: Nt
    integer :: ItS, ItE
  END TYPE T_NC
  TYPE (T_NC), allocatable :: NC(:)
  real(8), allocatable :: time2(:)
  integer, allocatable :: iNCt(:)
  integer, allocatable :: idt(:)
  integer :: iNCs, iNCe
  real(8) :: d_jdate
  real(8) :: d_jdate_Start, d_jdate_Ref
  integer :: jdate_Start, jdate_End, jdate_Ref
  integer :: jdate
  integer :: N_days
  integer :: iNC, iNCm
  integer :: Nt
  integer :: iHR

  character(256) :: NC_NAME(2) = (/ &
     "Uwind"                        &
    ,"Vwind"                        &
    /)

  real(8), parameter :: PI = 3.141592653589793d0

  real(8), allocatable :: yr(:, :)
  real(8), allocatable :: xr(:, :)
  real(8), allocatable :: angler(:, :)
#if defined   OUTPUT_SWAN_GRID  
  real(8), allocatable :: h(:, :)
  real(8), allocatable :: rmask(:, :)
#endif
  real(8), allocatable :: wind(:,:,:)
  real(8), allocatable :: wind2(:,:,:)
       
  integer :: Nxr, Nyr
  integer :: L, M  

  real(8) :: t
  real(8) :: time(1)
  
  character(256) :: IN_FILE
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: i,j,k

  character(15) :: YYYYMMDDphhmmss, YYYYMMDDphhmmss0
  character(8)  :: YYYYMMDD

  character(256) :: WIND_FILE, SERIES_FILE2
  character(len=*), parameter :: SERIES_FILE = "wind_series.dat"
  character(len=*), parameter :: BOTTOM_FILE = 'swan_Bottom.bot'
  character(len=*), parameter :: SWN_GRD_FILE = 'swan_grid.grd'
  character(len=*), parameter :: BLOCK_WIND_TXT = 'swan_in_wind.txt'
  character(len=*), parameter :: BLOCK_BOTTOM_TXT = 'swan_in_bottom.txt'
  
  integer :: ncid,var_id
  integer :: start3D(3), count3D(3)
  
  integer :: iparam,ifc
  logical :: step1
!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/roms2swan_1/NCnum
  namelist/roms2swan_2/ATM_FILE
  namelist/roms2swan_2/OUTPUT_Dir
  namelist/roms2swan_3/grid_size, angle


  ! Read parameters in namelist file
  
  read (5, nml=grd)
  rewind(5)
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
  rewind(5)
  read (5, nml=refdate)
  rewind(5)
  read (5, nml=roms2swan_1)
  allocate( ATM_FILE(NCnum) )
  allocate( NC(NCnum) )
  rewind(5)
  read (5, nml=roms2swan_2)
  rewind(5)
  read (5, nml=roms2swan_3)

     
!---- Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
  L = Nxr-1
  M = Nyr-1
  allocate( yr(0:L, 0:M) )
  allocate( xr(0:L, 0:M) )
#if defined UTM_COORD
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'x_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, xr) )
  call check( nf90_inq_varid(ncid, 'y_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, yr) )
#else
 call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, yr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, xr) )
#endif  
#if defined   OUTPUT_SWAN_GRID
  allocate( h(0:L, 0:M) )
  allocate( rmask(0:L, 0:M) )
  call check( nf90_inq_varid(ncid, 'h', var_id) )
  call check( nf90_get_var(ncid, var_id, h) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, rmask) )
#endif
  allocate( angler(0:L, 0:M) )
  call check( nf90_inq_varid(ncid, 'angle', var_id) )
  call check( nf90_get_var(ncid, var_id, angler) )

  ! Close NetCDF file
  call check( nf90_close(ncid) )

!---- Set SWAN grid coordinate ---------------------------------

  mxc = L
  myc = M
  alpc = angle * 180.0d0/PI
  dxinp = grid_size
  dyinp = grid_size
  xlenc = grid_size*dble(mxc)
  ylenc = grid_size*dble(myc)

  xpc = xr(0,0)
  ypc = yr(0,0)

  allocate( wind(0:mxc, 0:myc, 2) )
  allocate( wind2(0:mxc, 0:myc, 2) )

#if defined   OUTPUT_SWAN_GRID
!==== Output SWAN grid file ===============================
! output BOTTOM file
  do j=0, M
    do i=0, L
      if(rmask(i,j)==0.0d0) then
        h(i,j) = 9999.0d0
      endif 
    enddo
  enddo
  
  open(10,file = BOTTOM_FILE)
  do j=0,myc
    write(10,*) h(0:mxc,j)
  enddo
  close(10)

! output GRID file
  open(10,file = SWN_GRD_FILE)
  do j=0,M
    do i=0,L
      write(10,*) xr(i,j)
    enddo
  enddo
  do j=0,M
    do i=0,L
      write(10,*) yr(i,j)
    enddo
  enddo
  close(10)

  open(10,file=BLOCK_BOTTOM_TXT)
  write(10, '( "&& KEYWORDS TO CREATE AND READ COMPUTATIONAL GRID &&" )' )
  write(10, '( "CGRID REGULAR", 1x,f0.8, 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,i0, 1x,i0, " &" )' ) &
                                    xpc,     ypc,    alpc,   xlenc,   ylenc,   mxc,   myc
  write(10, '( "        CIRCLE 36 0.04 1.0 20",/ )' )

  write(10, '( "&& KEYWORDS TO CREATE AND READ BATHYMETRY GRID &&" )' )
  write(10, '( "INPGRID BOTTOM REGULAR", 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,i0, 1x,i0, 1x,f0.6, 1x,f0.6, " EXC 9.999000e+003" )' ) &
                                             xpc,     ypc,    alpc,   mxc,   myc,   dxinp,   dyinp
  write(10, '( "READINP BOTTOM  1 ", A, " 4 0 FREE" )' )  "'"//trim(BOTTOM_FILE)//"'"
  close(10) 
#endif  

!==== Read ATM file coordinates ==================================

!---- Read ROMS NetCDF atm file --------------------------------
  IN_FILE = trim(ATM_FILE(1))

!======= Set starting time and Ending time ========================

  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)

  jdate_End = jdate_Start + N_days
  
  write(*,*) jdate_Start, jdate_End, N_days

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref

! Check time
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    ! Open NetCDF file
    write(*,*) "OPEN: ", trim( ATM_FILE(iNC) )
    call check( nf90_open(trim( ATM_FILE(iNC) ), nf90_nowrite, ncid) )
    call get_dimension(ncid, 'time', NC(iNC)%Nt)
    write(*,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    allocate( time2(NC(iNC)%Nt) )
    call check( nf90_inq_varid(ncid, 'time', var_id) )
    call check( nf90_get_var(ncid, var_id, time2) )
    call check( nf90_close(ncid) )

    NC(iNC)%time_all = time2 + d_jdate_Ref ! julian date (day)

    deallocate(time2)
  end do

  iHR = int((NC(1)%time_all(2)-NC(1)%time_all(1)+1.0d-4)*24.0d0)

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

!  allocate(atm_time(Nt))
  allocate(idt(Nt))
  allocate(iNCt(Nt))

  i=1
  do iNC=iNCs,iNCe
    do k=0,NC(iNC)%ItE-NC(iNC)%ItS
      iNCt(i+k) = iNC
      idt(i+k)  = NC(iNC)%ItS + k
    enddo
    j=i+NC(iNC)%ItE-NC(iNC)%ItS
!    atm_time(i:j) = NC(iNC)%time_all( NC(iNC)%ItS : NC(iNC)%ItE )  ! julian date (day)
    i=j+1
  enddo

  write(*,*) "*************************************"
  
!==== LOOP set up ==========================================
  SERIES_FILE2 = trim(OUTPUT_Dir)//"/"//SERIES_FILE
  open(11, file=trim(SERIES_FILE2))

  step1=.true.

! LOOP1 START ================================================
  DO iNC=1,NCnum   
    IN_FILE = trim(ATM_FILE(iNC))
!--- LOOP2 START ---------------------------------------------------
    DO ifc=NC(iNC)%ItS,NC(iNC)%ItE
      t = NC(iNC)%time_all(ifc)
      jdate = int(t+1.0d-4)
      call cdate(jdate, iyear, imonth, iday)
      ihour = int( (t+1.0d-4 - dble(jdate)) * 24.0d0 )
      write (YYYYMMDD, '(I4.4, I2.2, I2.2)') iyear, imonth, iday
      write (YYYYMMDDphhmmss, '(I4.4, I2.2, I2.2, ".", I2.2, "0000")') iyear, imonth, iday, ihour
  
      if(ihour == 0) then
        WIND_FILE = trim(OUTPUT_Dir)//"/wind_"//YYYYMMDD//".dat"
        if(step1) then
          YYYYMMDDphhmmss0 = YYYYMMDDphhmmss
          step1=.false.
        else
          close(10)
        endif
        write(*,*) 'OPEN: '//trim(WIND_FILE)
        open(10, file=WIND_FILE)
        write(11,'( A )') trim(WIND_FILE)
      endif
  
      write (10, '( "Time = ", A)') YYYYMMDDphhmmss
      write (*,  '( "Time = ", A)') YYYYMMDDphhmmss
      write(*,*) "OPEN: ", trim( IN_FILE )
      !Open NetCDF file
      call check( nf90_open(trim( IN_FILE ), nf90_nowrite, ncid) )
      start3D = (/   1,   1, ifc /)
      count3D = (/ Nxr, Nyr,   1 /)

      DO iparam=1,2
        ! Get variable id
        write(*,*) "READ: ", trim( NC_NAME(iparam) )
        call check( nf90_inq_varid(ncid, trim( NC_NAME(iparam) ), var_id) ) 
        call check( nf90_get_var(ncid, var_id, wind(:,:,iparam), start=start3D, count=count3D) )
      END DO
      
      call check( nf90_close(ncid) ) 
      
#if defined NAUTICAL
      do j=0,myc
        do i=0,mxc
          wind2(i,j,1) = wind(i,j,1)*cos(angler(i,j)) - wind(i,j,2)*sin(angler(i,j)) 
          wind2(i,j,2) = wind(i,j,1)*sin(angler(i,j)) + wind(i,j,2)*cos(angler(i,j)) 
        enddo
      enddo
#else
      wind2(:,:,:) = wind(:,:,:)
#endif

      ! Write wind data
      write(*,*) "WRITE: Wind data to: "//trim( WIND_FILE )
      DO iparam=1,2
        write (10, '(A)') trim( NC_NAME(iparam) )
        do j=0,myc
          write(10,*) wind2(0:mxc,j,iparam)
        enddo
      END DO

    END DO
!-- LOOP2 END --------------------------------
  END DO
  close(10)
! LOOP1 END --------------------------------
  close(11) 

  open(10,file=BLOCK_WIND_TXT)
  write(10, '( "&& KEYWORD TO CREATE WIND GRID &&" )' )
  write(10, '( "INPGRID WIND REGULAR", 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,i0, 1x,i0, 1x,f0.6, 1x,f0.6, " EXC 9.999000e+003  &" )' ) &
                                           xpc,     ypc,    alpc,   mxc,   myc,   dxinp,   dyinp
  write(10, '( "        NONSTAT ", A, 1x,i0, " HR ", A )' ) YYYYMMDDphhmmss0, iHR, YYYYMMDDphhmmss
  write(10, '( "READINP WIND 1 SERIES ", A, " 4 0 1 1 FREE" )' ) "'"//trim(SERIES_FILE2)//"'"
  close(10) 
  
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM windROMS2SWAN
      
