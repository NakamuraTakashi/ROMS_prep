
!!!=== Copyright (c) 2023 Takashi NAKAMURA  =====

PROGRAM grdROMS2SWAN
  use netcdf
  use mod_roms_netcdf
 
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
  character(256) :: GRID_FILE
  integer :: NCnum
  character(256), allocatable :: ATM_FILE(:)
  character(256) :: OUTPUT_Dir
  real(8) :: grid_size, angle

  real(8), parameter :: PI = 3.141592653589793d0

  real(8), allocatable :: yr(:, :)
  real(8), allocatable :: xr(:, :)
  real(8), allocatable :: angler(:, :)
  real(8), allocatable :: h(:, :)
  real(8), allocatable :: rmask(:, :)
       
  integer :: Nxr, Nyr
  integer :: L, M  

  integer :: i,j,k

  character(len=*), parameter :: BOTTOM_FILE = 'swan_Bottom.bot'
  character(len=*), parameter :: SWN_GRD_FILE = 'swan_grid.grd'
  character(len=*), parameter :: BLOCK_BOTTOM_TXT = 'swan_in_bottom.txt'
  character(len=*), parameter :: BLOCK_NGRID_TXT = 'swan_in_NGRID.txt'
  
  integer :: ncid,var_id
  integer :: start3D(3), count3D(3)
!
!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/roms2swan_3/grid_size, angle


  ! Read parameters in namelist file
  read (5, nml=grd)
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

  allocate( h(0:L, 0:M) )
  allocate( rmask(0:L, 0:M) )
  call check( nf90_inq_varid(ncid, 'h', var_id) )
  call check( nf90_get_var(ncid, var_id, h) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, rmask) )

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

  open(10,file=BLOCK_NGRID_TXT)
  write(10, '( "NGRID ", A, 1x,f0.8, 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,f0.6, 1x,i0, 1x,i0 )' ) &
                  "'sname'",     xpc,     ypc,    alpc,   xlenc,   ylenc,   mxc,   myc
  write(10, '( "NESTOUT  ", A, 1x, A, " OUTPUT YYYYMMDD.000000 1 HR" )' )  &
                    "'sname'", "'fname.nst'"
  close(10) 
  
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM grdROMS2SWAN
      
