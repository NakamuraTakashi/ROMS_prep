
!!!=== Copyright (c) 2024 Takashi NAKAMURA  =====

PROGRAM grdROMS_add_ncparams
  use netcdf
  use mod_interpolation
  use mod_roms_netcdf
 
  implicit none
  
! --namelist grd parameters ----------------------------------------------
  character(256) :: GRID_FILE
! --namelist utm_zone parameters ------------------------------------------------
  integer :: izone, ispher
! -------------------------------------------------------------------------
  
  real(8), allocatable :: h(:,:)        ! depth (meter)
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: xr(:, :)
  real(8), allocatable :: yr(:, :)
  real(8), allocatable :: pm(:, :)
  real(8), allocatable :: pn(:, :)
  real(8), allocatable :: rmask(:, :)
  real(8), allocatable :: area(:, :)
  
  integer :: i,j

  integer :: Nxr, Nyr
  integer :: Nxp, Nyp
  integer :: Nxu, Nyu
  integer :: Nxv, Nyv
  integer :: L, M
 
  integer :: ncid,var_id
  integer :: status
  integer :: start2D(2), count2D(2), dim2Dids(2)
  character(256), allocatable :: varname(:)
  character(2) :: varnum
  integer :: xr_dimid, yr_dimid
  integer :: Nid

  real(8), allocatable :: data(:,:,:)
  integer :: ios

#if defined AQUACULTURE1
! --namelist aquaculture1 parameters -------------------------------------
  character(256) :: INPUT_CSV
  integer :: Nheader
  integer :: Nidaq
  real(8) :: aqunit_interval
  character(256) :: NC_VARNAME, NC_LONGNAME, NC_UNIT
! -------------------------------------------------------------------------
  integer :: Ndata
  character(10) :: str_idaq
  integer, allocatable :: idaq(:)
  real(8), allocatable :: lat_s(:), lat_e(:), lon_s(:), lon_e(:)
  real(8), allocatable :: y_s(:), y_e(:), x_s(:), x_e(:)
  real(8) :: y_au, x_au
  real(8) :: length_aq
  integer :: N_au
  integer :: Iau, Jau
  integer :: i_au
#elif defined AQUACULTURE2
! --namelist aquaculture1 parameters -------------------------------------
  character(256) :: INPUT_CSV
  integer :: Nheader
  character(256) :: NC_VARNAME, NC_LONGNAME, NC_UNIT
! -------------------------------------------------------------------------
  integer :: Ndata
  character(10) :: str_idaq
  integer, allocatable :: idaq(:)
  real(8), allocatable :: lat_c(:), lon_c(:), area_aq(:)
  real(8), allocatable :: y_c(:), x_c(:)
  integer :: Iau, Jau
#endif
#if defined OUTPUT_SWAN_GRID
!--SWAN parameters ------------------------------------------------------ 
  integer :: mxc, myc      ! number of meshes in computational grid
                           !   *this number is one less than the number of grid points in this domain!
                           !   *Maybe because indices start from 0???
!-------------------------------------------------------------------------------
  character(256) :: SWAN_INPGRID_PREFIX = 'aquaculture'
  character(256) :: SWAN_INPGRID_FILE
#endif

  namelist/grd/GRID_FILE
  namelist/utm_zone/izone, ispher
#if defined AQUACULTURE1
  namelist/aquaculture1/INPUT_CSV, Nheader, Nidaq, aqunit_interval
  namelist/aquaculture1/NC_VARNAME, NC_LONGNAME, NC_UNIT
#elif defined AQUACULTURE2
  namelist/aquaculture2/INPUT_CSV, Nheader
  namelist/aquaculture2/NC_VARNAME, NC_LONGNAME, NC_UNIT
#endif

! Read parameters in namelist file
  
  read (5, nml=grd)
  rewind(5)
  read (5, nml=utm_zone)
#if defined AQUACULTURE1
  rewind(5)
  read (5, nml=aquaculture1)
#elif defined AQUACULTURE2
  rewind(5)
  read (5, nml=aquaculture2)
#endif

!-Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )

  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)

  Nxu = Nxr-1
  Nyu = Nyr
  Nxv = Nxr
  Nyv = Nyr-1
  Nxp = Nxr-1
  Nyp = Nyr-1
  
  L = Nxr-1
  M = Nyr-1
  
  allocate( h(0:L, 0:M) )
  allocate( latr(0:L, 0:M) )
  allocate( lonr(1:L, 0:M) )
  allocate( xr(0:L, 0:M) )
  allocate( yr(0:L, 0:M) )
  allocate( pm(0:L, 0:M) )
  allocate( pn(0:L, 0:M) )
  allocate( rmask(0:L, 0:M) )
  allocate( area(0:L, 0:M) )
  
  ! Get variables
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr) )
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, h) )
  call check( nf90_inq_varid(ncid, 'pm', var_id) ) 
  call check( nf90_get_var(ncid, var_id, pm) )
  call check( nf90_inq_varid(ncid, 'pn', var_id) ) 
  call check( nf90_get_var(ncid, var_id, pn) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, rmask) )
#if defined UTM_COORD
  call check( nf90_inq_varid(ncid, 'x_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, xr) )
  call check( nf90_inq_varid(ncid, 'y_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, yr) )
#endif

  ! Close NetCDF file
  call check( nf90_close(ncid) )

  write(*,*) "CLOSE: ", trim( GRID_FILE )

  area(:,:)=1.0d0/pm(:,:)/pn(:,:)

!=========== AQUACULTURE1 ===============================
!---- Compute aquaculture data for ROMS grid netCDF file --------------------------------
#if defined AQUACULTURE1

  allocate( data(0:L,0:M, Nidaq) )
  data(:,:,:) = 0.0d0

  write(*,*) "OPEN: ", trim( INPUT_CSV )
  open(10,file=INPUT_CSV, status='old')

  !--- Count data -----
  Ndata = 0
  do i=1, Nheader
    read(10,*)
  end do
  do
    read(10,*,iostat=ios)
    if(ios==-1) exit
    Ndata = Ndata + 1
  end do
  allocate( idaq(Ndata), lat_s(Ndata), lat_e(Ndata), lon_s(Ndata), lon_e(Ndata) )
  allocate( x_s(Ndata), x_e(Ndata), y_s(Ndata), y_e(Ndata) )
  rewind(10) 
  !--- Read data -----
  do i=1, Nheader
    read(10,*)
  end do
  do i=1, Ndata
    ! read aquaculture1 data
    read(10,*) str_idaq, lat_s(i), lat_e(i), lon_s(i), lon_e(i)
    ! convert characters to integer
    read (str_idaq, *) idaq(i)
    ! convert lat lon coordinate to UTM coordinate
    call ll2utm(lat_s(i),lon_s(i),x_s(i),y_s(i),izone,ispher)
    call ll2utm(lat_e(i),lon_e(i),x_e(i),y_e(i),izone,ispher)

!    write(*,*) idaq(i), lat_s(i), lat_e(i), lon_s(i), lon_e(i) ! Debug
!    write(*,*) idaq(i), x_s(i), x_e(i), y_s(i), y_e(i) ! Debug
  end do
  close(10)
  write(*,*) "CLOSE: ", trim( INPUT_CSV )

  write(*,*) "Processing... "
!$omp parallel
!$omp do private(i,i_au,length_aq,N_au,y_au,x_au,Iau,Jau)
  do i=1, Ndata
    length_aq = sqrt( (x_s(i)-x_e(i))**2.0d0 + (y_s(i)-y_e(i))**2.0d0 )
    N_au = int( length_aq/aqunit_interval )
!    write(*,*) i,length_aq,N_au

    do i_au=1, N_au
      x_au= ( dble(i_au)*x_e(i) + dble(N_au-i_au)*x_s(i) )/dble(N_au)
      y_au= ( dble(i_au)*y_e(i) + dble(N_au-i_au)*y_s(i) )/dble(N_au)
      call nearest_id( 0, L, 0, M, xr, yr, rmask, x_au, y_au, Iau, Jau )
      if(idaq(i)==1) then
        data(Iau,Jau,1) = data(Iau,Jau,1) + 1.0d0/area(Iau,Jau)
      elseif(idaq(i)==2) then
        data(Iau,Jau,2) = data(Iau,Jau,2) + 1.0d0/area(Iau,Jau)
      elseif(idaq(i)==3) then
        data(Iau,Jau,3) = data(Iau,Jau,3) + 1.0d0/area(Iau,Jau)
      elseif(idaq(i)==4) then
        data(Iau,Jau,4) = data(Iau,Jau,4) + 1.0d0/area(Iau,Jau)
      elseif(idaq(i)==6) then
        data(Iau,Jau,3) = data(Iau,Jau,3) + 0.5d0/area(Iau,Jau)
        data(Iau,Jau,4) = data(Iau,Jau,4) + 0.5d0/area(Iau,Jau)
      endif
    end do
  end do
!$omp end do
!$omp end parallel

  Nid = Nidaq

!=========== AQUACULTURE2 (GINZAKE) ===============================
#elif defined AQUACULTURE2

  allocate( data(0:L,0:M, 1) )
  data(:,:,:) = 0.0d0

  write(*,*) "OPEN: ", trim( INPUT_CSV )
  open(10,file=INPUT_CSV, status='old')

  !--- Count data -----
  Ndata = 0
  do i=1, Nheader
    read(10,*)
  end do
  do
    read(10,*,iostat=ios)
    if(ios==-1) exit
    Ndata = Ndata + 1
  end do
  allocate( idaq(Ndata), lat_c(Ndata), lon_c(Ndata), area_aq(Ndata) )
  allocate( x_c(Ndata), y_c(Ndata) )
  rewind(10) 
  !--- Read data -----
  do i=1, Nheader
    read(10,*)
  end do
  do i=1, Ndata
    ! read aquaculture1 data
    read(10,*) str_idaq, lat_c(i), lon_c(i), area_aq(i)
    ! convert characters to integer
    read (str_idaq, *) idaq(i)
!    write(*,*) idaq(i), lat_c(i), lon_c(i), area_aq(i)
    ! convert lat lon coordinate to UTM coordinate
    call ll2utm(lat_c(i),lon_c(i),x_c(i),y_c(i),izone,ispher)

  end do
  close(10)
  write(*,*) "CLOSE: ", trim( INPUT_CSV )

  write(*,*) "Processing... "
!$omp parallel
!$omp do private(i,Iau,Jau)
  do i=1, Ndata
    call nearest_id( 0, L, 0, M, xr, yr, rmask, x_c(i), y_c(i), Iau, Jau )
    if(idaq(i)==5) then
      data(Iau,Jau,1) = data(Iau,Jau,1) + area_aq(i)/area(Iau,Jau)
    endif
  end do
!$omp end do
!$omp end parallel

  Nid = 1

#endif
!==================================================================
  ! Set zero at boundary
  data(:,0,:) = 0.0d0
  data(:,M,:) = 0.0d0
  data(0,:,:) = 0.0d0
  data(L,:,:) = 0.0d0

!===== Add parameters into ROMS grid netCDF file ===========================
  write(*,*) "OPEN: ", trim( GRID_FILE )

  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), NF90_WRITE, ncid) )

  call check( nf90_inq_dimid(ncid, 'xi_rho' , xr_dimid) )
  call check( nf90_inq_dimid(ncid, 'eta_rho', yr_dimid) )
  dim2Dids = (/ xr_dimid, yr_dimid /)

  allocate( varname(Nid) )

  do j=1,Nid
    if(Nid==1) then
      varname(j) = trim( NC_VARNAME )
    else
      write(varnum,'(I2.2)') j
      varname(j) = trim( NC_VARNAME )//varnum
    endif
    status = nf90_inq_varid(ncid, trim( varname(j) ), var_id)
    if (status /= nf90_noerr) then
      write(*,*) 'Add variable: ', trim( varname(j) )
      call check( nf90_redef(ncid) )
      call check( nf90_def_var(ncid, trim( varname(j) ), NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LONGNAME ) ) )
      call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT ) ) )
!      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      call check( nf90_enddef(ncid) )
    else
      write(*,*) 'Found variable: ', trim( varname(j) )
    endif
  enddo

!  start2D = (/ 1,  1 /)
!  count2D = (/ Nxr, Nyr /)
  do j=1,Nid
    write(*,*)  'Write: ', trim(varname(j))
    call check( nf90_inq_varid(ncid, varname(j), var_id) )
!    call check( nf90_put_var(ncid, var_id, data, start = start2D, count = count2D) )
    call check( nf90_put_var(ncid, var_id, data(:,:,j)) )
   
  enddo
  call check( nf90_close(ncid) )

! ===== Output SWAN grid file =========================
#if defined OUTPUT_SWAN_GRID
  mxc = L
  myc = M
  
  do i=1,Nid
    write(varnum,'(I2.2)') i
    SWAN_INPGRID_FILE = trim( SWAN_INPGRID_PREFIX )//"_"//varnum//".dat"
    write(*,*) "WRITE: ", trim( SWAN_INPGRID_FILE )  
    open(10,file = SWAN_INPGRID_FILE)
    do j=0,myc
      write(10,*) data(0:mxc,j,i)
    enddo
    close(10)
  enddo

#endif  

  write(*,*) 'FINISH!!'
      
END PROGRAM grdROMS_add_ncparams
      
