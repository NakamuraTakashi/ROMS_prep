
!!!=== Copyright (c) 2019-2020 Takashi NAKAMURA  =====

PROGRAM grdROMS_LandMasking
  use netcdf
  use mod_utility
  use mod_roms_netcdf
 
  implicit none
  
! -------------------------------------------------------------------------
  character(256) :: GRID_FILE
  real(8) :: hmin
  real(8) :: rx0max
  integer :: I_ocn, J_ocn
! -------------------------------------------------------------------------
  integer :: start2D(2), count2D(2)
  
  real(8), allocatable :: h(:,:)        ! depth (meter)
  real(8), allocatable :: rmask(:, :)
  real(8), allocatable :: pmask(:, :)
  real(8), allocatable :: umask(:, :)
  real(8), allocatable :: vmask(:, :)
  
  integer :: i,j

  integer :: Nxr, Nyr
  integer :: Nxp, Nyp
  integer :: Nxu, Nyu
  integer :: Nxv, Nyv
  integer :: L, M
 
  integer :: ncid,var_id
  
  namelist/grd/GRID_FILE
  namelist/bath/hmin
  namelist/bath/rx0max
  namelist/bath/I_ocn, J_ocn
  ! Read parameters in namelist file
  
  read (*, nml=grd)
  read (*, nml=bath)
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
  
  allocate( rmask(0:L, 0:M) )
  allocate( umask(1:L, 0:M) )
  allocate( vmask(0:L, 1:M) )
  allocate( pmask(1:L, 1:M) )
  allocate( h(0:L, 0:M) )
 
  
  ! Get variables
  start2D = (/ 1,  1 /)
  count2D = (/ Nxr, Nyr /)
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, h) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )

  write(*,*) "CLOSE: ", trim( GRID_FILE )
    
!---- Compute Lat, Lon fields for ROMS grid netCDF file --------------------------------
  
  ! Compute land mask
  rmask(:,:) = 1.0d0
  pmask(:,:) = 1.0d0
  umask(:,:) = 1.0d0
  vmask(:,:) = 1.0d0

  CALL land_masking(Nxr, Nyr, h, hmin, rmask)
#if defined BATH_SMOOTHING
  ! Smoothing Bathymetry
  CALL bathy_smooth(Nxr, Nyr, rx0max, rmask, h)
#endif

  CALL isolated_water_masking(Nxr, Nyr, I_ocn, J_ocn, rmask)

  CALL uvp_masks(Nxr, Nyr, rmask, umask, vmask, pmask)

  do i=0, L
    do j=0, M
      if(rmask(i,j) == 0.0d0) then
        h(i,j) = -10.0d0
      endif
      if(h(i,j) == 0.0d0) then
        h(i,j) = 0.001d0
      endif
    enddo
  enddo
      
!---- Write ROMS grid netCDF file --------------------------------
  start2D = (/ 1,  1 /)
  count2D = (/ Nxr, Nyr /)
  call writeNetCDF_2d('h', trim( GRID_FILE )       &
    , Nxr, Nyr, h, start2D, count2D )
  call writeNetCDF_2d('mask_rho', trim( GRID_FILE )   &
    , Nxr, Nyr, rmask, start2D, count2D  )
  
  count2D = (/ Nxp, Nyp /)
  call writeNetCDF_2d('mask_psi', trim( GRID_FILE )   &
        , Nxp, Nyp, pmask, start2D, count2D )
  
  count2D = (/ Nxu, Nyu /)
  call writeNetCDF_2d('mask_u', trim( GRID_FILE )    &
        , Nxu, Nyu, umask, start2D, count2D )
     
  count2D = (/ Nxv, Nyv /)
  call writeNetCDF_2d('mask_v', trim( GRID_FILE )    &
        , Nxv, Nyv, vmask, start2D, count2D )

  write(*,*) 'FINISH!!'
      
END PROGRAM grdROMS_LandMasking
      
