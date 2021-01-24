
!!!=== Copyright (c) 2019-2021 Takashi NAKAMURA  =====

PROGRAM grdROMS_update
  use netcdf
  use mod_utility
  use mod_roms_netcdf
 
  implicit none
  
! -------------------------------------------------------------------------
  character(256) :: GRID_FILE, GRID_FILE2
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

#if defined GRID_FINE2COARSE
  real(8), allocatable :: h_rg(:,:)        ! depth (meter)
  character(256) :: parent_grid
  integer :: parent_Imin, parent_Imax
  integer :: parent_Jmin, parent_Jmax
  integer :: refine_factor
  character(256) :: F2C_OUT_FILE
  integer :: Nxr_rg, Nyr_rg
  integer :: L_rg, M_rg
  integer :: i_rg, j_rg, i0_rg, j0_rg
#endif

  namelist/grd/GRID_FILE
#if defined BATH_SMOOTHING
  namelist/bath_smooth/rx0max
#endif
  namelist/land_mask/hmin
  namelist/land_mask/I_ocn, J_ocn
#if defined GRID_FINE2COARSE
  namelist/refinement/parent_grid
  namelist/refinement/parent_Imin, parent_Imax
  namelist/refinement/parent_Jmin, parent_Jmax
  namelist/refinement/refine_factor
  namelist/fine2coarse/F2C_OUT_FILE
#endif

  ! Read parameters in namelist file
  
  read (*, nml=grd)
#if defined BATH_SMOOTHING
  rewind(5)
  read (5, nml=bath_smooth)
#endif
  rewind(5)
  read (5, nml=land_mask)
#if defined GRID_FINE2COARSE
  rewind(5)
  read (5, nml=refinement)
  rewind(5)
  read (5, nml=fine2coarse)
#endif

#if defined GRID_FINE2COARSE
!-Read ROMS fine grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )

  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr_rg)
  call get_dimension(ncid, 'eta_rho', Nyr_rg)
 
  L_rg = Nxr_rg-1
  M_rg = Nyr_rg-1
  
  allocate( h_rg(0:L_rg, 0:M_rg) )
  
  ! Get variables
  start2D = (/ 1,  1 /)
  count2D = (/ Nxr_rg, Nyr_rg /)
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, h_rg) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )

  write(*,*) "CLOSE: ", trim( GRID_FILE )

#endif

#if defined GRID_FINE2COARSE 
!-Copy and set coarse ROMS grid netCDF file ----------------
  write(*,*) "CLEATE: ", trim( F2C_OUT_FILE ) 
  call system("cp "//trim( parent_grid )//" "//trim( F2C_OUT_FILE ))
  GRID_FILE2 = F2C_OUT_FILE

#else
!-Set ROMS grid netCDF file --------------------------------
  GRID_FILE2 = GRID_FILE
#endif

!-Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE2 )

  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE2 ), nf90_nowrite, ncid) )
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

  write(*,*) "CLOSE: ", trim( GRID_FILE2 )

#if defined GRID_FINE2COARSE 
!-Merge coase grid with fine grid ---------------------------------------

  i0_rg = int( refine_factor/2 +1 )
  j0_rg = int( refine_factor/2 +1 )

  do i=parent_Imin, parent_Imax-1
    do j=parent_Jmin, parent_Jmax-1
      if( i<0 .or. i>L .or. j<0 .or. j>M ) cycle
      i_rg = i0_rg + refine_factor*(i-parent_Imin)
      j_rg = j0_rg + refine_factor*(j-parent_Jmin)

      h(i,j) = h_rg(i_rg,j_rg)
      
    enddo
  enddo  

#endif   
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
  call writeNetCDF_2d('h', trim( GRID_FILE2 )       &
    , Nxr, Nyr, h, start2D, count2D )
  call writeNetCDF_2d('mask_rho', trim( GRID_FILE2 )   &
    , Nxr, Nyr, rmask, start2D, count2D  )
  
  count2D = (/ Nxp, Nyp /)
  call writeNetCDF_2d('mask_psi', trim( GRID_FILE2 )   &
        , Nxp, Nyp, pmask, start2D, count2D )
  
  count2D = (/ Nxu, Nyu /)
  call writeNetCDF_2d('mask_u', trim( GRID_FILE2 )    &
        , Nxu, Nyu, umask, start2D, count2D )
     
  count2D = (/ Nxv, Nyv /)
  call writeNetCDF_2d('mask_v', trim( GRID_FILE2 )    &
        , Nxv, Nyv, vmask, start2D, count2D )

  write(*,*) 'FINISH!!'
      
END PROGRAM grdROMS_update
      
