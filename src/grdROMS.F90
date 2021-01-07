
!!!=== Copyright (c) 2018-2021 Takashi NAKAMURA  =====

PROGRAM grdROMS
  use netcdf
  use mod_utility
  use mod_roms_netcdf
  use mod_interpolation
 
  implicit none
  
! -------------------------------------------------------------------------
  character(256) :: BATH_FILE
  character(256) :: GRID_FILE
  real(8) :: s_lat, e_lat, s_lon, e_lon
  real(8) :: RESOLUTION, angle
  real(8) :: hmin
  real(8) :: rx0max
  integer :: spherical
  integer :: I_ocn, J_ocn
#if defined GRID_REFINEMENT
  character(256) :: parent_grid
  integer :: parent_Imin, parent_Imax
  integer :: parent_Jmin, parent_Jmax
  integer :: refine_factor
  integer, parameter :: Nmerge = 15
  real(8), parameter :: fmerge(Nmerge) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0  &
       , 1.0d0, 1.0d0, 0.9d0, 0.8d0, 0.7d0, 0.6d0, 0.5d0, 0.4d0, 0.3d0  &
       , 0.2d0, 0.1d0 /)
#endif
! ---------------------------------------------------------------------

  real(8), allocatable :: lat_all(:), lon_all(:)
  real(8), allocatable :: lat(:), lon(:), depth(:,:)
  integer :: start1D(1), count1D(1)
  integer :: start2D(2), count2D(2)
  
  real(8), allocatable :: h(:,:)        ! depth (meter)
  real(8), allocatable :: hraw(:,:)        ! depth (meter)
  real(8), allocatable :: f(:,:)
  real(8), allocatable :: pm(:,:)
  real(8), allocatable :: pn(:,:)
  real(8), allocatable :: dndx(:,:)
  real(8), allocatable :: dmde(:,:)
  real(8), allocatable :: angler(:,:)
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: rmask(:, :)
  real(8), allocatable :: latp(:, :)
  real(8), allocatable :: lonp(:, :)
  real(8), allocatable :: pmask(:, :)
  real(8), allocatable :: latu(:, :)
  real(8), allocatable :: lonu(:, :)
  real(8), allocatable :: umask(:, :)
  real(8), allocatable :: latv(:, :)
  real(8), allocatable :: lonv(:, :)
  real(8), allocatable :: vmask(:, :)
#if defined UTM_COORD
  real(8), allocatable :: xr(:, :)
  real(8), allocatable :: yr(:, :)
  real(8), allocatable :: xp(:, :)
  real(8), allocatable :: yp(:, :)
  real(8), allocatable :: xu(:, :)
  real(8), allocatable :: yu(:, :)
  real(8), allocatable :: xv(:, :)
  real(8), allocatable :: yv(:, :)
#endif
#if defined GRID_REFINEMENT
  real(8), allocatable :: hr(:,:)        ! depth (meter)
  real(8), allocatable :: rmask_r(:,:)        ! depth (meter)
  real(8), allocatable :: xr0(:,:) 
  real(8), allocatable :: yr0(:,:)
   
  real(8), allocatable :: Xd(:)
  real(8), allocatable :: Yd(:)
  real(8), allocatable :: P_h(:,:)        ! depth (meter)
#if defined UTM_COORD
  real(8), allocatable :: P_xr(:, :)
  real(8), allocatable :: P_yr(:, :)
#endif
  real(8), allocatable :: P_latr(:, :)
  real(8), allocatable :: P_lonr(:, :)
  real(8), allocatable :: P_rmask(:, :)
  real(8), allocatable :: P_latp(:, :)
  real(8), allocatable :: P_lonp(:, :)
  real(8), allocatable :: P_latu(:, :)
  real(8), allocatable :: P_lonu(:, :)
  real(8), allocatable :: P_latv(:, :)
  real(8), allocatable :: P_lonv(:, :)
  real(8), allocatable :: P_angler(:, :)
#endif

  real(8) :: d_lat, d_lon
  real(8) :: latr_min, latr_max, lonr_min, lonr_max
  real(8) :: sf,off

  integer :: i,j
  
  integer :: N_lat_all, N_lon_all
  integer :: N_lat, N_lon
  integer :: id_slat, id_slon, id_elat, id_elon
  
  integer :: Nxr, Nyr, Nxp, Nyp, Nxu, Nyu, Nxv, Nyv
#if defined GRID_REFINEMENT
  integer :: P_Nxr, P_Nyr, P_Nxp, P_Nyp, P_Nxu, P_Nyu, P_Nxv, P_Nyv
  integer :: ip,jp,is,js
#endif
     
  integer :: ncid,var_id
  
#if defined UTM_COORD
  real(8) :: conv
  integer :: izone,ispher
#endif

  namelist/grd/GRID_FILE
#if defined GEBCO2ROMS
  namelist/gebco/BATH_FILE
#elif defined EMODNET2ROMS
  namelist/emodnet/BATH_FILE
#endif
#if defined BATH_SMOOTHING
  namelist/bath_smooth/rx0max
#endif
  namelist/land_mask/hmin
  namelist/land_mask/I_ocn, J_ocn
#if defined GRID_REFINEMENT
  namelist/refinement/parent_grid
  namelist/refinement/parent_Imin, parent_Imax
  namelist/refinement/parent_Jmin, parent_Jmax
  namelist/refinement/refine_factor
#else
  namelist/grd_setting_ll/s_lat, e_lat, s_lon, e_lon
  namelist/grd_setting_ll/RESOLUTION, angle
#endif
#if defined UTM_COORD
  namelist/utm_zone/izone,ispher
#endif
  namelist/hcoord/spherical

  ! Read parameters in namelist file
  
  read (5, nml=grd)
#if defined GEBCO2ROMS
  rewind(5)
  read (5, nml=gebco)
#elif defined EMODNET2ROMS
  rewind(5)
  read (5, nml=emodnet)
#endif
#if defined BATH_SMOOTHING
  rewind(5)
  read (5, nml=bath_smooth)
#endif
  rewind(5)
  read (5, nml=land_mask)
#if defined GRID_REFINEMENT
  rewind(5)
  read (5, nml=refinement)
#else
  rewind(5)
  read (5, nml=grd_setting_ll)
#endif
#if defined UTM_COORD
  rewind(5)
  read (5, nml=utm_zone)
#endif
  rewind(5)
  read (5, nml=hcoord)

!---- Compute Lat, Lon fields for ROMS grid netCDF file --------------------------------
#if defined GRID_REFINEMENT
  Nxr = (parent_Imax-parent_Imin)*refine_factor+2
  Nyr = (parent_Jmax-parent_Jmin)*refine_factor+2 
#else      
!  Nxr = int( (e_lon - s_lon)*RESOLUTION ) + 1
!  Nyr = int( (e_lat - s_lat)*RESOLUTION ) + 1
  call grid_size_rectangular(s_lon, s_lat, e_lon, e_lat, RESOLUTION, angle, Nxr, Nyr)
#endif
  Nxp = Nxr-1
  Nyp = Nyr-1
  Nxu = Nxr-1
  Nyu = Nyr
  Nxv = Nxr
  Nyv = Nyr-1
!
  write(*,*) 'CHECK:',Nxr, Nyr
  
  allocate(h(Nxr, Nyr))
  allocate(hraw(Nxr, Nyr))
  allocate(f(Nxr, Nyr))
  allocate(pm(Nxr, Nyr))
  allocate(pn(Nxr, Nyr))
  allocate(dndx(Nxr, Nyr))
  allocate(dmde(Nxr, Nyr))
  allocate(angler(Nxr, Nyr))
  allocate(latr(Nxr, Nyr))
  allocate(lonr(Nxr, Nyr))
  allocate(rmask(Nxr, Nyr))
  allocate(latp(Nxp, Nyp))
  allocate(lonp(Nxp, Nyp))
  allocate(pmask(Nxp, Nyp))
  allocate(latu(Nxu, Nyu))
  allocate(lonu(Nxu, Nyu))
  allocate(umask(Nxu, Nyu))
  allocate(latv(Nxv, Nyv))
  allocate(lonv(Nxv, Nyv))
  allocate(vmask(Nxv, Nyv))
#if defined UTM_COORD
  allocate(xr(Nxr, Nyr))
  allocate(yr(Nxr, Nyr))
  allocate(xp(Nxp, Nyp))
  allocate(yp(Nxp, Nyp))
  allocate(xu(Nxu, Nyu))
  allocate(yu(Nxu, Nyu))
  allocate(xv(Nxv, Nyv))
  allocate(yv(Nxv, Nyv))
#endif
#if defined GRID_REFINEMENT
  allocate(hr(Nxr, Nyr))
  allocate(rmask_r(Nxr, Nyr))
  allocate(xr0(Nxr, Nyr))
  allocate(yr0(Nxr, Nyr))
#endif

#if defined GRID_REFINEMENT
!-Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim(parent_grid)
  
  ! Open NetCDF grid file
  call check( nf90_open( trim(parent_grid), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  P_Nxr)
  call get_dimension(ncid, 'eta_rho', P_Nyr)
  call get_dimension(ncid, 'xi_u',  P_Nxu)
  call get_dimension(ncid, 'eta_u', P_Nyu)
  call get_dimension(ncid, 'xi_v',  P_Nxv)
  call get_dimension(ncid, 'eta_v', P_Nyv)
# if defined UTM_COORD
  allocate(P_xr(P_Nxr, P_Nyr))
  allocate(P_yr(P_Nxr, P_Nyr))
# endif
  allocate(P_latr(P_Nxr, P_Nyr))
  allocate(P_lonr(P_Nxr, P_Nyr))
  allocate(P_latu(P_Nxu, P_Nyu))
  allocate(P_lonu(P_Nxu, P_Nyu))
  allocate(P_latv(P_Nxv, P_Nyv))
  allocate(P_lonv(P_Nxv, P_Nyv))
  allocate(P_h(P_Nxr, P_Nyr))
  allocate(P_angler(P_Nxr, P_Nyr))
  allocate(P_rmask(P_Nxr, P_Nyr))
  allocate(Xd(P_Nxr), Yd(P_Nyr))
        
  ! Get variable id
# if defined UTM_COORD
  call check( nf90_inq_varid(ncid, 'x_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, P_xr) )
  call check( nf90_inq_varid(ncid, 'y_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, P_yr) )
# endif
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, P_latr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, P_lonr) )
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, P_h) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, P_rmask) )
  call check( nf90_inq_varid(ncid, 'angle', var_id) ) 
  call check( nf90_get_var(ncid, var_id, P_angler) )

  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
  write(*,*) "CLOSE: ", trim(parent_grid)

  do i=1, P_Nxr
    Xd(i)=dble(i-1)
  enddo
  do j=1, P_Nyr
    Yd(j)=dble(j-1)
  enddo
  do i=1, Nxr
    xr0(i,:)=dble(parent_Imin-1) + dble(i-1+refine_factor/2)/dble(refine_factor)
  enddo
  do j=1, Nyr
    yr0(:,j)=dble(parent_Jmin-1) + dble(j-1+refine_factor/2)/dble(refine_factor)
  enddo

# if defined UTM_COORD
  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_xr , 0.0d0, 1.0d9                  &
                    , Nxr, Nyr, xr0, yr0, xr )
  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_yr , 0.0d0, 1.0d9                  &
                    , Nxr, Nyr, xr0, yr0, yr )
  ! UTM -> lat lon
  do i=1,Nxr
    do j=1,Nyr
      call utm2ll( xr(i,j), yr(i,j), izone, ispher         &
                 , latr(i,j), lonr(i,j), conv)
    enddo
  enddo

# else

  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_latr , -10000.0d0, 10000.0d0       &
                    , Nxr, Nyr, xr0, yr0, latr )
  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_lonr , -10000.0d0, 10000.0d0       &
                    , Nxr, Nyr, xr0, yr0, lonr )
# endif

  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_h , -9.0d0, 10000.0d0              &
                    , Nxr, Nyr, xr0, yr0, hr )
  call LinearInterpolation2D_grid2(P_Nxr, P_Nyr, Xd, Yd    &
                    , P_angler , -10.0d0, 10.0d0           &
                    , Nxr, Nyr, xr0, yr0, angler )

  do ip=parent_Imin+1, parent_Imax
    do jp=parent_Jmin+1, parent_Jmax
      is = (ip-parent_Imin-1)*refine_factor+2
      js = (jp-parent_Jmin-1)*refine_factor+2
      do i=is, is+refine_factor
        do j=js, js+refine_factor
          rmask_r(i,j) = P_rmask(ip,jp)
        enddo
      enddo
    enddo
  enddo
  rmask_r(1,:)=rmask_r(2,:)
  rmask_r(Nxr,:)=rmask_r(Nxr-1,:)
  rmask_r(:,1)=rmask_r(:,2)
  rmask_r(:,Nyr)=rmask_r(:,Nyr-1)

#else      

!  d_lat = 1.0d0/RESOLUTION
!  d_lon = 1.0d0/RESOLUTION
!
!  do i=1, Nxr
!    lonr(i,:)=s_lon + d_lon*dble(i-1)
!  enddo
!  do j=1, Nyr
!    latr(:,j)=s_lat + d_lat*dble(j-1)
!  enddo
  call grid_gen_rectangular( s_lon, s_lat, RESOLUTION, angle &
                           , Nxr, Nyr, lonr, latr )
  
  angler(:,:) = angle

#endif

#if defined UTM_COORD
  do i=1, Nxu
    do j=1, Nyu
      xu(i,j)=0.5d0*(xr(i,j)+xr(i+1,j))
      yu(i,j)=0.5d0*(yr(i,j)+yr(i+1,j))
      call utm2ll( xu(i,j), yu(i,j), izone, ispher           &
                 , latu(i,j), lonu(i,j), conv )
    enddo      
  enddo
  do i=1, Nxv
    do j=1, Nyv
      xv(i,j)=0.5d0*(xr(i,j)+xr(i,j+1))
      yv(i,j)=0.5d0*(yr(i,j)+yr(i,j+1))
      call utm2ll( xv(i,j), yv(i,j), izone, ispher           &
                 , latv(i,j), lonv(i,j), conv )
    enddo      
  enddo
  do i=1, Nxp
    do j=1, Nyp
      xp(i,j)=0.5d0*(xu(i,j)+xu(i,j+1))
      yp(i,j)=0.5d0*(yv(i,j)+yv(i+1,j))
      call utm2ll( xp(i,j), yp(i,j), izone, ispher           &
                 , latp(i,j), lonp(i,j), conv )
    enddo      
  enddo

  ! compute pm, pn
  do i=2, Nxr-1
    do j=1, Nyr
      pm(i,j) = 1.0d0/(xu(i,j)-xu(i-1,j))
    enddo
  enddo
  pm(1,:) = pm(2,:) 
  pm(Nxr,:) = pm(Nxr-1,:) 

  do i=1, Nxr
    do j=2, Nyr-1
      pn(i,j) = 1.0d0/(yv(i,j)-yv(i,j-1))
    enddo
  enddo
  pn(:,1) = pn(:,2) 
  pn(:,Nyr) = pn(:,Nyr-1)

#else
  do i=1, Nxu
    do j=1, Nyu
      lonu(i,j)=0.5d0*(lonr(i,j)+lonr(i+1,j))
      latu(i,j)=0.5d0*(latr(i,j)+latr(i+1,j))
    enddo      
  enddo
  do i=1, Nxv
    do j=1, Nyv
      lonv(i,j)=0.5d0*(lonr(i,j)+lonr(i,j+1))
      latv(i,j)=0.5d0*(latr(i,j)+latr(i,j+1))
    enddo      
  enddo
  do i=1, Nxp
    do j=1, Nyp
      lonp(i,j)=0.5d0*(lonu(i,j)+lonu(i,j+1))
      latp(i,j)=0.5d0*(latv(i,j)+latv(i+1,j))
    enddo      
  enddo

  ! compute pm, pn
  do i=2, Nxr-1
    do j=1, Nyr
      pm(i,j) = 1.0d0/distance( latr(i,j), lonu(i-1,j)   &
                              , latr(i,j), lonu(i,j))
    enddo
  enddo
  pm(1,:) = pm(2,:) 
  pm(Nxr,:) = pm(Nxr-1,:) 

  do i=1, Nxr
    do j=2, Nyr-1
      pn(i,j) = 1.0d0/distance( latv(i,j-1), lonr(i,j)   &
                              , latv(i,j), lonr(i,j))
    enddo
  enddo
  pn(:,1) = pn(:,2) 
  pn(:,Nyr) = pn(:,Nyr-1)

#endif

  ! Coriolis parameter
  do i=1, Nxr
    do j=1, Nyr
      f(i,j) = Coriolis(latr(i,j))
    enddo
  enddo
  dndx(:,:) = 0.0d0
  dmde(:,:) = 0.0d0

#if defined GEBCO2ROMS || defined EMODNET2ROMS
!---- Read GEBCO grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( BATH_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( BATH_FILE ), nf90_nowrite, ncid) )
# if defined GEBCO2ROMS
  ! Get dimension data
  call get_dimension(ncid, 'lat', N_lat_all)
  call get_dimension(ncid, 'lon', N_lon_all)
  ! Allocate variable
  allocate(lat_all(N_lat_all))
  allocate(lon_all(N_lon_all))
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat', var_id) )
  call check( nf90_get_var(ncid, var_id, lat_all) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) )
  call check( nf90_get_var(ncid, var_id, lon_all) )

# elif defined EMODNET2ROMS
  call get_dimension(ncid, 'LINES', N_lat_all)
  call get_dimension(ncid, 'COLUMNS', N_lon_all)
  ! Allocate variable
  allocate(lat_all(N_lat_all))
  allocate(lon_all(N_lon_all))
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'LINES', var_id) )
  call check( nf90_get_var(ncid, var_id, lat_all) )
  call check( nf90_inq_varid(ncid, 'COLUMNS', var_id) )
  call check( nf90_get_var(ncid, var_id, lon_all) )
# endif

  ! Trim GEBCO/EMODnet bathymetry
# if defined GRID_REFINEMENT
!  s_lat = latr(1,1)
!  e_lat = latr(Nxr,Nyr)
!  s_lon = lonr(1,1)
!  e_lon = lonr(Nxr,Nyr)
# endif
!
!  do i=2, N_lat_all
!    if(s_lat < lat_all(i)) then
!      id_slat = i-1
!      exit
!    endif
!  enddo
!  do i=1, N_lat_all
!    if(e_lat < lat_all(i)) then
!      id_elat = i
!      exit
!    endif
!  enddo
!  do i=2, N_lon_all
!    if(s_lon < lon_all(i)) then
!      id_slon = i-1
!      exit
!    endif
!  enddo
!  do i=1, N_lon_all
!    if(e_lon < lon_all(i)) then
!      id_elon = i
!      exit
!    endif
!  enddo

  call min_max_2D( Nxr, Nyr, latr, latr_min, latr_max)
  call min_max_2D( Nxr, Nyr, lonr, lonr_min, lonr_max)
  call seek_id_range( N_lat_all, lat_all, latr_min, latr_max, id_slat, id_elat)
  call seek_id_range( N_lon_all, lon_all, lonr_min, lonr_max, id_slon, id_elon)
  
  N_lon = id_elon-id_slon+1
  N_lat = id_elat-id_slat+1
  allocate(depth(N_lon,N_lat))
  start2D = (/ id_slon, id_slat /)
  count2D = (/ N_lon  , N_lat   /)
  write(*,*) 'CHECK:',id_slon,id_elon, id_slat,id_elat, N_lon  , N_lat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'CHECK:',lon_all(id_slon),lon_all(id_elon), lat_all(id_slat),lat_all(id_elat)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Get variable id
# if defined GEBCO2ROMS
  call check( nf90_inq_varid(ncid, 'elevation', var_id) )
  call check( nf90_get_var(ncid, var_id, depth, start=start2D, count=count2D) )

# elif defined EMODNET2ROMS
  call check( nf90_inq_varid(ncid, 'DEPTH', var_id) )
  call check( nf90_get_var(ncid, var_id, depth, start=start2D, count=count2D) )
  call check( nf90_get_att(ncid, var_id, 'scale_factor', sf) )
  call check( nf90_get_att(ncid, var_id, 'add_offset', off) )
  depth = depth*sf+off
# endif

  ! Close NetCDF file
  call check( nf90_close(ncid) )
  write(*,*) "CLOSE: ", trim( BATH_FILE )
  
  depth(:,:) = -1.0d0*depth(:,:)
  
  allocate(lon(N_lon))
  allocate(lat(N_lat))
  do i=1, N_lon
    lon(i) = lon_all(id_slon+i-1)
  enddo
  do i=1, N_lat
    lat(i) = lat_all(id_slat+i-1)
  enddo
  write(*,*) 'CHECK:',lon(1),lon(N_lon), lat(1),lat(N_lat), depth(1,1),depth(N_lon,N_lat)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'CHECK:',lonr(1,1),lonr(Nxr, Nyr),latr(1,1),latr(Nxr, Nyr)

! --- Linear interporation --------------------------------

  write(*,*) 'Linear Interporation: h'
  call LinearInterpolation2D_grid2(N_lon, N_lat, lon, lat      &
                  , depth , -10000.0d0, 10000.0d0              &
                  , Nxr, Nyr, lonr, latr, h )

  write(*,*) '*** SUCCESS Linear Interporation'
  
#else

  h(:,:) = hr(:,:)
  
#endif

  hraw(:,:) = h(:,:)

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

#if defined GRID_REFINEMENT
  ! Merge GEBCO grid with refined ROMS grid
  do j=1, Nmerge
    do i=1+j-1, Nxr-j+1
      h(i,j) = fmerge(j)*hr(i,j) + (1.0d0-fmerge(j))*h(i,j)
      h(i,Nyr-j+1) = fmerge(j)*hr(i,Nyr-j+1) + (1.0d0-fmerge(j))*h(i,Nyr-j+1)
    enddo
  enddo
  do i=1, Nmerge
    do j=1+i-1, Nyr-i+1
      h(i,j) = fmerge(i)*hr(i,j) + (1.0d0-fmerge(i))*h(i,j)
      h(Nxr-i+1,j) = fmerge(i)*hr(Nxr-i+1,j) + (1.0d0-fmerge(i))*h(Nxr-i+1,j)
    enddo
  enddo
  do j=1, refine_factor/2+1
    rmask(:,j) = rmask_r(:,j)
    rmask(:,Nyr-j+1) = rmask_r(:,Nyr-j+1)
  enddo
  do i=1, refine_factor/2+1
    rmask(i,:) = rmask_r(i,:)
    rmask(Nxr-i+1,:) = rmask_r(Nxr-i+1,:)
  enddo
  
#endif

!---------------------------------------------------------------------  

  CALL isolated_water_masking(Nxr, Nyr, I_ocn, J_ocn, rmask)

  CALL uvp_masks(Nxr, Nyr, rmask, umask, vmask, pmask)

  do i=1, Nxr
    do j=1, Nyr
      if(rmask(i,j) == 0.0d0) then
        h(i,j) = -10.0d0
      endif
      if(h(i,j) == 0.0d0) then
        h(i,j) = 0.001d0
      endif
    enddo
  enddo
  
!
!---- Create the ROMS grid netCDF file --------------------------------
      
  call createNetCDFgrd( trim( GRID_FILE ), Nxr, Nyr )
  
!---- Write ROMS grid netCDF file --------------------------------

  call check( nf90_open(trim( GRID_FILE ), NF90_WRITE, ncid) )
! Put global attribute
  call check( nf90_redef(ncid) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'Bath_grid', BATH_FILE) )
#if defined GRID_REFINEMENT
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_grid', parent_grid) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_Imin', parent_Imin) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_Imax', parent_Imax) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_Jmin', parent_Jmin) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_Jmax', parent_Jmax) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, 'refine_factor', refine_factor) )
#endif
  call check( nf90_enddef(ncid) )

  call check( nf90_inq_varid(ncid, 'spherical', var_id) )
  call check( nf90_put_var(ncid, var_id, spherical) )
  ! call check( nf90_inq_varid(ncid, 'xl', var_id) )
  ! call check( nf90_put_var(ncid, var_id, xl) )
  ! call check( nf90_inq_varid(ncid, 'el', var_id) )
  ! call check( nf90_put_var(ncid, var_id, el) )
  call check( nf90_close(ncid) )

  start2D = (/ 1,  1 /)
  count2D = (/ Nxr, Nyr /)
  call writeNetCDF_2d( 'h', trim( GRID_FILE )        &
        , Nxr, Nyr, h, start2D, count2D )
  call writeNetCDF_2d( 'f', trim( GRID_FILE )        &
        , Nxr, Nyr, f, start2D, count2D )
  call writeNetCDF_2d( 'pm', trim( GRID_FILE )       &
        , Nxr, Nyr, pm, start2D, count2D )
  call writeNetCDF_2d( 'pn', trim( GRID_FILE )       &
        , Nxr, Nyr, pn, start2D, count2D )
  call writeNetCDF_2d( 'dndx', trim( GRID_FILE )     &
        , Nxr, Nyr, dndx, start2D, count2D )
  call writeNetCDF_2d( 'dmde', trim( GRID_FILE )     &
        , Nxr, Nyr, dmde, start2D, count2D )
  call writeNetCDF_2d( 'angle', trim( GRID_FILE )    &
        , Nxr, Nyr, angler, start2D, count2D )
  call writeNetCDF_2d( 'lat_rho', trim( GRID_FILE )  &
        , Nxr, Nyr, latr, start2D, count2D )
  call writeNetCDF_2d( 'lon_rho', trim( GRID_FILE )  &
        , Nxr, Nyr, lonr, start2D, count2D )
  call writeNetCDF_2d( 'mask_rho', trim( GRID_FILE ) &
        , Nxr, Nyr, rmask, start2D, count2D )
  
  call writeNetCDF_2d( 'hraw', trim( GRID_FILE )     &
        , Nxr, Nyr, hraw, start2D, count2D )

#if defined GRID_REFINEMENT
  call writeNetCDF_2d( 'mask_rho_coarse', trim( GRID_FILE )  &
        , Nxr, Nyr, rmask_r, start2D, count2D )
#endif
#if defined UTM_COORD
  call writeNetCDF_2d( 'x_rho', trim( GRID_FILE )     &
        , Nxr, Nyr, xr, start2D, count2D )
  call writeNetCDF_2d( 'y_rho', trim( GRID_FILE )     &
        , Nxr, Nyr, yr, start2D, count2D )
#endif

  count2D = (/ Nxp, Nyp /)
  call writeNetCDF_2d( 'lat_psi', trim( GRID_FILE )   &
        , Nxp, Nyp, latp, start2D, count2D )
  call writeNetCDF_2d( 'lon_psi', trim( GRID_FILE )   &
        , Nxp, Nyp, lonp, start2D, count2D )
  call writeNetCDF_2d( 'mask_psi', trim( GRID_FILE )  &
        , Nxp, Nyp, pmask, start2D, count2D )
#if defined UTM_COORD
  call writeNetCDF_2d( 'x_psi', trim( GRID_FILE )     &
        , Nxp, Nyp, xp, start2D, count2D )
  call writeNetCDF_2d( 'y_psi', trim( GRID_FILE )     &
        , Nxp, Nyp, yp, start2D, count2D )
#endif
  
  count2D = (/ Nxu, Nyu /)
  call writeNetCDF_2d( 'lat_u', trim( GRID_FILE )      &
        , Nxu, Nyu, latu, start2D, count2D )
  call writeNetCDF_2d( 'lon_u', trim( GRID_FILE )      &
        , Nxu, Nyu, lonu, start2D, count2D )
  call writeNetCDF_2d( 'mask_u', trim( GRID_FILE )     &
        , Nxu, Nyu, umask, start2D, count2D )
#if defined UTM_COORD
  call writeNetCDF_2d( 'x_u', trim( GRID_FILE )        &
        , Nxu, Nxu, xu, start2D, count2D )
  call writeNetCDF_2d( 'y_u', trim( GRID_FILE )        &
        , Nxu, Nyu, yu, start2D, count2D )
#endif
     
  count2D = (/ Nxv, Nyv /)
  call writeNetCDF_2d( 'lat_v', trim( GRID_FILE )      &
        , Nxv, Nyv, latv, start2D, count2D )
  call writeNetCDF_2d( 'lon_v', trim( GRID_FILE )      &
        , Nxv, Nyv, lonv, start2D, count2D )
  call writeNetCDF_2d( 'mask_v', trim( GRID_FILE )     &
        , Nxv, Nyv, vmask, start2D, count2D )
#if defined UTM_COORD
  call writeNetCDF_2d( 'x_v', trim( GRID_FILE )        &
        , Nxv, Nxv, xv, start2D, count2D )
  call writeNetCDF_2d( 'y_v', trim( GRID_FILE )        &
        , Nxv, Nyv, yv, start2D, count2D )
#endif

  write(*,*) 'FINISH!!'
      
END PROGRAM grdROMS
      
