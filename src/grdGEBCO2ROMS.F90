
!!!=== Copyright (c) 2018-2019 Takashi NAKAMURA  =====

    PROGRAM grdGEBCO2ROMS
      use netcdf
      use mod_utility
      use mod_roms_netcdf
      use mod_interpolation
     
      implicit none
      
! -------------------------------------------------------------------------
      character(256) :: BATH_FILE
      character(256) :: GRID_FILE
      real(8) :: s_lat, e_lat, s_lon, e_lon
      real(8) :: RESOLUTION
      real(8) :: hmin
      real(8) :: rx0max
      integer :: mode
      integer :: spherical
      integer :: I_ocn, J_ocn
#if defined GRID_REFINEMENT
      character(256) :: parent_grid
      integer :: parent_Imin, parent_Imax
      integer :: parent_Jmin, parent_Jmax
      integer :: refine_factor
      integer, parameter :: Nmerge = 15
      real(8), parameter :: fmerge(Nmerge) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0  &
     &     , 1.0d0, 1.0d0, 0.9d0, 0.8d0, 0.7d0, 0.6d0, 0.5d0, 0.4d0, 0.3d0  &
     &     , 0.2d0, 0.1d0 /)
#endif
! -------------------------------------------------------------------------

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
      real(8), allocatable :: angle(:,:)
      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      real(8), allocatable :: mask_rho(:, :)
      real(8), allocatable :: lat_psi(:, :)
      real(8), allocatable :: lon_psi(:, :)
      real(8), allocatable :: mask_psi(:, :)
      real(8), allocatable :: lat_u(:, :)
      real(8), allocatable :: lon_u(:, :)
      real(8), allocatable :: mask_u(:, :)
      real(8), allocatable :: lat_v(:, :)
      real(8), allocatable :: lon_v(:, :)
      real(8), allocatable :: mask_v(:, :)
#if defined GRID_REFINEMENT
      real(8), allocatable :: hr(:,:)        ! depth (meter)
      real(8), allocatable :: mask_rho_r(:,:)        ! depth (meter)
      real(8), allocatable :: xr(:,:) 
      real(8), allocatable :: yr(:,:)
       
      real(8), allocatable :: Xd(:)
      real(8), allocatable :: Yd(:)
      real(8), allocatable :: P_h(:,:)        ! depth (meter)
      real(8), allocatable :: P_lat_rho(:, :)
      real(8), allocatable :: P_lon_rho(:, :)
      real(8), allocatable :: P_mask_rho(:, :)
      real(8), allocatable :: P_lat_psi(:, :)
      real(8), allocatable :: P_lon_psi(:, :)
      real(8), allocatable :: P_lat_u(:, :)
      real(8), allocatable :: P_lon_u(:, :)
      real(8), allocatable :: P_lat_v(:, :)
      real(8), allocatable :: P_lon_v(:, :)
#endif

      real(8) :: d_lat, d_lon
      
      integer :: i,j
      
      integer :: N_lat_all, N_lon_all
      integer :: N_lat, N_lon
      integer :: id_slat, id_slon, id_elat, id_elon
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_xi_psi, N_eta_psi
      integer :: N_xi_u, N_eta_u
      integer :: N_xi_v, N_eta_v
#if defined GRID_REFINEMENT
      integer :: P_N_xi_rho, P_N_eta_rho
      integer :: P_N_xi_psi, P_N_eta_psi
      integer :: P_N_xi_u, P_N_eta_u
      integer :: P_N_xi_v, P_N_eta_v
      integer :: ip,jp,is,js
#endif
     
      integer :: ncid,var_id
      
      namelist/grd/GRID_FILE
      namelist/bath/BATH_FILE
      namelist/bath/s_lat, e_lat, s_lon, e_lon
      namelist/bath/RESOLUTION
      namelist/bath/hmin
      namelist/bath/rx0max
      namelist/bath/I_ocn, J_ocn
#if defined GRID_REFINEMENT
      namelist/refinement/parent_grid
      namelist/refinement/parent_Imin, parent_Imax
      namelist/refinement/parent_Jmin, parent_Jmax
      namelist/refinement/refine_factor
#endif
      namelist/intpmode/mode
      namelist/hcoord/spherical
      
      ! Read parameters in namelist file
      
      read (*, nml=grd)
      read (*, nml=bath)
#if defined GRID_REFINEMENT
      read (*, nml=refinement)
#endif
      read (*, nml=intpmode)
      read (*, nml=hcoord)

!---- Compute Lat, Lon fields for ROMS grid netCDF file --------------------------------
#if defined GRID_REFINEMENT
      N_xi_rho  = (parent_Imax-parent_Imin)*refine_factor+2
      N_eta_rho = (parent_Jmax-parent_Jmin)*refine_factor+2 
#else      
      N_xi_rho  = int( (e_lon - s_lon)*RESOLUTION ) + 1
      N_eta_rho = int( (e_lat - s_lat)*RESOLUTION ) + 1
#endif
      N_xi_psi  = N_xi_rho-1
      N_eta_psi = N_eta_rho-1
      N_xi_u    = N_xi_rho-1
      N_eta_u   = N_eta_rho
      N_xi_v    = N_xi_rho
      N_eta_v   = N_eta_rho-1
!
      write(*,*) 'CHECK:',N_xi_rho, N_eta_rho
      
      allocate(h(N_xi_rho, N_eta_rho))
      allocate(hraw(N_xi_rho, N_eta_rho))
      allocate(f(N_xi_rho, N_eta_rho))
      allocate(pm(N_xi_rho, N_eta_rho))
      allocate(pn(N_xi_rho, N_eta_rho))
      allocate(dndx(N_xi_rho, N_eta_rho))
      allocate(dmde(N_xi_rho, N_eta_rho))
      allocate(angle(N_xi_rho, N_eta_rho))
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(mask_rho(N_xi_rho, N_eta_rho))
      allocate(lat_psi(N_xi_psi, N_eta_psi))
      allocate(lon_psi(N_xi_psi, N_eta_psi))
      allocate(mask_psi(N_xi_psi, N_eta_psi))
      allocate(lat_u(N_xi_u, N_eta_u))
      allocate(lon_u(N_xi_u, N_eta_u))
      allocate(mask_u(N_xi_u, N_eta_u))
      allocate(lat_v(N_xi_v, N_eta_v))
      allocate(lon_v(N_xi_v, N_eta_v))
      allocate(mask_v(N_xi_v, N_eta_v))
#if defined GRID_REFINEMENT
      allocate(hr(N_xi_rho, N_eta_rho))
      allocate(mask_rho_r(N_xi_rho, N_eta_rho))
      allocate(xr(N_xi_rho, N_eta_rho))
      allocate(yr(N_xi_rho, N_eta_rho))
#endif

#if defined GRID_REFINEMENT
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", trim(parent_grid)
      
      ! Open NetCDF grid file
      call check( nf90_open( trim(parent_grid), nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  P_N_xi_rho)
      call get_dimension(ncid, 'eta_rho', P_N_eta_rho)
      call get_dimension(ncid, 'xi_u',  P_N_xi_u)
      call get_dimension(ncid, 'eta_u', P_N_eta_u)
      call get_dimension(ncid, 'xi_v',  P_N_xi_v)
      call get_dimension(ncid, 'eta_v', P_N_eta_v)
      
      allocate(P_lat_rho(P_N_xi_rho, P_N_eta_rho))
      allocate(P_lon_rho(P_N_xi_rho, P_N_eta_rho))
      allocate(P_lat_u(P_N_xi_u, P_N_eta_u))
      allocate(P_lon_u(P_N_xi_u, P_N_eta_u))
      allocate(P_lat_v(P_N_xi_v, P_N_eta_v))
      allocate(P_lon_v(P_N_xi_v, P_N_eta_v))
      allocate(P_h(P_N_xi_rho, P_N_eta_rho))
      allocate(P_mask_rho(P_N_xi_rho, P_N_eta_rho))
      allocate(Xd(P_N_xi_rho), Yd(P_N_eta_rho))
            
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, P_lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, P_lon_rho) )
!      call check( nf90_inq_varid(ncid, 'lat_u', var_id) ) ! latitude at U-points (degree_east)
!      call check( nf90_get_var(ncid, var_id, P_lat_u) )
!      call check( nf90_inq_varid(ncid, 'lon_u', var_id) ) ! longitude at U-points (degree_east)
!      call check( nf90_get_var(ncid, var_id, P_lon_u) )
!      call check( nf90_inq_varid(ncid, 'lat_v', var_id) ) ! latitude at V-points (degree_east)
!      call check( nf90_get_var(ncid, var_id, P_lat_v) )
!      call check( nf90_inq_varid(ncid, 'lon_v', var_id) ) ! longitude at V-points (degree_east)
!      call check( nf90_get_var(ncid, var_id, P_lon_v) )
      call check( nf90_inq_varid(ncid, 'h', var_id) ) 
      call check( nf90_get_var(ncid, var_id, P_h) )
      call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
      call check( nf90_get_var(ncid, var_id, P_mask_rho) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", trim(parent_grid)

      do i=1, P_N_xi_rho
        Xd(i)=dble(i-1)
      enddo
      do j=1, P_N_eta_rho
        Yd(j)=dble(j-1)
      enddo
      do i=1, N_xi_rho
        xr(i,:)=dble(parent_Imin-1) + dble(i-1+refine_factor/2)/dble(refine_factor)
      enddo
      do j=1, N_eta_rho
        yr(:,j)=dble(parent_Jmin-1) + dble(j-1+refine_factor/2)/dble(refine_factor)
      enddo

      call LinearInterpolation2D_grid2(P_N_xi_rho, P_N_eta_rho, Xd, Yd    &
     &                  , P_lat_rho , -10000.0d0, 10000.0d0               &
     &                  , N_xi_rho, N_eta_rho, xr, yr, lat_rho )
      call LinearInterpolation2D_grid2(P_N_xi_rho, P_N_eta_rho, Xd, Yd    &
     &                  , P_lon_rho , -10000.0d0, 10000.0d0               &
     &                  , N_xi_rho, N_eta_rho, xr, yr, lon_rho )

      call LinearInterpolation2D_grid2(P_N_xi_rho, P_N_eta_rho, Xd, Yd    &
     &                  , P_h , -9.0d0, 10000.0d0                         &
     &                  , N_xi_rho, N_eta_rho, xr, yr, hr )

     do ip=parent_Imin+1, parent_Imax
      do jp=parent_Jmin+1, parent_Jmax
        is = (ip-parent_Imin-1)*refine_factor+2
        js = (jp-parent_Jmin-1)*refine_factor+2
        do i=is, is+refine_factor
          do j=js, js+refine_factor
            mask_rho_r(i,j) = P_mask_rho(ip,jp)
          enddo
        enddo
      enddo
    enddo
    mask_rho_r(1,:)=mask_rho_r(2,:)
    mask_rho_r(N_xi_rho,:)=mask_rho_r(N_xi_rho-1,:)
    mask_rho_r(:,1)=mask_rho_r(:,2)
    mask_rho_r(:,N_eta_rho)=mask_rho_r(:,N_eta_rho-1)

#else      

      d_lat = 1.0d0/RESOLUTION
      d_lon = 1.0d0/RESOLUTION

      do i=1, N_xi_rho
        lon_rho(i,:)=s_lon + d_lon*dble(i-1)
      enddo
      do j=1, N_eta_rho
        lat_rho(:,j)=s_lat + d_lat*dble(j-1)
      enddo
#endif

      do i=1, N_xi_u
        do j=1, N_eta_u
          lon_u(i,j)=0.5d0*(lon_rho(i,j)+lon_rho(i+1,j))
          lat_u(i,j)=0.5d0*(lat_rho(i,j)+lat_rho(i+1,j))
        enddo      
      enddo
      do i=1, N_xi_v
        do j=1, N_eta_v
          lon_v(i,j)=0.5d0*(lon_rho(i,j)+lon_rho(i,j+1))
          lat_v(i,j)=0.5d0*(lat_rho(i,j)+lat_rho(i,j+1))
        enddo      
      enddo
      do i=1, N_xi_psi
        do j=1, N_eta_psi
          lon_psi(i,j)=0.5d0*(lon_u(i,j)+lon_u(i,j+1))
          lat_psi(i,j)=0.5d0*(lat_v(i,j)+lat_v(i+1,j))
        enddo      
      enddo
      ! Coriolis parameter
      do i=1, N_xi_rho
        do j=1, N_eta_rho
          f(i,j) = Coriolis(lat_rho(i,j))
        enddo
      enddo
      ! compute pm, pn
!      do i=1, N_xi_rho
!        do j=1, N_eta_rho
!          pm(i,j) = 1.0d0/distance( lat_rho(i,j),                &
!     &                              lon_rho(i,j)-0.5d0*d_lon,    &
!     &                              lat_rho(i,j),                &
!     &                              lon_rho(i,j)+0.5d0*d_lon )
!          pn(i,j) = 1.0d0/distance( lat_rho(i,j)-0.5d0*d_lat,    &
!     &                              lon_rho(i,j),                &
!     &                              lat_rho(i,j)+0.5d0*d_lat,    &
!     &                              lon_rho(i,j)             )
!        enddo
!      enddo
      do i=2, N_xi_rho-1
        do j=1, N_eta_rho
          pm(i,j) = 1.0d0/distance( lat_rho(i,j), lon_u(i-1,j)   &
     &                            , lat_rho(i,j), lon_u(i,j))
        enddo
      enddo
      pm(1,:) = pm(2,:) 
      pm(N_xi_rho,:) = pm(N_xi_rho-1,:) 

      do i=1, N_xi_rho
        do j=2, N_eta_rho-1
          pn(i,j) = 1.0d0/distance( lat_v(i,j-1), lon_rho(i,j)   &
     &                            , lat_v(i,j), lon_rho(i,j))
        enddo
      enddo
      pn(:,1) = pn(:,2) 
      pn(:,N_eta_rho) = pn(:,N_eta_rho-1) 
           
      dndx(:,:) = 0.0d0
      dmde(:,:) = 0.0d0
      angle(:,:) = 0.0d0

!---- Read GEBCO grid netCDF file --------------------------------
      write(*,*) "OPEN: ", trim( BATH_FILE )
      
      ! Open NetCDF grid file
      call check( nf90_open(trim( BATH_FILE ), nf90_nowrite, ncid) )
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
      
      ! Trim GEBCO bathymetry
#if defined GRID_REFINEMENT
      s_lat = lat_rho(1,1)
      e_lat = lat_rho(N_xi_rho,N_eta_rho)
      s_lon = lon_rho(1,1)
      e_lon = lon_rho(N_xi_rho,N_eta_rho)
#endif

      do i=2, N_lat_all
        if(s_lat < lat_all(i)) then
          id_slat = i-1
          exit
        endif
      enddo
      do i=1, N_lat_all
        if(e_lat < lat_all(i)) then
          id_elat = i
          exit
        endif
      enddo
      do i=2, N_lon_all
        if(s_lon < lon_all(i)) then
          id_slon = i-1
          exit
        endif
      enddo
      do i=1, N_lon_all
        if(e_lon < lon_all(i)) then
          id_elon = i
          exit
        endif
      enddo
      
      N_lon = id_elon-id_slon+1
      N_lat = id_elat-id_slat+1
      allocate(depth(N_lon,N_lat))
      start2D = (/ id_slon, id_slat /)
      count2D = (/ N_lon  , N_lat   /)
      write(*,*) 'CHECK:',id_slon,id_elon, id_slat,id_elat, N_lon  , N_lat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'CHECK:',lon_all(id_slon),lon_all(id_elon), lat_all(id_slat),lat_all(id_elat)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'elevation', var_id) )
      call check( nf90_get_var(ncid, var_id, depth, start=start2D, count=count2D) )
     
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
      write(*,*) 'CHECK:',lon_rho(1,1),lon_rho(N_xi_rho, N_eta_rho),lat_rho(1,1),lat_rho(N_xi_rho, N_eta_rho)
      
! --- Linear interporation --------------------------------
      
      if (mode==1) then
      
        write(*,*) 'Linear Interporation: h'
        call LinearInterpolation2D_grid2(N_lon, N_lat, lon, lat      &
     &                  , depth , -10000.0d0, 10000.0d0              &
     &                  , N_xi_rho, N_eta_rho, lon_rho, lat_rho, h )

        write(*,*) '*** SUCCESS Linear Interporation'
        
      else if (mode ==2) then
      
      end if

      hraw(:,:) = h(:,:)
      
      ! Compute land mask
      mask_rho(:,:) = 1.0d0
      mask_psi(:,:) = 1.0d0
      mask_u(:,:) = 1.0d0
      mask_v(:,:) = 1.0d0
      
      CALL land_masking(N_xi_rho, N_eta_rho, h, hmin, mask_rho)

      ! Smoothing Bathymetry
      CALL bathy_smooth(N_xi_rho, N_eta_rho, rx0max, mask_rho, h)
    
#if defined GRID_REFINEMENT
      ! Merge GEBCO grid with refined ROMS grid
      do j=1, Nmerge
        do i=1+j-1, N_xi_rho-j+1
          h(i,j) = fmerge(j)*hr(i,j) + (1.0d0-fmerge(j))*h(i,j)
          h(i,N_eta_rho-j+1) = fmerge(j)*hr(i,N_eta_rho-j+1) + (1.0d0-fmerge(j))*h(i,N_eta_rho-j+1)
        enddo
      enddo
      do i=1, Nmerge
        do j=1+i-1, N_eta_rho-i+1
          h(i,j) = fmerge(i)*hr(i,j) + (1.0d0-fmerge(i))*h(i,j)
          h(N_xi_rho-i+1,j) = fmerge(i)*hr(N_xi_rho-i+1,j) + (1.0d0-fmerge(i))*h(N_xi_rho-i+1,j)
        enddo
      enddo
      do j=1, refine_factor/2+1
        mask_rho(:,j) = mask_rho_r(:,j)
        mask_rho(:,N_eta_rho-j+1) = mask_rho_r(:,N_eta_rho-j+1)
      enddo
      do i=1, refine_factor/2+1
        mask_rho(i,:) = mask_rho_r(i,:)
        mask_rho(N_xi_rho-i+1,:) = mask_rho_r(N_xi_rho-i+1,:)
      enddo
      
#endif

      CALL isolated_water_masking(N_xi_rho, N_eta_rho, 50, 50, mask_rho)

      CALL uvp_masks(N_xi_rho, N_eta_rho, mask_rho, mask_u, mask_v, mask_psi)

      do i=1, N_xi_rho
        do j=1, N_eta_rho
          if(mask_rho(i,j) == 0.0d0) then
            h(i,j) = -10.0d0
          endif
          if(h(i,j) == 0.0d0) then
            h(i,j) = 0.001d0
          endif
        enddo
      enddo
      
!
!---- Create the ROMS grid netCDF file --------------------------------
          
      call createNetCDFgrd(                &
!          input parameters
     &        trim( GRID_FILE )                     &
     &      , N_xi_rho, N_eta_rho          &   
     &      )
      
!---- Write ROMS grid netCDF file --------------------------------

      call check( nf90_open(trim( GRID_FILE ), NF90_WRITE, ncid) )
!     Put global attribute
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
      count2D = (/ N_xi_rho, N_eta_rho /)
      call writeNetCDF_2d(                      &
!          input parameters
     &        'h'                               &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , h                                 &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'f'                               &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , f                                 &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'pm'                              &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , pm                                &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'pn'                              &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , pn                                &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'dndx'                            &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , dndx                              &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'dmde'                            &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , dmde                              &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'angle'                           &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , angle                             &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lat_rho'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , lat_rho                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lon_rho'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , lon_rho                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'mask_rho'                        &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , mask_rho                          &
     &      , start2D, count2D                  &
     &      )
     
      call writeNetCDF_2d(                      &
!          input parameters
     &        'hraw'                            &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , hraw                              &
     &      , start2D, count2D                  &
     &      )

#if defined GRID_REFINEMENT
      call writeNetCDF_2d(                      &
!          input parameters
     &        'mask_rho_coarse'                 &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_rho, N_eta_rho               &
     &      , mask_rho_r                        &
     &      , start2D, count2D                  &
     &      )
#endif
      
      count2D = (/ N_xi_psi, N_eta_psi /)
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lat_psi'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , lat_psi                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lon_psi'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , lon_psi                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'mask_psi'                        &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , mask_psi                          &
     &      , start2D, count2D                  &
     &      )
      
      count2D = (/ N_xi_psi, N_eta_psi /)
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lat_psi'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , lat_psi                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lon_psi'                         &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , lon_psi                           &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'mask_psi'                        &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_psi, N_eta_psi               &
     &      , mask_psi                          &
     &      , start2D, count2D                  &
     &      )
      
      count2D = (/ N_xi_u, N_eta_u /)
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lat_u'                           &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_u, N_eta_u                   &
     &      , lat_u                             &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lon_u'                           &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_u, N_eta_u                   &
     &      , lon_u                             &
     &      , start2D, count2D                  &
     &      )
     call writeNetCDF_2d(                       &
!          input parameters
     &        'mask_u'                          &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_u, N_eta_u                   &
     &      , mask_u                            &
     &      , start2D, count2D                  &
     &      )
         
      count2D = (/ N_xi_v, N_eta_v /)
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lat_v'                           &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_v, N_eta_v                   &
     &      , lat_v                             &
     &      , start2D, count2D                  &
     &      )
      call writeNetCDF_2d(                      &
!          input parameters
     &        'lon_v'                           &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_v, N_eta_v                   &
     &      , lon_v                             &
     &      , start2D, count2D                  &
     &      )
     call writeNetCDF_2d(                       &
!          input parameters
     &        'mask_v'                          &
     &      , trim( GRID_FILE )                 &
     &      , N_xi_v, N_eta_v                   &
     &      , mask_v                            &
     &      , start2D, count2D                  &
     &      )

      write(*,*) 'FINISH!!'
      
    END PROGRAM grdGEBCO2ROMS
      
