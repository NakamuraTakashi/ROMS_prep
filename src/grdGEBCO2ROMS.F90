
!!!=== Copyright (c) 2018 Takashi NAKAMURA  =====

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
! -------------------------------------------------------------------------

      real(8), allocatable :: lat_all(:), lon_all(:)
      real(8), allocatable :: lat(:), lon(:), depth(:,:)
      integer :: start1D(1), count1D(1)
      integer :: start2D(2), count2D(2)
      
      real(8), allocatable :: h(:,:)        ! depth (meter)
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

      real(8) :: d_lat, d_lon
      
      integer :: i,j
      
      integer :: N_lat_all, N_lon_all
      integer :: N_lat, N_lon
      integer :: id_slat, id_slon, id_elat, id_elon
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_xi_psi, N_eta_psi
      integer :: N_xi_u, N_eta_u
      integer :: N_xi_v, N_eta_v
      
      integer :: ncid,var_id
      
      namelist/grd/GRID_FILE
      namelist/bath/BATH_FILE
      namelist/bath/s_lat, e_lat, s_lon, e_lon
      namelist/bath/RESOLUTION
      namelist/bath/hmin
      namelist/bath/rx0max
      namelist/intpmode/mode
      namelist/hcoord/spherical
      
      ! Read parameters in namelist file
      
      read (*, nml=grd)
      read (*, nml=bath)
      read (*, nml=intpmode)
      read (*, nml=hcoord)

!---- Compute Lat, Lon fields fot ROMS grid netCDF file --------------------------------
      N_xi_rho  = int( (e_lon - s_lon)*RESOLUTION ) + 1
      N_eta_rho = int( (e_lat - s_lat)*RESOLUTION ) + 1
      N_xi_psi  = N_xi_rho-1
      N_eta_psi = N_eta_rho-1
      N_xi_u    = N_xi_rho-1
      N_eta_u   = N_eta_rho
      N_xi_v    = N_xi_rho
      N_eta_v   = N_eta_rho-1
!
      write(*,*) 'CHECK:',N_xi_rho, N_eta_rho
      
      allocate(h(N_xi_rho, N_eta_rho))
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

      d_lat = 1.0d0/RESOLUTION
      d_lon = 1.0d0/RESOLUTION

      do i=1, N_xi_rho
        lon_rho(i,:)=s_lon + d_lon*dble(i-1)
        lon_v(i,:)=s_lon + d_lon*dble(i-1)
      enddo
      do i=1, N_xi_rho-1
        lon_psi(i,:)=s_lon + d_lon*(dble(i-1)+0.5d0)
        lon_u(i,:)=s_lon + d_lon*(dble(i-1)+0.5d0)
      enddo
      do j=1, N_eta_rho
        lat_rho(:,j)=s_lat + d_lat*dble(j-1)
        lat_u(:,j)=s_lat + d_lat*dble(j-1)
      enddo
      do j=1, N_eta_rho-1
        lat_psi(:,j)=s_lat + d_lat*(dble(j-1)+0.5d0)
        lat_v(:,j)=s_lat + d_lat*(dble(j-1)+0.5d0)
      enddo
      
      ! Coriolis parameter
      do i=1, N_xi_rho
        do j=1, N_eta_rho
          f(i,j) = Coriolis(lat_rho(i,j))
        enddo
      enddo
      ! compute pm, pn
      do i=1, N_xi_rho
        do j=1, N_eta_rho
          pm(i,j) = 1.0d0/distance( lat_rho(i,j),                &
     &                              lon_rho(i,j)-0.5d0*d_lon,    &
     &                              lat_rho(i,j),                &
     &                              lon_rho(i,j)+0.5d0*d_lon )
          pn(i,j) = 1.0d0/distance( lat_rho(i,j)-0.5d0*d_lat,    &
     &                              lon_rho(i,j),                &
     &                              lat_rho(i,j)+0.5d0*d_lat,    &
     &                              lon_rho(i,j)             )
        enddo
      enddo
      
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
      CALL readNetCDF_2d(          &
!        input parameters
     &      ncid                   &
     &    , 'elevation'            &
     &    , N_lon, N_lat           &
     &    , start2D, count2D       &
!        output parameters
     &    , depth                  &
     & )
     
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
     &                  , depth , -5000.0d0, 8000.0d0                &
     &                  , N_xi_rho, N_eta_rho, lon_rho, lat_rho, h )

        write(*,*) '*** SUCCESS Linear Interporation'
        
      else if (mode ==2) then
      
      end if
      
      ! Compute land mask
      mask_rho(:,:) = 1.0d0
      mask_psi(:,:) = 1.0d0
      mask_u(:,:) = 1.0d0
      mask_v(:,:) = 1.0d0
      
      CALL land_masking(N_xi_rho, N_eta_rho, h, hmin, mask_rho)
      CALL uvp_masks(N_xi_rho, N_eta_rho, mask_rho, mask_u, mask_v, mask_psi)

      do i=1, N_xi_rho
        do j=1, N_eta_rho
          if(mask_rho(i,j) == 0.0d0) then
            h(i,j) = -10.0d0
          endif
        enddo
      enddo
      
      ! Smoothing Bathymetry
      write(*,*) 'Bathymetry smoothing'
      CALL bathy_smooth(N_xi_rho, N_eta_rho, rx0max, mask_rho, h)
      
!
!---- Create the ROMS grid netCDF file --------------------------------
          
      call createNetCDFgrd(                &
!          input parameters
     &        trim( GRID_FILE )                     &
     &      , N_xi_rho, N_eta_rho          &   
     &      )
      
!---- Write ROMS grid netCDF file --------------------------------

      call check( nf90_open(trim( GRID_FILE ), NF90_WRITE, ncid) )
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
      
