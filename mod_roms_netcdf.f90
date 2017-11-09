
!!!=== Copyright (c) 2014-2017 Takashi NAKAMURA  =====

!!!**** ROMS netCDF MODULE ************************************

  MODULE mod_roms_netcdf
  
    use netcdf
  
    CONTAINS
    
!**** create initial conditions NetCDF file **********************************************

      SUBROUTINE createNetCDFini(   &
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
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: ocean_time_dimid
      integer :: dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid2, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid2, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid2, 'ocean_time', NF90_UNLIMITED, ocean_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid2, 'spherical', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid2, 'Vtransform', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid2, 'Vstretching', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid2, 'theta_s', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid2, 'theta_b', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid2, 'Tcline', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid2, 'hc', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )


      call check( nf90_def_var(ncid2, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 's_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid2, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid2, 'ocean_time', NF90_DOUBLE, ocean_time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      dim3Dids = (/ xi_rho_dimid, eta_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'zeta', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_u_dimid, eta_u_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'ubar', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertically integrated u-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_v_dimid, eta_v_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'vbar', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertically integrated v-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_rho_dimid, eta_rho_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'temp', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid2, 'salt', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'psu') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_u_dimid, eta_u_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'u', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_v_dimid, eta_v_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid2, 'v', NF90_DOUBLE, dim4Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum component') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'ocean_time') )

  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFini
      
!**** create boundary conditions NetCDF file **********************************************

      SUBROUTINE createNetCDFbry(   &
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
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: bry_time_dimid
      integer :: dim2Dids(2), dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid2) )

      call check( nf90_def_dim(ncid2, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid2, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid2, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid2, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid2, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid2, 'bry_time', NF90_UNLIMITED, bry_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid2, 'spherical', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid2, var_id2, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid2, 'Vtransform', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid2, 'Vstretching', NF90_INT, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid2, 'theta_s', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid2, 'theta_b', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid2, 'Tcline', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid2, 'hc', NF90_DOUBLE, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )


      call check( nf90_def_var(ncid2, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 's_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'positive',  'up' ) )
      call check( nf90_put_att(ncid2, var_id2, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid2, var_id2, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid2, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid2, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid2, var_id2, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid2, var_id2, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid2, 'bry_time', NF90_DOUBLE, bry_time_dimid, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'open boundary conditions time') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )

      dim2Dids = (/ eta_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'zeta_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'zeta_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'zeta_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'zeta_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'free-surface northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'ubar_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'ubar_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'ubar_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'ubar_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'vbar_west', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'vbar_east', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim2Dids = (/ xi_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'vbar_south', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'vbar_north', NF90_DOUBLE, dim2Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', '2D v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'u_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'u_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'u_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'u_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'v_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'v_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'v_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'v_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'temp_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'temp_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'temp_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'temp_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'potential temperature northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'units',     'Celsius') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'salt_west', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity western boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'salt_east', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity eastern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid2, 'salt_south', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity southern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )
      call check( nf90_def_var(ncid2, 'salt_north', NF90_DOUBLE, dim3Dids, var_id2) )
      call check( nf90_put_att(ncid2, var_id2, 'long_name', 'salinity northern boundary condition') )
      call check( nf90_put_att(ncid2, var_id2, 'time',      'bry_time') )


  ! End define mode.
      call check( nf90_enddef(ncid2) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFbry
      
!**** writeNetCDF_1d **********************************************
      
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
      
      write(*,*) "WRITE ", NCNAME," to ", trim( OUT_FILE )
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start1D, count = count1D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_1d
      
!**** writeNetCDF_2d **********************************************
      
      SUBROUTINE writeNetCDF_2d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm                 &
     &    , data                   &
     &    , start2D, count2D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm
      real(8), intent( in) :: data(Im, Jm)
      integer, intent( in) :: start2D(2), count2D(2)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", trim( OUT_FILE )
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start2D, count = count2D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_2d
      
!**** writeNetCDF_3d **********************************************
      
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
      
      write(*,*) "WRITE ", NCNAME," to ", trim( OUT_FILE )
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_3d
      
!**** writeNetCDF_4d **********************************************
      
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
      
      write(*,*) "WRITE ", NCNAME," to ", trim( OUT_FILE )
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start4D, count = count4D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_4d
      
!**** readNetCDF_3d **********************************************
      
      SUBROUTINE readNetCDF_3d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nt             &
     &    , start3D, count3D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nt
      integer, intent( in) :: start3D(3), count3D(3)
      real(8), intent(out) :: data(Im, Jm, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start3D, count=count3D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do

      data(:,:,:)=data(:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_3d
      
!**** readNetCDF_4d **********************************************
      
      SUBROUTINE readNetCDF_4d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start4D, count4D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start4D(4), count4D(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start4D, count=count4D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do
      
      data(:,:,:,:)=data(:,:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_4d
      
!**** NetCDF utility **********************************************
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension

! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          err_flag = 1
!          stop "Stopped"
      end if
      
      END SUBROUTINE check2
      
  END MODULE mod_roms_netcdf
      
! -------------------------------------------------------------------------
