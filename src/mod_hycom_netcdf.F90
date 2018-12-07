!!!=== Copyright (c) 2014-2018 Takashi NAKAMURA  =====

!!!**** ROMS netCDF MODULE ************************************

  MODULE mod_roms_netcdf
    
    use netcdf
    implicit none  
  
    CONTAINS

!**** create initial conditions NetCDF file **********************************************

      SUBROUTINE createNetCDFgrd(   &
!        input parameters
     &      OUT_FILE                &
     &    , Im, Jm                  &   
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm
      
      integer :: ncid,var_id
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_psi_dimid, eta_psi_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: dim2Dids(2)
      
!---- Create the ROMS grid netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid, 'xi_psi', Im-1, xi_psi_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid, 'eta_psi', Jm-1, eta_psi_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Jm-1, eta_v_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid, 'spherical', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid, var_id, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid, var_id, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid, 'xl', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'domain length in the XI-direction') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'el', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'domain length in the ETA-direction') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      dim2Dids = (/ xi_rho_dimid, eta_rho_dimid /)

      call check( nf90_def_var(ncid, 'h', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'Bathymetry at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'f', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'Coriolis parameter at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'second-1') )

      call check( nf90_def_var(ncid, 'pm', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'curvilinear coordinate metric in XI') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter-1') )

      call check( nf90_def_var(ncid, 'pn', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'curvilinear coordinate metric in ETA') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter-1') )

      call check( nf90_def_var(ncid, 'dndx', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'xi derivative of inverse metric factor pn') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'dmde', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'eta derivative of inverse metric factor pm') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'angle', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'angle between XI-axis and EAST') )
      call check( nf90_put_att(ncid, var_id, 'units',     'radians') )
      
      call check( nf90_def_var(ncid, 'x_rho', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'x location of RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'y_rho', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'y location of RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'lat_rho', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'latitude of RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_north') )

      call check( nf90_def_var(ncid, 'lon_rho', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'longitude of RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_east') )
      
      call check( nf90_def_var(ncid, 'mask_rho', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'mask on RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )


      dim2Dids = (/ xi_psi_dimid, eta_psi_dimid /)
      
      call check( nf90_def_var(ncid, 'x_psi', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'x location of PSI-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'y_psi', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'y location of PSI-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'lat_psi', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'latitude of PSI-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_north') )

      call check( nf90_def_var(ncid, 'lon_psi', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'longitude of PSI-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_east') )
      
      call check( nf90_def_var(ncid, 'mask_psi', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'mask on PSI-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      dim2Dids = (/ xi_u_dimid, eta_u_dimid /)
      
      call check( nf90_def_var(ncid, 'x_u', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'x location of U-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'y_u', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'y location of U-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'lat_u', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'latitude of U-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_north') )

      call check( nf90_def_var(ncid, 'lon_u', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'longitude of U-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_east') )
      
      call check( nf90_def_var(ncid, 'mask_u', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'mask on U-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      dim2Dids = (/ xi_v_dimid, eta_v_dimid /)
      
      call check( nf90_def_var(ncid, 'x_v', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'x location of V-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'y_v', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'y location of V-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      call check( nf90_def_var(ncid, 'lat_v', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'latitude of V-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_north') )

      call check( nf90_def_var(ncid, 'lon_v', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'longitude of V-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'degree_east') )
      
      call check( nf90_def_var(ncid, 'mask_v', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'mask on V-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

      dim2Dids = (/ xi_rho_dimid, eta_rho_dimid /)

      call check( nf90_def_var(ncid, 'p_coral', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'coral coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_def_var(ncid, 'p_coral2', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'coral2 coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_def_var(ncid, 'p_algae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'algal coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_def_var(ncid, 'p_seagrass', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'seagrass coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_def_var(ncid, 'p_sand', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'sand coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_def_var(ncid, 'sgd_src', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'source points of submarine groundwater discharge') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFgrd
      
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
      
      integer :: ncid,var_id
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: ocean_time_dimid
      integer :: dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ocean_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid, 'spherical', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid, var_id, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid, var_id, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid, 'Vtransform', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid, 'Vstretching', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid, 'theta_s', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid, 'theta_b', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid, 'Tcline', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid, 'hc', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )


      call check( nf90_def_var(ncid, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 's_w', NF90_DOUBLE, s_w_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid, 'ocean_time', NF90_DOUBLE, ocean_time_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      dim3Dids = (/ xi_rho_dimid, eta_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'zeta', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_u_dimid, eta_u_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'ubar', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated u-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim3Dids = (/ xi_v_dimid, eta_v_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'vbar', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated v-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_rho_dimid, eta_rho_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'temp', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'salt', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity') )
      call check( nf90_put_att(ncid, var_id, 'units',     'psu') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_u_dimid, eta_u_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'u', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xi_v_dimid, eta_v_dimid, s_rho_dimid, ocean_time_dimid /)

      call check( nf90_def_var(ncid, 'v', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
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
      
      integer :: ncid,var_id
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: xi_u_dimid, eta_u_dimid
      integer :: xi_v_dimid, eta_v_dimid
      integer :: s_rho_dimid, s_w_dimid
      integer :: bry_time_dimid
      integer :: dim2Dids(2), dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Im-1, xi_u_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Im, xi_v_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Jm, eta_u_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Jm-1, eta_v_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nz, s_rho_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nz+1, s_w_dimid) )
      call check( nf90_def_dim(ncid, 'bry_time', NF90_UNLIMITED, bry_time_dimid) )
      
    ! Define the netCDF variables.
      call check( nf90_def_var(ncid, 'spherical', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'grid type logical switch') )
      call check( nf90_put_att(ncid, var_id, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid, var_id, 'flag_meanings', 'Cartesian spherical' ) )
      
      call check( nf90_def_var(ncid, 'Vtransform', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertical terrain-following transformation equation') )

      call check( nf90_def_var(ncid, 'Vstretching', NF90_INT, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertical terrain-following stretching function') )

      call check( nf90_def_var(ncid, 'theta_s', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate surface control parameter') )

      call check( nf90_def_var(ncid, 'theta_b', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate bottom control parameter') )

      call check( nf90_def_var(ncid, 'Tcline', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate surface/bottom layer width') )

      call check( nf90_def_var(ncid, 'hc', NF90_DOUBLE, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate parameter, critical depth') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )


      call check( nf90_def_var(ncid, 's_rho', NF90_DOUBLE, s_rho_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 's_w', NF90_DOUBLE, s_w_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 'Cs_r', NF90_DOUBLE, s_rho_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid, 'Cs_w', NF90_DOUBLE, s_w_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )


      call check( nf90_def_var(ncid, 'bry_time', NF90_DOUBLE, bry_time_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'open boundary conditions time') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      dim2Dids = (/ eta_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'zeta_west', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'zeta_east', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim2Dids = (/ xi_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'zeta_south', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'zeta_north', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'ubar_west', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'ubar_east', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim2Dids = (/ xi_u_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'ubar_south', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'ubar_north', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      
      dim2Dids = (/ eta_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'vbar_west', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'vbar_east', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim2Dids = (/ xi_v_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'vbar_south', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'vbar_north', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )

      dim3Dids = (/ eta_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'u_west', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'u_east', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim3Dids = (/ xi_u_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'u_south', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'u_north', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )

      dim3Dids = (/ eta_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'v_west', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'v_east', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim3Dids = (/ xi_v_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'v_south', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'v_north', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'temp_west', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'temp_east', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'temp_south', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'temp_north', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )

      dim3Dids = (/ eta_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'salt_west', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity western boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'salt_east', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity eastern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      dim3Dids = (/ xi_rho_dimid, s_rho_dimid, bry_time_dimid /)
      call check( nf90_def_var(ncid, 'salt_south', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity southern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
      call check( nf90_def_var(ncid, 'salt_north', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity northern boundary condition') )
      call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )


  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
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
      
!
!
!**** readNetCDF_1d **********************************************
      
      SUBROUTINE readNetCDF_1d(    &
!    input parameters
     &    ncid, NCNAME, Im         &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im
      real(8), intent(out) :: data(Im)
      
      integer, parameter :: Num_try = 10
      integer :: var_id
      integer :: status
      integer :: itry
      
! --- Read NetCDF file ------------------------
      
      write(*,*) 'READ: ', NCNAME
      ! Get variable id
      do itry=1,Num_try
        status = nf90_inq_varid(ncid, NCNAME, var_id)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 1: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** READ FAILED: Retry!'
      end do        
      write(*,*) '*** SUCCESS'

    END SUBROUTINE readNetCDF_1d

!**** readNetCDF_2d **********************************************
      
      SUBROUTINE readNetCDF_2d(    &
!        input parameters
     &      ncid                   &
     &    , NCNAME                 &
     &    , Im, Jm                 &
     &    , start, count           &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm
      integer, intent( in) :: start(2), count(2)
      real(8), intent(out) :: data(Im, Jm)
      
      integer, parameter :: Num_try = 10
      integer :: var_id
      integer :: err_flag, status
      integer :: itry
      real(8) :: sf, off
      
! --- Read NetCDF file ------------------------
      
      write(*,*) 'READ: ', NCNAME
      ! Get variable id
      do itry=1,Num_try
        status = nf90_inq_varid(ncid, NCNAME, var_id)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 1: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data, start=start, count=count)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** READ FAILED: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'scale_factor', sf)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 2: Retry!'
      end do          
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'add_offset', off)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 3: Retry!'
      end do          
        
      data = data*sf+off
      write(*,*) '*** SUCCESS'

      END SUBROUTINE readNetCDF_2d

!**** readNetCDF_3d **********************************************
      
      SUBROUTINE readNetCDF_3d(    &
!        input parameters
     &      ncid                   &
     &    , NCNAME                 &
     &    , Im, Jm, Nt             &
     &    , start, count           &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nt
      integer, intent( in) :: start(3), count(3)
      real(8), intent(out) :: data(Im, Jm, Nt)
      
      integer, parameter :: Num_try = 10
      integer :: var_id
      integer :: err_flag, status
      integer :: itry
      real(8) :: sf, off
      
      data(:,:,:)=-9999.0d0
      
! --- Read NetCDF file ------------------------
      
      write(*,*) 'READ: ', NCNAME
      ! Get variable id
      do itry=1,Num_try
        status = nf90_inq_varid(ncid, NCNAME, var_id)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 1: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data, start=start, count=count)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr .and. data(Im,Jm,Nt)/=-9999.0d0) exit
        if (itry== Num_try) stop
        write(*,*) '*** READ FAILED: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'scale_factor', sf)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 2: Retry!'
      end do          
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'add_offset', off)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 3: Retry!'
      end do          
        
      data = data*sf+off
      write(*,*) '*** SUCCESS'

      END SUBROUTINE readNetCDF_3d
      
!**** readNetCDF_4d **********************************************
      
      SUBROUTINE readNetCDF_4d(    &
!        input parameters
     &      ncid                   &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start, count           &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start(4), count(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer, parameter :: Num_try = 10
      integer :: var_id
      integer :: err_flag, status
      integer :: itry
      real(8) :: sf, off
      
! --- Read NetCDF file ------------------------
      
      write(*,*) 'READ: ', NCNAME
      ! Get variable id
      do itry=1,Num_try
        status = nf90_inq_varid(ncid, NCNAME, var_id)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 1: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data, start=start, count=count)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** READ FAILED: Retry!'
      end do        
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'scale_factor', sf)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 2: Retry!'
      end do          
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'add_offset', off)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 3: Retry!'
      end do          
        
      data = data*sf+off
      write(*,*) '*** SUCCESS'

      END SUBROUTINE readNetCDF_4d
!
!**** readNetCDF_4d ver 2 **********************************************
      
      SUBROUTINE readNetCDF_4d_2(  &
!        input parameters
     &      ncid                   &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start, count           &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start(4), count(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer, parameter :: Num_try = 10
      integer :: start2(4), count2(4)
      integer :: var_id
      integer :: err_flag, status
      integer :: itry
      integer :: k
      real(8) :: sf, off
      
      start2 = start
      count2 = count
      count2(3)=1
      
      data(:,:,:,:) = -9999.0d0
      
! --- Read NetCDF file ------------------------
      
      write(*,*) 'READ: ', NCNAME
      ! Get variable id
      do itry=1,Num_try
        status = nf90_inq_varid(ncid, NCNAME, var_id)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 1: Retry!'
      end do        
      do k=1,Nz
        start2(3)=k
        do itry=1,Num_try
          status = nf90_get_var(ncid, var_id, data(:,:,k,:), start=start2, count=count2)
          write(*,*) trim(nf90_strerror(status)), k
          if (status == nf90_noerr .and. data(Im,Jm,k,Nt)/=-9999.0d0) exit
          if (itry== Num_try) stop
          write(*,*) '*** READ FAILED: Retry!'
        end do
      end do        
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'scale_factor', sf)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 2: Retry!'
      end do          
      do itry=1,Num_try
        status = nf90_get_att(ncid, var_id, 'add_offset', off)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** FAILED 3: Retry!'
      end do          
        
      data = data*sf+off
      write(*,*) '*** SUCCESS'

      END SUBROUTINE readNetCDF_4d_2
      
!**** NetCDF utility **********************************************
            
      SUBROUTINE try_nf_open(NC_FILE, nf90_open_mode, ncid)
      
      character(len=*),  intent( in) :: NC_FILE
      integer,           intent( in) :: nf90_open_mode
      integer,           intent(out) :: ncid

      integer, parameter :: Num_try = 10
      integer :: status
      integer :: itry
      
      do itry=1,Num_try
        status = nf90_open(NC_FILE, nf90_open_mode, ncid)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop
        write(*,*) '*** OPEN FAILED: Retry!'
      end do
      
      END SUBROUTINE try_nf_open
      
! -------------------------------------------------------------------------
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer, parameter :: Num_try = 10
      integer :: dimid
      integer :: status
      integer :: itry
      
      do itry=1,Num_try
        status = nf90_inq_dimid(ncid, name, dimid)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry== Num_try) stop 
      end do
      do itry=1,Num_try
        status = nf90_inquire_dimension(ncid, dimid, len=dim)
        write(*,*) trim(nf90_strerror(status))
        if (status == nf90_noerr) exit
        if (itry == Num_try) stop 
      end do
      
      END SUBROUTINE get_dimension


! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status
      
!      print *, trim(nf90_strerror(status))
      if (status /= nf90_noerr) then 
        write(*,*) trim(nf90_strerror(status))
        stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
        write(*,*) trim(nf90_strerror(status))
        err_flag = 1
!        stop "Stopped"
      end if
      
      END SUBROUTINE check2
      
  END MODULE mod_roms_netcdf
      
! -------------------------------------------------------------------------
