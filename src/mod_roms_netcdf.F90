
!!!=== Copyright (c) 2014-2021 Takashi NAKAMURA  =====

!!!**** ROMS netCDF MODULE ************************************

MODULE mod_roms_netcdf
  
  use netcdf
  implicit none  

  integer, parameter :: N_var = 36

  integer, parameter :: Ndom = 2     ! Number of Dissolved Organic Matter types
  integer, parameter :: Npom  = 2    ! Number of Particulate Organic Matter types
  integer, parameter :: Nphy  = 3    ! Number of Phytoplankton types
  integer, parameter :: Nzoo = 1     ! Number of Zooplankton types
  integer, parameter :: Npim   = 1   ! Number of Particulate Inorganic Matter types

  character(256), parameter :: VAR_NAME(N_var) = (/ &
     "zeta          " &  !  1
    ,"u             " &  !  2 
    ,"v             " &  !  3
    ,"ubar          " &  !  4 
    ,"vbar          " &  !  5 
    ,"temp          " &  !  6 
    ,"salt          " &  !  7 
    ,"mud_          " &  !  8 
    ,"sand_         " &  !  9
    ,"TIC           " &  ! 10
    ,"alkalinity    " &  ! 11
    ,"oxygen        " &  ! 12
    ,"DOC_          " &  ! 13
    ,"POC_          " &  ! 14
    ,"phytoplankton_" &  ! 15
    ,"zooplankton_  " &  ! 16
    ,"PIC_          " &  ! 17
    ,"NO3           " &  ! 18
    ,"NH4           " &  ! 19
    ,"PO4           " &  ! 20
    ,"DON_          " &  ! 21
    ,"PON_          " &  ! 22
    ,"DOP_          " &  ! 23
    ,"POP_          " &  ! 24
    ,"TI13C         " &  ! 25
    ,"DO13C_        " &  ! 26
    ,"PO13C_        " &  ! 27
    ,"phyt13C_      " &  ! 28
    ,"zoop13C_      " &  ! 29
    ,"PI13C_        " &  ! 30
    ,"15NO3         " &  ! 31
    ,"15NH4         " &  ! 32
    ,"DO15N_        " &  ! 33
    ,"PO15N_        " &  ! 34
    ,"phyt15N_      " &  ! 35
    ,"zoop15N_      " &  ! 36
     /)
  character(256), parameter :: VAR_LONG_NAME(N_var) = (/ &
     "free-surface                                " &
    ,"u-momentum component                        " &
    ,"v-momentum component                        " &
    ,"vertically integrated u-momentum component  " &
    ,"vertically integrated v-momentum component  " &
    ,"potential temperature                       " &
    ,"salinity                                    " &
    ,"suspended cohesive sediment                 " &
    ,"suspended noncohesive sediment              " &
    ,"total inorganic carbon                      " &
    ,"total alkalinity                            " &
    ,"dissolved oxygen concentration              " &
    ,"dissolved organic carbon                    " &
    ,"particulate organic carbon                  " &
    ,"phytoplankton                               " &
    ,"zooplankton                                 " &
    ,"particulate inorganic carbon                " &
    ,"nitrate concentration                       " &
    ,"ammonium concentration                      " &
    ,"phosphate concentration                     " &
    ,"dissolved organic nitrogen concentration    " &
    ,"particulate organic nitrogen concentration  " &
    ,"dissolved organic phosphorus concentration  " &
    ,"particulate organic phosphorus concentration" &
    ,"carbon 13 of total inorganic carbon         " &
    ,"carbon 13 of DOC                            " &
    ,"carbon 13 of POC                            " &
    ,"carbon 13 of phytoplankton                  " &
    ,"carbon 13 of zooplankton                    " &
    ,"carbon 13 of PIC                            " &
    ,"nitrogen 15 of nitrate                      " &
    ,"nitrogen 15 of ammonium                     " &
    ,"nitrogen 15 of DON                          " &
    ,"nitrogen 15 of PON                          " &
    ,"nitrogen 15 of phytoplankton                " &
    ,"nitrogen 15 of zooplankton                  " &
    /)
  character(256), parameter :: VAR_UNIT(N_var) = (/ &
     "meter           " &
    ,"meter second-1  " &
    ,"meter second-1  " &
    ,"meter second-1  " &
    ,"meter second-1  " &
    ,"Celsius         " &
    ,"nondimensional  " &
    ,"kilogram meter-3" &
    ,"kilogram meter-3" &
    ,"umol kg-1       " &
    ,"umol kg-1       " &
    ,"oxygen          " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    ,"umol L-1        " &
    /)
  character(256), parameter :: BRY_NAME(4) =  &
    (/ "south", "north", "west ", "east " /)
  character(256), parameter :: BRY_LONGNAME(4) = (/ &
      "southern boundary condition" &
    , "northern boundary condition" &
    , "western boundary condition " &
    , "eastern boundary condition " &
    /)

! ----- HYCOM OpenDAP URL -----

#if defined GOFS_31
# if defined ANALYSIS_Y
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2018-12-04 to present) -----
!        ** GLBy0.08 grid is 0.08 deg lon x 0.04 deg lat that covers 80 S to 90 N.
  integer, parameter :: NCnum   = 1
  character(54) :: HYCOM_FILE(NCnum) = (/                           &
       "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0"  &  ! lon(4500): 0 - 359.92
    /)
# elif defined ANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2014-07-01 to 2020-02-18) -----
!        ** grid is 0.08 deg lon x 0.08 deg lat between 40S-40N. 
!           Poleward of 40S/40N, the grid is 0.08 deg lon x 0.04 deg lat. It spans 80S to 90N.
  integer, parameter :: NCnum   = 6
  character(53) :: HYCOM_FILE(NCnum) = (/                          &
       "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3"  &  ! lon(4500): -180 - 179.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.2"  &  ! lon(4500): -180 - 179.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.8"  &  ! lon(4500): 0 - 359.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7"  &  ! lon(4500): -180 - 179.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9"  &  ! lon(4500): 0 - 359.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0"  &  ! lon(4500): 0 - 359.92
    /)
# elif defined REANALYSIS
! ----- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Renalysis (since 1994-01-01 to 2015-12-31) -----
  integer, parameter :: NCnum   = 22
  character(63) :: HYCOM_FILE(NCnum) = (/                                    &
       "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1994"  &  ! lon(4500): -180 - 179.92
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1995"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1996"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1997"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1998"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/1999"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2000"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2001"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2002"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2003"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2004"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2005"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2006"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2007"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2008"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2009"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2010"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2011"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2012"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2013"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2014"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2015"  &  ! lon(4500): -180 - 179.92
    /)
# endif
#elif defined GOFS_30
# if defined ANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to 2018-11-20) -----
  integer, parameter :: NCnum   = 4
  character(53) :: HYCOM_FILE(NCnum) = (/                          &
       "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2"  &
    /)
# elif defined REANALYSIS
! -----  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31) -----
  integer, parameter :: NCnum   = 2
  character(53) :: HYCOM_FILE(NCnum) = (/                          &
       "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0"  &
      ,"http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"  &
    /)
# endif
#endif

  CONTAINS

!**** create initial conditions NetCDF file **********************************************

    SUBROUTINE createNetCDFgrd(   &
!      input parameters
          OUT_FILE                &
        , Nxr, Nyr                &   
    )
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      integer, intent( in) :: Nxr, Nyr
      
      integer :: ncid,var_id
      integer :: xr_dimid, yr_dimid
      integer :: xi_psi_dimid, eta_psi_dimid
      integer :: xu_dimid, yu_dimid
      integer :: xv_dimid, yv_dimid
      integer :: dim2Dids(2)
      
!---- Create the ROMS grid netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xr_dimid) )
      call check( nf90_def_dim(ncid, 'xi_psi', Nxr-1, xi_psi_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Nxr-1, xu_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Nxr, xv_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Nyr, yr_dimid) )
      call check( nf90_def_dim(ncid, 'eta_psi', Nyr-1, eta_psi_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Nyr, yu_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Nyr-1, yv_dimid) )
      
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

      dim2Dids = (/ xr_dimid, yr_dimid /)

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

      call check( nf90_def_var(ncid, 'hraw', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'Raw bathymetric data at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )

#if defined GRID_REFINEMENT
      call check( nf90_def_var(ncid, 'mask_rho_coarse', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'mask on RHO-points of coarse grid') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
#endif

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

      dim2Dids = (/ xu_dimid, yu_dimid /)
      
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

      dim2Dids = (/ xv_dimid, yv_dimid /)
      
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

      dim2Dids = (/ xr_dimid, yr_dimid /)

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

!************************************************
! Add variables of vertical, terrain-following coordinates 
! transformation equation and stretching function.

    SUBROUTINE addNetCDFvertical( ncid, zr_dimid, zw_dimid )
                               
!    input parameters
      integer, intent( in) :: ncid
      integer, intent( in) :: zr_dimid, zw_dimid

      integer :: var_id
     
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


      call check( nf90_def_var(ncid, 's_rho', NF90_DOUBLE, zr_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 's_w', NF90_DOUBLE, zw_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'positive',  'up' ) )
      call check( nf90_put_att(ncid, var_id, 'standard_name', 'ocean_s_coordinate_g1') )
      call check( nf90_put_att(ncid, var_id, 'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc') )

      call check( nf90_def_var(ncid, 'Cs_r', NF90_DOUBLE, zr_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at RHO-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )

      call check( nf90_def_var(ncid, 'Cs_w', NF90_DOUBLE, zw_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'S-coordinate stretching curves at W-points') )
      call check( nf90_put_att(ncid, var_id, 'valid_min',  -1.0d0 ) )
      call check( nf90_put_att(ncid, var_id, 'valid_max',   0.0d0 ) )

    END SUBROUTINE addNetCDFvertical
      
!**** create initial conditions NetCDF file **********************************************

    SUBROUTINE createNetCDFini(   &
!        input parameters
            OUT_FILE              &
          , TIME_ATT              &  
          , Nxr, Nyr, Nzr, Nt     &   
      )
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Nxr, Nyr, Nzr, Nt
      
      integer :: ncid,var_id
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xr_dimid, yr_dimid
      integer :: xu_dimid, yu_dimid
      integer :: xv_dimid, yv_dimid
      integer :: zr_dimid, zw_dimid
      integer :: t_dimid
      integer :: dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xr_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Nxr-1, xu_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Nxr, xv_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Nyr, yr_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Nyr, yu_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Nyr-1, yv_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nzr, zr_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nzr+1, zw_dimid) )
      call check( nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, t_dimid) )
      
    ! Define the netCDF variables.
      call addNetCDFvertical( ncid, zr_dimid, zw_dimid )

      call check( nf90_def_var(ncid, 'ocean_time', NF90_DOUBLE, t_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      dim3Dids = (/ xr_dimid, yr_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'zeta', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim3Dids = (/ xu_dimid, yu_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'ubar', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated u-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim3Dids = (/ xv_dimid, yv_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'vbar', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated v-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xr_dimid, yr_dimid, zr_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'temp', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature') )
      call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'salt', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity') )
      call check( nf90_put_att(ncid, var_id, 'units',     'psu') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xu_dimid, yu_dimid, zr_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'u', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      dim4Dids = (/ xv_dimid, yv_dimid, zr_dimid, t_dimid /)

      call check( nf90_def_var(ncid, 'v', NF90_DOUBLE, dim4Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum component') )
      call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFini
!
!**** Create initial conditions NetCDF file **********************************************

    SUBROUTINE createNetCDFini2( &
!        input parameters
            IN_FILE              &
          , OUT_FILE             &
          , TIME_ATT             &  
          , Nxr, Nyr, Nzr, Nt    &   
          , name_flag            &   
      )
                               
!    input parameters
      character(len=*),  intent( in) :: IN_FILE
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Nxr, Nyr, Nzr, Nt
      integer, intent( in) :: name_flag( N_var )
      character(256) :: varname
      character(2) :: varnum

      integer :: ncid,var_id, ncid2,var_id2
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xr_dimid, yr_dimid
      integer :: xu_dimid, yu_dimid
      integer :: xv_dimid, yv_dimid
      integer :: zr_dimid, zw_dimid
      integer :: t_dimid
      integer :: status
      integer, allocatable :: dimids(:)
      integer :: i,j

!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", trim( OUT_FILE )

      call check( nf90_create( trim( OUT_FILE ), nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xr_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Nxr-1, xu_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Nxr, xv_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Nyr, yr_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Nyr, yu_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Nyr-1, yv_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nzr, zr_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nzr+1, zw_dimid) )
      call check( nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, t_dimid) )
      
    ! Define the netCDF variables.
      call addNetCDFvertical( ncid, zr_dimid, zw_dimid )

      call check( nf90_def_var(ncid, 'ocean_time', NF90_DOUBLE, t_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

!      call check( nf90_redef(ncid) )

      call check( nf90_open( trim(IN_FILE), nf90_nowrite, ncid2) )
      do i=1, N_var
        
        if (name_flag( i ) == 0 ) cycle

        if     ( i == 1 ) then  ! zeta
          allocate(dimids(3))
          dimids = (/ xr_dimid, yr_dimid, t_dimid /)
        elseif ( i == 2 ) then  ! u
          allocate(dimids(4))
          dimids = (/ xu_dimid, yu_dimid, zr_dimid, t_dimid /)
        elseif ( i == 3 ) then  ! v
          allocate(dimids(4))
          dimids = (/ xv_dimid, yv_dimid, zr_dimid, t_dimid /)
        elseif ( i == 4 ) then  ! ubar
          allocate(dimids(3))
          dimids = (/ xu_dimid, yu_dimid, t_dimid /)
        elseif ( i == 5 ) then  ! vbar
          allocate(dimids(3))
          dimids = (/ xv_dimid, yv_dimid, t_dimid /)
        else
          allocate(dimids(4))   ! tracers
          dimids = (/ xr_dimid, yr_dimid, zr_dimid, t_dimid /)
        endif

        if( i == 4 ) then  ! ubar
          write(*,*) 'Add variable: ', trim( VAR_NAME(i) )
          call check( nf90_def_var(ncid, trim( VAR_NAME(i) ), NF90_DOUBLE, dimids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated u-momentum component') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
        else if( i == 5 ) then  ! vbar
          write(*,*) 'Add variable: ', trim( VAR_NAME(i) )
          call check( nf90_def_var(ncid, trim( VAR_NAME(i) ), NF90_DOUBLE, dimids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'vertically integrated v-momentum component') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
        else if( i==8  .or. i==9  .or. i==13 .or. i==14 .or. i==15 .or. i==16 .or.  &
                 i==17 .or. i==21 .or. i==22 .or. i==23 .or. i==24 .or. i==26 .or.  &
                 i==27 .or. i==28 .or. i==29 .or. i==30 .or. i==33 .or. i==34 .or.  &
                 i==35 .or. i==36  ) then  ! mud_, sand_, etc...
           do j=1,99
            write(varnum,'(I2.2)') j
            varname = trim( VAR_NAME(i) )//varnum
            status = nf90_inq_varid(ncid2, trim( varname ), var_id2)
            if (status /= nf90_noerr) exit
            write(*,*) 'Add variable: ', trim( varname )
            call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dimids, var_id) )
            call check( nf90_copy_att(ncid2, var_id2, 'long_name', ncid, var_id) )
            call check( nf90_copy_att(ncid2, var_id2, 'units', ncid, var_id) )
            call check( nf90_copy_att(ncid2, var_id2, 'time', ncid, var_id) )
          enddo
        else
          write(*,*) 'Add variable: ', trim( VAR_NAME(i) )
          call check( nf90_inq_varid(ncid2, trim( VAR_NAME(i) ), var_id2) )
          call check( nf90_def_var(ncid, trim( VAR_NAME(i) ), NF90_DOUBLE, dimids, var_id) )
          call check( nf90_copy_att(ncid2, var_id2, 'long_name', ncid, var_id) )
        !  call check( nf90_copy_att(ncid2, var_id2, 'units', ncid, var_id) )
          status = nf90_copy_att(ncid2, var_id2, 'units', ncid, var_id) 
          call check( nf90_copy_att(ncid2, var_id2, 'time', ncid, var_id) )
        endif
        deallocate(dimids)
      enddo

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFini2
             
!**** create boundary conditions NetCDF file **********************************************

    SUBROUTINE createNetCDFbry(  &
!        input parameters
            OUT_FILE             &
          , TIME_ATT             &  
          , Nxr, Nyr, Nzr, Nt    &   
          , name_flag            &   
          , snwe_flag            &   
      )
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Nxr, Nyr, Nzr, Nt
      integer, intent( in) :: name_flag( 7 )
      integer, intent( in) :: snwe_flag(4)
      
      integer :: ncid,var_id
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xr_dimid, yr_dimid
      integer :: xu_dimid, yu_dimid
      integer :: xv_dimid, yv_dimid
      integer :: zr_dimid, zw_dimid
      integer :: bry_time_dimid
      integer :: dim2Dids(2), dim3Dids(3), dim4Dids(4)
      
!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xr_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Nxr-1, xu_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Nxr, xv_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Nyr, yr_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Nyr, yu_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Nyr-1, yv_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nzr, zr_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nzr+1, zw_dimid) )
      call check( nf90_def_dim(ncid, 'bry_time', NF90_UNLIMITED, bry_time_dimid) )
      
    ! Define the netCDF variables.
      call addNetCDFvertical( ncid, zr_dimid, zw_dimid )

      call check( nf90_def_var(ncid, 'bry_time', NF90_DOUBLE, bry_time_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'open boundary conditions time') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      if( name_flag(1) == 1 ) then
        dim2Dids = (/ yr_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'zeta_west', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'zeta_east', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim2Dids = (/ xr_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'zeta_south', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'zeta_north', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'free-surface northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

      if( name_flag(4) == 1 ) then
        dim2Dids = (/ yu_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'ubar_west', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'ubar_east', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim2Dids = (/ xu_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'ubar_south', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'ubar_north', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D u-momentum northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

      if( name_flag(5) == 1 ) then
        dim2Dids = (/ yv_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'vbar_west', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'vbar_east', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim2Dids = (/ xv_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'vbar_south', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'vbar_north', NF90_DOUBLE, dim2Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', '2D v-momentum northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif         
      endif

      if( name_flag(2) == 1 ) then
        dim3Dids = (/ yu_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'u_west', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'u_east', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim3Dids = (/ xu_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'u_south', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'u_north', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'u-momentum northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

      if( name_flag(3) == 1 ) then
        dim3Dids = (/ yv_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'v_west', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'v_east', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim3Dids = (/ xv_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'v_south', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'v_north', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'v-momentum northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

      if( name_flag(6) == 1 ) then
        dim3Dids = (/ yr_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'temp_west', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'temp_east', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim3Dids = (/ xr_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'temp_south', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'temp_north', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'potential temperature northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'units',     'Celsius') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

      if( name_flag(7) == 1 ) then
        dim3Dids = (/ yr_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(3) == 1 ) then
          call check( nf90_def_var(ncid, 'salt_west', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity western boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(4) == 1 ) then
          call check( nf90_def_var(ncid, 'salt_east', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity eastern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        dim3Dids = (/ xr_dimid, zr_dimid, bry_time_dimid /)
        if( snwe_flag(1) == 1 ) then
          call check( nf90_def_var(ncid, 'salt_south', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity southern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
        if( snwe_flag(2) == 1 ) then
          call check( nf90_def_var(ncid, 'salt_north', NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', 'salinity northern boundary condition') )
          call check( nf90_put_att(ncid, var_id, 'time',      'bry_time') )
        endif
      endif

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFbry

!**** Create boundary conditions NetCDF file **********************************************

    SUBROUTINE createNetCDFbry2( &
!        input parameters
            IN_FILE              &
          , OUT_FILE             &
          , TIME_ATT             &  
          , Nxr, Nyr, Nzr, Nt    &   
          , name_flag            &   
          , snwe_flag            &   
      )
                               
!    input parameters
      character(len=*),  intent( in) :: IN_FILE
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Nxr, Nyr, Nzr, Nt
      integer, intent( in) :: name_flag( N_var )
      integer, intent( in) :: snwe_flag(4)
      character(256) :: varname, varname2, longname
      character(2) :: varnum

      integer :: ncid,var_id, ncid2,var_id2
!      integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
      integer :: xr_dimid, yr_dimid
      integer :: xu_dimid, yu_dimid
      integer :: xv_dimid, yv_dimid
      integer :: zr_dimid, zw_dimid
      integer :: t_dimid
      integer :: status
      integer, allocatable :: dimids(:)
      integer :: i,j
      integer :: ibry

!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", trim( OUT_FILE )

      call check( nf90_create( trim( OUT_FILE ), nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Nxr, xr_dimid) )
      call check( nf90_def_dim(ncid, 'xi_u', Nxr-1, xu_dimid) )
      call check( nf90_def_dim(ncid, 'xi_v', Nxr, xv_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Nyr, yr_dimid) )
      call check( nf90_def_dim(ncid, 'eta_u', Nyr, yu_dimid) )
      call check( nf90_def_dim(ncid, 'eta_v', Nyr-1, yv_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nzr, zr_dimid) )
      call check( nf90_def_dim(ncid, 's_w', Nzr+1, zw_dimid) )
      call check( nf90_def_dim(ncid, 'bry_time', NF90_UNLIMITED, t_dimid) )
      
    ! Define the netCDF variables.
      call addNetCDFvertical( ncid, zr_dimid, zw_dimid )

      call check( nf90_def_var(ncid, 'bry_time', NF90_DOUBLE, t_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

!      call check( nf90_redef(ncid) )

      call check( nf90_open( trim(IN_FILE), nf90_nowrite, ncid2) )
      do i=1, N_var
        
        if ( name_flag( i ) == 0 ) cycle

        do ibry=1,4

          if ( snwe_flag( ibry ) == 0 ) cycle

          if     ( i == 1 ) then  ! zeta
            allocate(dimids(2))
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xr_dimid, t_dimid /)
            else
              dimids = (/ yr_dimid, t_dimid /)
            endif
          elseif ( i == 2 ) then  ! u
            allocate(dimids(3))
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xu_dimid, zr_dimid, t_dimid /)
            else
              dimids = (/ yu_dimid, zr_dimid, t_dimid /)
            endif
          elseif ( i == 3 ) then  ! v
            allocate(dimids(3))
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xv_dimid, zr_dimid, t_dimid /)
            else
              dimids = (/ yv_dimid, zr_dimid, t_dimid /)
            endif
          elseif ( i == 4 ) then  ! ubar
            allocate(dimids(2))
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xu_dimid, t_dimid /)
            else
              dimids = (/ yu_dimid, t_dimid /)
            endif
          elseif ( i == 5 ) then  ! vbar
            allocate(dimids(2))
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xv_dimid, t_dimid /)
            else
              dimids = (/ yv_dimid, t_dimid /)
            endif
          else
            allocate(dimids(3))   ! tracers
            if(ibry==1 .or. ibry==2 ) then
              dimids = (/ xr_dimid, zr_dimid, t_dimid /)
            else
              dimids = (/ yr_dimid, zr_dimid, t_dimid /)
            endif
          endif
          if( i == 4 .or. i ==5 ) then  ! ubar, vbar
            varname = trim( VAR_NAME(i) )//'_'//trim( BRY_NAME(ibry) )
            write(*,*) 'Add variable: ', trim( varname )
            call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dimids, var_id) )
            if(i==4) then
              longname = '2D u-momentum '//BRY_LONGNAME(ibry)
            elseif(i==5) then
              longname = '2D v-momentum '//BRY_LONGNAME(ibry)
            else
              longname = trim(longname)//' '//BRY_LONGNAME(ibry)
            endif
            call check( nf90_put_att(ncid, var_id, 'long_name', longname) )
            call check( nf90_put_att(ncid, var_id, 'units',     'meter second-1') )
            call check( nf90_put_att(ncid, var_id, 'time', 'bry_time') )
          else if( i==8  .or. i==9  .or. i==13 .or. i==14 .or. i==15 .or. i==16 .or.  &
                   i==17 .or. i==21 .or. i==22 .or. i==23 .or. i==24 .or. i==26 .or.  &
                   i==27 .or. i==28 .or. i==29 .or. i==30 .or. i==33 .or. i==34 .or.  &
                   i==35 .or. i==36  ) then  ! mud_, sand_, etc...
            do j=1,99
              write(varnum,'(I2.2)') j
              varname2 = trim( VAR_NAME(i) )//varnum
              varname = trim( VAR_NAME(i) )//trim( BRY_NAME(ibry) )//'_'//varnum
              status = nf90_inq_varid(ncid2, trim( varname2 ), var_id2)
              if (status /= nf90_noerr) exit
              write(*,*) 'Add variable: ', trim( varname )
              call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dimids, var_id) )
              call check( nf90_get_att(ncid2, var_id2, 'long_name', longname) )
              longname = trim(longname)//' '//BRY_LONGNAME(ibry)
              call check( nf90_put_att(ncid, var_id, 'long_name', longname) )
              call check( nf90_copy_att(ncid2, var_id2, 'units', ncid, var_id) )
              call check( nf90_put_att(ncid, var_id, 'time', 'bry_time') )
            enddo
          else
            varname = trim( VAR_NAME(i) )//'_'//trim( BRY_NAME(ibry) )
            write(*,*) 'Add variable: ', trim( varname )
            call check( nf90_inq_varid(ncid2, trim( VAR_NAME(i) ), var_id2) )
            call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dimids, var_id) )
            call check( nf90_get_att(ncid2, var_id2, 'long_name', longname) )
            if(i==2) then
              longname = 'u-momentum '//BRY_LONGNAME(ibry)
            elseif(i==3) then
              longname = 'v-momentum '//BRY_LONGNAME(ibry)
            else
              longname = trim(longname)//' '//BRY_LONGNAME(ibry)
            endif
            call check( nf90_put_att(ncid, var_id, 'long_name', longname) )
            status = nf90_copy_att(ncid2, var_id2, 'units', ncid, var_id) 
            call check( nf90_put_att(ncid, var_id, 'time', 'bry_time') )
          endif
          deallocate(dimids)

        enddo
      enddo

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      call check( nf90_close(ncid2) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFbry2

!**** Create river forcing NetCDF file **********************************************

    SUBROUTINE createNetCDFriver( &
!        input parameters
            OUT_FILE              &
          , GLOBAL_ATT            &  
          , TIME_ATT              &  
          , Nriv, Nzr             &   
          , name_flag             &   
          , NCS, NNS              &   
      )
                               
!    input parameters
      character(len=*), intent( in) :: OUT_FILE
      character(len=*), intent( in) :: GLOBAL_ATT
      character(len=*), intent( in) :: TIME_ATT
      integer, intent( in) :: Nriv, Nzr
      integer, intent( in) :: name_flag( N_var )
      integer, intent( in) :: NCS  ! Number of cohesive (mud) sediment tracers  !!! changed for Ando-kun's simulation
      integer, intent( in) :: NNS  ! Number of non-cohesive (sand) sediment tracers

      character(256) :: varname, lvarname
      character(2) :: varnum

      integer :: ncid,var_id
      integer :: zr_dimid
      integer :: riv_dimid
      integer :: t_dimid
      integer :: dim2Dids(2), dim3Dids(3)
      integer :: i,j
      integer :: Ntype

!---- Create the ROMS initial condition netCDF file --------------------------------

      write(*,*) "CREATE: ", trim( OUT_FILE )

      call check( nf90_create( trim( OUT_FILE ), nf90_clobber, ncid) )

      call check( nf90_put_att(ncid, NF90_GLOBAL, 'rivers', trim( GLOBAL_ATT )) )

      call check( nf90_def_dim(ncid, 'river', Nriv, riv_dimid) )
      call check( nf90_def_dim(ncid, 's_rho', Nzr, zr_dimid) )
      call check( nf90_def_dim(ncid, 'river_time', NF90_UNLIMITED, t_dimid) )

      call check( nf90_def_var(ncid, 'river', NF90_DOUBLE, riv_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river runoff identification number') )

      call check( nf90_def_var(ncid, 'river_time', NF90_DOUBLE, t_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river runoff time') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      call check( nf90_def_var(ncid, 'river_direction', NF90_DOUBLE, riv_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river runoff direction') )
      call check( nf90_put_att(ncid, var_id, 'flag_values', '0, 1') )
      call check( nf90_put_att(ncid, var_id, 'flag_meanings', 'flow across u-face, flow across v-face') )
      call check( nf90_put_att(ncid, var_id, 'LwSrc_True', 'flag not used') )

      call check( nf90_def_var(ncid, 'river_Xposition', NF90_DOUBLE, riv_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river XI-position') )
      call check( nf90_put_att(ncid, var_id, 'LuvSrc_meaning', 'i point index of U or V face source/sink') )
      call check( nf90_put_att(ncid, var_id, 'LwSrc_meaning', 'i point index of RHO center source/sink') )

      call check( nf90_def_var(ncid, 'river_Eposition', NF90_DOUBLE, riv_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river ETA-position') )
      call check( nf90_put_att(ncid, var_id, 'LuvSrc_meaning', 'j point index of U or V face source/sink') )
      call check( nf90_put_att(ncid, var_id, 'LwSrc_meaning', 'j point index of RHO center source/sink') )

      dim2Dids = (/ riv_dimid, zr_dimid /)
      call check( nf90_def_var(ncid, 'river_Vshape', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river runoff mass transport vertical profile') )
      call check( nf90_put_att(ncid, var_id, 'requires', 'must sum to 1 over s_rho') )

      dim2Dids = (/ riv_dimid, t_dimid /)
      call check( nf90_def_var(ncid, 'river_transport', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'river ETA-position') )
      call check( nf90_put_att(ncid, var_id, 'units', 'meter3 second-1' ) )
      call check( nf90_put_att(ncid, var_id, 'LuvSrc_meaning', 'j point index of U or V face source/sink') )
      call check( nf90_put_att(ncid, var_id, 'LwSrc_meaning', 'j point index of RHO center source/sink') )  

      dim3Dids = (/ riv_dimid, zr_dimid, t_dimid /)

      do i=6, N_var

        if (name_flag( i ) == 0 ) cycle
        
        if( i==8  .or. i==9  .or. i==13 .or. i==14 .or. i==15 .or. i==16 .or.  &
            i==17 .or. i==21 .or. i==22 .or. i==23 .or. i==24 .or. i==26 .or.  &
            i==27 .or. i==28 .or. i==29 .or. i==30 .or. i==33 .or. i==34 .or.  &
            i==35 .or. i==36  ) then  ! mud_, sand_, etc...
          
          if( i==8 ) Ntype = NCS
          if( i==9 ) Ntype = NNS
          if( i==13 .or. i==21 .or. i==23 .or. i==26 .or. i==33  ) Ntype = Ndom
          if( i==14 .or. i==22 .or. i==24 .or. i==27 .or. i==34  ) Ntype = Npom
          if( i==15 .or. i==28 .or. i==35  ) Ntype = Nphy
          if( i==16 .or. i==29 .or. i==36  ) Ntype = Nzoo
          if( i==17  ) Ntype = Npim

          do j=1,Ntype
            write(varnum,'(I2.2)') j
            varname = 'river_'//trim( VAR_NAME(i) )//varnum
            lvarname = 'river runoff '//trim( VAR_NAME(i) )//varnum

            write(*,*) 'Add variable: ', trim( varname )
            call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dim3Dids, var_id) )
            call check( nf90_put_att(ncid, var_id, 'long_name', trim( lvarname ) ) )
            call check( nf90_put_att(ncid, var_id, 'units', trim( VAR_UNIT(i) ) ) )
          enddo
        else
          varname = 'river_'//trim( VAR_NAME(i) )
          lvarname = 'river runoff '//trim( VAR_NAME(i) )

          write(*,*) 'Add variable: ', trim( varname  )       
          call check( nf90_def_var(ncid, trim( varname ), NF90_DOUBLE, dim3Dids, var_id) )
          call check( nf90_put_att(ncid, var_id, 'long_name', trim( lvarname ) ) )
          call check( nf90_put_att(ncid, var_id, 'units', trim( VAR_UNIT(i) ) ) )
        endif
      enddo

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

    END SUBROUTINE createNetCDFriver
                   
!**** create HYCOM NetCDF **********************************************

  SUBROUTINE createNetCDF_HYCOM(   &
!    input parameters
        OUT_FILE             &
      , TIME_ATT             &  
      , Im, Jm, Nz, Nt       &   
  )
                           
! input parameters
    character(len=*),  intent( in) :: OUT_FILE
    character(len=*),  intent( in) :: TIME_ATT
    integer, intent( in) :: Im, Jm, Nz, Nt
    
    integer :: ncid2,var_id2
    integer :: lat_dimid, lon_dimid, depth_dimid, time_dimid
    integer :: dim3Dids(3), dim4Dids(4)
  
!---- Create the extracted HYCOM netCDF file --------------------------------

    write(*,*) "CREATE: ", trim( OUT_FILE )
  
    call check( nf90_create(trim( OUT_FILE ), nf90_clobber, ncid2) )
  
    call check( nf90_def_dim(ncid2, 'lat', Jm, lat_dimid) )
    call check( nf90_def_dim(ncid2, 'lon', Im, lon_dimid) )
    call check( nf90_def_dim(ncid2, 'depth',Nz, depth_dimid) )
    call check( nf90_def_dim(ncid2, 'time', NF90_UNLIMITED, time_dimid) )
  
    dim3Dids = (/ lon_dimid, lat_dimid, time_dimid /)
    dim4Dids = (/ lon_dimid, lat_dimid, depth_dimid, time_dimid /)
    
!   Define the netCDF variables.
    call check( nf90_def_var(ncid2, 'time', NF90_DOUBLE, time_dimid, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Valid Time') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     TIME_ATT ) )
  
    call check( nf90_def_var(ncid2, 'depth', NF90_DOUBLE, depth_dimid, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Depth') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'meter' ) )
  
    call check( nf90_def_var(ncid2, 'lon', NF90_DOUBLE, lon_dimid, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Longitude') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'degrees_east') )
  
    call check( nf90_def_var(ncid2, 'lat', NF90_DOUBLE, lat_dimid, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Latitude') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'degrees_north' ) )
  
    call check( nf90_def_var(ncid2, 'surf_el', NF90_DOUBLE, dim3Dids, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Water Surface Elevation') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'meter') )
    call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )
  
    call check( nf90_def_var(ncid2, 'water_temp', NF90_DOUBLE, dim4Dids, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Water Temperature') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'degC') )
    call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )
  
    call check( nf90_def_var(ncid2, 'salinity', NF90_DOUBLE, dim4Dids, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Salinity') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'psu') )
    call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )
  
    call check( nf90_def_var(ncid2, 'water_u', NF90_DOUBLE, dim4Dids, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Eastward Water Velocity') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'm/s') )
    call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )
  
    call check( nf90_def_var(ncid2, 'water_v', NF90_DOUBLE, dim4Dids, var_id2) )
    call check( nf90_put_att(ncid2, var_id2, 'long_name', 'Northward Water Velocity') )
    call check( nf90_put_att(ncid2, var_id2, 'units',     'm/s') )
    call check( nf90_put_att(ncid2, var_id2, 'time',      'time') )
  
    ! End define mode.
    call check( nf90_enddef(ncid2) )
    call check( nf90_close(ncid2) )
    write(*,*) '*** SUCCESS'
  
  END SUBROUTINE createNetCDF_HYCOM
     
!**** writeNetCDF_1d **********************************************
      
    SUBROUTINE writeNetCDF_1d(   &
!        input parameters
            NCNAME                 &
          , OUT_FILE               &
          , Nxr                     &
          , data                   &
          , start1D, count1D       &
      )
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Nxr 
      real(8), intent( in) :: data(Nxr )
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
            NCNAME               &
          , OUT_FILE             &
          , Nxr, Nyr             &
          , data                 &
          , start2D, count2D     &
      )
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Nxr, Nyr
      real(8), intent( in) :: data(Nxr, Nyr)
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
            NCNAME               &
          , OUT_FILE             &
          , Nxr, Nyr, Nt         &
          , data                 &
          , start3D, count3D     &
      )
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Nxr, Nyr, Nt 
      real(8), intent( in) :: data(Nxr, Nyr, Nt )
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
            NCNAME               &
          , OUT_FILE             &
          , Nxr, Nyr, Nzr, Nt    &
          , data                 &
          , start4D, count4D     &
      )
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Nxr, Nyr, Nzr, Nt
      real(8), intent( in) :: data(Nxr, Nyr, Nzr, Nt)
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
!**** readNetCDF_2d **********************************************
      
  SUBROUTINE readNetCDF_2d(    &
!      input parameters
          ncid                   &
        , NCNAME                 &
        , Im, Jm                 &
        , start, count           &
!      output parameters
        , data                   &
    )
                               
!    input parameters
    integer, intent( in) :: ncid
    character(len=*), intent( in) :: NCNAME
    integer, intent( in) :: Im, Jm
    integer, intent( in) :: start(2), count(2)
    real(8), intent(out) :: data(Im, Jm)
    
    integer, parameter :: Num_try = 30
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
#if !defined HYCOM_LOCAL        
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
#endif
    write(*,*) '*** SUCCESS'

  END SUBROUTINE readNetCDF_2d

!**** readNetCDF_3d **********************************************
      
  SUBROUTINE readNetCDF_3d_hycom(    &
!    input parameters
        ncid                   &
      , NCNAME                 &
      , Im, Jm, Nt             &
      , start, count           &
!    output parameters
      , data                   &
    )
                           
!    input parameters
    integer, intent( in) :: ncid
    character(len=*), intent( in) :: NCNAME
    integer, intent( in) :: Im, Jm, Nt
    integer, intent( in) :: start(3), count(3)
    real(8), intent(out) :: data(Im, Jm, Nt)
    
    integer, parameter :: Num_try = 30
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
#if !defined HYCOM_LOCAL        
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
#endif
    write(*,*) '*** SUCCESS'

  END SUBROUTINE readNetCDF_3d_hycom

!
!**** readNetCDF_4d ver 2 **********************************************
      
  SUBROUTINE readNetCDF_4d_2(  &
!    input parameters
        ncid                 &
      , NCNAME               &
      , Nxr, Nyr, Nzr, Nt    &
      , start, count         &
!    output parameters
      , data                 &
    )
                               
!    input parameters
    integer, intent( in) :: ncid
    character(len=*), intent( in) :: NCNAME
    integer, intent( in) :: Nxr, Nyr, Nzr, Nt
    integer, intent( in) :: start(4), count(4)
    real(8), intent(out) :: data(Nxr, Nyr, Nzr, Nt)
    
    integer, parameter :: Num_try = 50
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
    do k=1,Nzr
      start2(3)=k
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data(:,:,k,:), start=start2, count=count2)
        write(*,*) trim(nf90_strerror(status)), k
        if (status == nf90_noerr .and. data(Nxr,Nyr,k,Nt)/=-9999.0d0) exit
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

!**** readNetCDF_4d ver 3 **********************************************
      
  SUBROUTINE readNetCDF_4d_hycom(  &
!    input parameters
        ncid                 &
      , NCNAME               &
      , Nxr, Nyr, Nzr, Nt    &
      , start, count         &
!    output parameters
      , data                 &
    )
                               
!    input parameters
    integer, intent( in) :: ncid
    character(len=*), intent( in) :: NCNAME
    integer, intent( in) :: Nxr, Nyr, Nzr, Nt
    integer, intent( in) :: start(4), count(4)
    real(8), intent(out) :: data(Nxr, Nyr, Nzr, Nt)
    
    integer, parameter :: Num_try = 50
    integer :: start2(4), count2(4)
    integer :: var_id
    integer :: err_flag, status
    integer :: itry
    integer :: i,j,k,l
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
    do k=1,Nzr
      start2(3)=k
      do itry=1,Num_try
        status = nf90_get_var(ncid, var_id, data(:,:,Nzr-k+1,:), start=start2, count=count2)
        write(*,*) trim(nf90_strerror(status)), k
        if (status == nf90_noerr .and. data(Nxr,Nyr,Nzr-k+1,Nt)/=-9999.0d0) exit
        if (itry== Num_try) stop
        write(*,*) '*** READ FAILED: Retry!'
      end do
    end do        
#if !defined HYCOM_LOCAL        
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
#endif

    do i=1,Nxr
      do j=1,Nyr
        do k=Nzr-1,1,-1
          do l=1,Nt
            if ( data(i,j,k,l) <-9.0d0 ) then
              data(i,j,k,l) = data(i,j,k+1,l)
            endif
          end do
        end do
      end do
    end do

    write(*,*) '*** SUCCESS'

  END SUBROUTINE readNetCDF_4d_hycom

!**** readNetCDF_4d ver 3 **********************************************
      
  SUBROUTINE readNetCDF_4d_hycom_fast(  &
!    input parameters
        ncid                 &
      , NCNAME               &
      , Nxr, Nyr, Nzr, Nt    &
      , start, count         &
!    output parameters
      , data                 &
    )
                               
!    input parameters
    integer, intent( in) :: ncid
    character(len=*), intent( in) :: NCNAME
    integer, intent( in) :: Nxr, Nyr, Nzr, Nt
    integer, intent( in) :: start(4), count(4)
    real(8), intent(out) :: data(Nxr, Nyr, Nzr, Nt)

    real(8) :: data_tmp(Nxr, Nyr, Nzr, Nt)
    integer, parameter :: Num_try = 50
    integer :: start2(4), count2(4)
    integer :: var_id
    integer :: err_flag, status
    integer :: itry
    integer :: i,j,k,l
    real(8) :: sf, off
    
    start2 = start
    count2 = count
!    count2(3)=1
    
    data_tmp(:,:,:,:) = -9999.0d0
      
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
      status = nf90_get_var(ncid, var_id, data_tmp(:,:,:,:), start=start2, count=count2)
      write(*,*) trim(nf90_strerror(status))
      if (status == nf90_noerr .and. data_tmp(Nxr,Nyr,Nzr,Nt)/=-9999.0d0) exit
!      if (status == nf90_noerr ) exit
      if (itry== Num_try) stop
      write(*,*) '*** READ FAILED: Retry!'
    end do

#if !defined HYCOM_LOCAL            
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

    do k=1,Nzr
      data(:,:,k,:) = data_tmp(:,:,Nzr-k+1,:)*sf+off
    end do      
!    data = data*sf+off
#else
    do k=1,Nzr
      data(:,:,k,:) = data_tmp(:,:,Nzr-k+1,:)
    end do      
#endif

    do i=1,Nxr
      do j=1,Nyr
        do k=Nzr-1,1,-1
          do l=1,Nt
            if ( data(i,j,k,l) <-9.0d0 ) then
              data(i,j,k,l) = data(i,j,k+1,l)
            endif
          end do
        end do
      end do
    end do

    write(*,*) '*** SUCCESS'

  END SUBROUTINE readNetCDF_4d_hycom_fast

!**** readNetCDF_1d **********************************************
      
    SUBROUTINE readNetCDF_1d(   &
!         input parameters
           ncid, NCNAME, Nxr    &
!         output parameters
           , data               &
       )
                               
!    input parameters
      integer, intent( in) :: ncid
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Nxr
      real(8), intent(out) :: data(Nxr)
      
      integer, parameter :: Num_try = 50
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

!**** NetCDF utility **********************************************
            
  SUBROUTINE try_nf_open(NC_FILE, nf90_open_mode, ncid)
    
    character(len=*),  intent( in) :: NC_FILE
    integer,           intent( in) :: nf90_open_mode
    integer,           intent(out) :: ncid

    integer, parameter :: Num_try = 30
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

    integer, parameter :: Num_try = 50
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
    
!    print *, trim(nf90_strerror(status))
    if (status /= nf90_noerr) then 
      write(*,*) trim(nf90_strerror(status))
      stop "Stopped"
    end if
    
  END SUBROUTINE check    
      
END MODULE mod_roms_netcdf
      
! -------------------------------------------------------------------------
