
!!!=== Copyright (c) 2014-2021 Takashi NAKAMURA  ===== 

#if defined HYCOM_MODEL || defined JCOPE_MODEL
# undef WET_DRY
#endif
#if defined ROMS_MODEL || defined JCOPE_MODEL
# define ARAKAWA_C_GRID
#endif
#if defined HYCOM_LOCAL
# undef GOFS_31
# undef GOFS_30
# undef ANALYSIS_Y
# undef ANALYSIS
# undef REANALYSIS
# undef SKIP_CHECK_TIME
#endif 

PROGRAM iniOCN2ROMS
  use netcdf
  use mod_utility
  use mod_roms_netcdf
  use mod_calendar
  use mod_interpolation
#if defined JCOPE_MODEL
  use mod_jcope
#endif
 
  implicit none
      
! ---------------------------------------------------------------------
  integer :: INIyear, INImonth, INIday, INIhour
  integer :: Ryear, Rmonth, Rday
  integer :: itime
  integer :: mode
  character(256) :: GRID_FILE
#if defined ROMS_MODEL
  character(256) :: parent_grid
  integer :: parent_Imin, parent_Imax
  integer :: parent_Jmin, parent_Jmax
  integer :: refine_factor
  character(256) :: ROMS_HISFILE
  integer :: romsvar(N_var)
#else
  integer :: romsvar(7)
#endif

  character(256) :: INI_prefix
  integer :: N_s_rho
  integer :: Nzr
  integer :: spherical
  integer :: Vtransform, Vstretching
  real(8) :: THETA_S , THETA_B, TCLINE, DCRIT
  integer :: Nzr_dg
  integer :: Vtransform_dg, Vstretching_dg
 
  character(33) :: TIME_ATT  = "seconds since 2000-01-01 00:00:00"
  character(15) :: INI_suffix   = "_20000101.00.nc"
  character(256) :: INI_FILE
  
  real(8), allocatable :: time_all(:)
  integer :: start1D(1), count1D(1)
  integer :: start2D(2), count2D(2)
  integer :: start3D(3), count3D(3)
  integer :: start4D(4), count4D(4)
  
  real(8), allocatable :: h(:,:)        ! depth (meter)
  real(8), allocatable :: rmask(:,:)    ! land mask
  real(8), allocatable :: umask(:,:)    ! land mask
  real(8), allocatable :: vmask(:,:)    ! land mask
  real(8), allocatable :: latr(:, :)
  real(8), allocatable :: lonr(:, :)
  real(8), allocatable :: latu(:, :)
  real(8), allocatable :: lonu(:, :)
  real(8), allocatable :: latv(:, :)
  real(8), allocatable :: lonv(:, :)
  real(8), allocatable :: cosAu(:,:)  ! angle differece between lat lon and ROMS coordinates (radian)
  real(8), allocatable :: sinAu(:,:)  ! angle differece between lat lon and ROMS coordinates (radian)
  real(8), allocatable :: cosAv(:,:)  ! angle differece between lat lon and ROMS coordinates (radian)
  real(8), allocatable :: sinAv(:,:)  ! angle differece between lat lon and ROMS coordinates (radian)
  real(8) :: hc       
  real(8), allocatable :: sc_w(:)       
  real(8), allocatable :: sc_r(:)  
  real(8), allocatable :: Cs_w(:)       
  real(8), allocatable :: Cs_r(:)  
  real(8), allocatable :: z_r(:,:,:)
  real(8), allocatable :: z_w(:,:,:)
  real(8), allocatable :: z_u(:,:,:)
  real(8), allocatable :: z_v(:,:,:)
  real(8), allocatable :: dz_u(:,:,:)
  real(8), allocatable :: dz_v(:,:,:)
  real(8), allocatable :: hu(:,:)
  real(8), allocatable :: hv(:,:)
  
  real(8), allocatable :: zeta(:,:,:)   ! free-surface (meter)
  real(8), allocatable :: ubar(:,:,:)   ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar(:,:,:)   ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u(:,:,:,:)    ! u-momentum component (meter second-1)
  real(8), allocatable :: v(:,:,:,:)    ! v-momentum component (meter second-1)
  real(8), allocatable :: t(:,:,:,:)    ! tracer 
  real(8), allocatable :: ull(:,:,:,:)  ! u-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: vll(:,:,:,:)  ! v-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: uu(:,:,:,:)
  real(8), allocatable :: uv(:,:,:,:)
  real(8), allocatable :: vu(:,:,:,:)
  real(8), allocatable :: vv(:,:,:,:)

  real(8), allocatable :: h_dg(:,:)      ! depth (meter) of donor grid
  real(8), allocatable :: rmask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: umask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: vmask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: pmask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: latr_dg(:,:)
  real(8), allocatable :: lonr_dg(:,:)
  real(8), allocatable :: latu_dg(:,:)
  real(8), allocatable :: lonu_dg(:,:)
  real(8), allocatable :: latv_dg(:,:)
  real(8), allocatable :: lonv_dg(:,:)
  real(8), allocatable :: cosAu_dg(:,:)  ! angle differece between lat lon and donor ROMS coordinates (radian)
  real(8), allocatable :: sinAu_dg(:,:)  ! angle differece between lat lon and donor ROMS coordinates (radian)
  real(8), allocatable :: cosAv_dg(:,:)  ! angle differece between lat lon and donor ROMS coordinates (radian)
  real(8), allocatable :: sinAv_dg(:,:)  ! angle differece between lat lon and donor ROMS coordinates (radian)
  real(8) :: hc_dg       
  real(8), allocatable :: sc_w_dg(:)       
  real(8), allocatable :: sc_r_dg(:)  
  real(8), allocatable :: Cs_w_dg(:)       
  real(8), allocatable :: Cs_r_dg(:)  
  real(8), allocatable :: z_r_dg(:,:,:)
  real(8), allocatable :: z_w_dg(:,:,:)
  real(8), allocatable :: z_u_dg(:,:,:)
  real(8), allocatable :: z_v_dg(:,:,:)

  real(8) :: ocean_time(1)                 ! Ocean time (sec)
  real(8), allocatable :: zeta_dg(:,:,:)   ! free-surface (meter)
  real(8), allocatable :: ubar_dg(:,:,:)   ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar_dg(:,:,:)   ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u_dg(:,:,:,:)    ! u-momentum component (meter second-1)
  real(8), allocatable :: v_dg(:,:,:,:)    ! v-momentum component (meter second-1)
  real(8), allocatable :: t_dg(:,:,:,:)    ! tracer 
  real(8), allocatable :: ull_dg(:,:,:,:)  ! u-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: vll_dg(:,:,:,:)  ! v-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: ullu_dg(:,:,:,:)
  real(8), allocatable :: ullv_dg(:,:,:,:)
  real(8), allocatable :: vllu_dg(:,:,:,:)
  real(8), allocatable :: vllv_dg(:,:,:,:)
#if defined WET_DRY
  real(8), allocatable :: rmask_wet_dg(:,:,:)
  real(8), allocatable :: umask_wet_dg(:,:,:)
  real(8), allocatable :: vmask_wet_dg(:,:,:)
#endif
  integer, allocatable :: ID_cnt2Dr(:,:)
  real(8), allocatable :: w_cnt2Dr (:,:)
  integer, allocatable :: ID_cnt3Dr(:,:)
  real(8), allocatable :: w_cnt3Dr (:,:)
  integer, allocatable :: ID_cnt2Du(:,:)
  real(8), allocatable :: w_cnt2Du (:,:)
  integer, allocatable :: ID_cnt2Dv(:,:)
  real(8), allocatable :: w_cnt2Dv (:,:)
  integer, allocatable :: ID_cnt3Du(:,:)
  real(8), allocatable :: w_cnt3Du (:,:)
  integer, allocatable :: ID_cnt3Dv(:,:)
  real(8), allocatable :: w_cnt3Dv (:,:)

  integer :: i,j,k
  integer :: idt,incdf
  real(8) :: d_lat,d_lon
  real(8) :: cff
  
  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh
  character(11) :: YYYYMMDDpHH
  character(12) :: YYYYMMDDhhmm
  character(8) :: YYYYMMDD
  
  integer :: Nxr, Nyr, Nxu, Nyu, Nxv, Nyv
  integer :: L, M, N
  integer :: Ldg, Mdg, Ndg
  integer :: irg, jrg, krg
  integer :: Idg, Jdg, Kdg
  integer :: Nxr_dg, Nyr_dg, Nxu_dg, Nyu_dg, Nxv_dg, Nyv_dg
  integer :: Nt
  integer :: Irdg_min, Irdg_max, Jrdg_min, Jrdg_max  ! Minimum and maximum ID of donor grid index.
  integer :: Iudg_min, Iudg_max, Judg_min, Judg_max  ! Minimum and maximum ID of donor grid index.
  integer :: Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max  ! Minimum and maximum ID of donor grid index.

  real(8) :: d_jdate
  integer :: jdate_Start, jdate_Ref
  real(8) :: d_jdate_Start, d_jdate_Ref
  integer :: ncid,var_id
  integer :: ncid2,var_id2

  character(256) :: varname
  character(2) :: varnum
  integer :: status, access

#if defined HYCOM_MODEL

  TYPE T_NC
    real(8), pointer :: time_all(:)
    integer :: Nt
    integer :: ItS
    integer :: ItE
  END TYPE T_NC
!  TYPE (T_NC) :: NC(NCnum)
  TYPE (T_NC), allocatable :: NC(:)

  real(8), allocatable :: lat(:), lon(:), depth(:),tau(:)
  real(8), allocatable :: time2(:)
  character(256) :: HY_VAR_NAME
  integer :: iNC
  integer :: jdate_20000101
  real(8) :: d_jdate_20000101
# if defined HYCOM_LOCAL
  integer :: NCnum
  character(256), allocatable :: HYCOM_FILE(:)
# endif
#endif

#if defined JCOPE_MODEL
  character(256) :: JCOPE_info_dir, JCOPE_data_dir
  character(256) :: JCOPE_data_file
  real(8), allocatable :: jcope_data(:,:,:)
  real(8), allocatable :: z(:,:,:), zz(:,:,:), dz(:,:,:)
  real(8), allocatable :: lat(:), lon(:)
  character(20) :: JCOPE_sufix
#endif

#if defined NAOTIDE || defined NAOTIDEJ
  real(8), allocatable :: zeta_tide(:,:)   ! free-surface by tide (meter)

! ----- Configulation for NAOTIDE ---------------------------------------
! -----< Set mode, location, epoch, data interval, and output file >-----
! Mode selection
  integer, parameter :: itmode   = 2
!      Set itmode = 1 to calculate geocentric tide
!                 = 2 to calculate pure ocean tide 
!                        with respect to ocean floor
!                 = 3 to calculate radial loading tide 
  integer, parameter :: lpmode  = 1
!      Set lpmode = 1 to use long-period ocean tide map of Takanezawa
!                 = 2 to use equilibrium tide (Valid for itmode = 1,2)
!
  real(8) :: hsp, hlp, jtime
  Logical :: Ldata
  real(8) :: Smjd
#endif

  namelist/grd/GRID_FILE
#if defined ROMS_MODEL
  namelist/refinement/parent_grid
  namelist/refinement/parent_Imin, parent_Imax
  namelist/refinement/parent_Jmin, parent_Jmax
  namelist/refinement/refine_factor
  namelist/roms2roms/ROMS_HISFILE, romsvar
#else
  namelist/ocn2roms/romsvar
#endif
#if defined JCOPE_MODEL
  namelist/jcope/JCOPE_info_dir, JCOPE_data_dir
#endif
#if defined HYCOM_MODEL
# if defined HYCOM_LOCAL
  namelist/hycom_local1/NCnum
  namelist/hycom_local2/HYCOM_FILE
# endif
#endif
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/ini/INI_prefix
  namelist/ini/itime
  namelist/ini/INIyear, INImonth, INIday, INIhour
  namelist/hcoord/spherical
  namelist/zcoord/N_s_rho
  namelist/zcoord/Vtransform, Vstretching
  namelist/zcoord/THETA_S, THETA_B, TCLINE, DCRIT

  ! Read parameters in namelist file
  
  read (5, nml=grd)
  rewind(5)
#if defined ROMS_MODEL
  read (5, nml=refinement)
  rewind(5)
  read (5, nml=roms2roms)
  rewind(5)
#else
  read (5, nml=ocn2roms)
  rewind(5)
#endif
#if defined JCOPE_MODEL
  read (5, nml=jcope)
  rewind(5)
#endif
  read (5, nml=refdate)
  rewind(5)
  read (5, nml=ini)
  rewind(5)
  read (5, nml=hcoord)
  rewind(5)
  read (5, nml=zcoord)
#if defined HYCOM_MODEL
# if defined HYCOM_LOCAL
  rewind(5)
  read (5, nml=hycom_local1)
  allocate( HYCOM_FILE(NCnum) )
  rewind(5)
  read (5, nml=hycom_local2)
# endif
  allocate( NC(NCnum) )
#endif         
!-Modify time-unit description ---------------------------------
  
  write (YYYY, "(I4.4)") Ryear
  write (MM, "(I2.2)") Rmonth
  write (DD, "(I2.2)") Rday
  TIME_ATT(15:18)=YYYY
  TIME_ATT(20:21)=MM
  TIME_ATT(23:24)=DD
  
!-Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
  Nxu  = Nxr-1
  Nyu = Nyr
  Nxv  = Nxr
  Nyv = Nyr-1
  Nzr = N_s_rho

  L = Nxr-1
  M = Nyr-1
  N = Nzr
  
  allocate( latr(0:L, 0:M) )
  allocate( lonr(0:L, 0:M) )
  allocate( latu(1:L, 0:M) )
  allocate( lonu(1:L, 0:M) )
  allocate( latv(0:L, 1:M) )
  allocate( lonv(0:L, 1:M) )
  allocate( rmask(0:L, 0:M) )
  allocate( umask(1:L, 0:M) )
  allocate( vmask(0:L, 1:M) )
  allocate( cosAu(1:L, 0:M) )
  allocate( sinAu(1:L, 0:M) )
  allocate( cosAv(0:L, 1:M) )
  allocate( sinAv(0:L, 1:M) )
  allocate( h(0:L, 0:M) )
  allocate( sc_w(0:N) )
  allocate( sc_r(1:N) )
  allocate( Cs_w(0:N) )
  allocate( Cs_r(1:N) )
  allocate( z_r(0:L, 0:M, 1:N) )
  allocate( z_w(0:L, 0:M, 0:N) )
  allocate( z_u(1:L, 0:M, 1:N) )
  allocate( z_v(0:L, 1:M, 1:N) )
  allocate( dz_u(1:L, 0:M, 1:N) )
  allocate( dz_v(0:L, 1:M, 1:N) )
  allocate( hu(1:L, 0:M) )
  allocate( hv(0:L, 1:M) )

  allocate( zeta(0:L, 0:M, 1) )
  allocate( ubar(1:L, 0:M, 1) )
  allocate( vbar(0:L, 1:M, 1) )
  allocate( u(1:L, 0:M, 1:N, 1) )
  allocate( v(0:L, 1:M, 1:N, 1) )
  allocate( ull(1:L, 0:M, 1:N, 1) )
  allocate( uu (1:L, 0:M, 1:N, 1) )
  allocate( vu (1:L, 0:M, 1:N, 1) )
  allocate( vll(0:L, 1:M, 1:N, 1) )
  allocate( uv (0:L, 1:M, 1:N, 1) )
  allocate( vv (0:L, 1:M, 1:N, 1) )
  allocate( t(0:L, 0:M, 1:N, 1) )

#if defined NAOTIDE || defined NAOTIDEJ
  allocate( zeta_tide(0:L, 0:M) )  
#endif

  ! Get variables
  start2D = (/ 1,  1 /)
  count2D = (/ Nxr, Nyr /)
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latr(0:L, 0:M)) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr) )
  call check( nf90_inq_varid(ncid, 'lat_u', var_id) ) ! latitude at U-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latu) )
  call check( nf90_inq_varid(ncid, 'lon_u', var_id) ) ! longitude at U-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonu) )
  call check( nf90_inq_varid(ncid, 'lat_v', var_id) ) ! latitude at V-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latv) )
  call check( nf90_inq_varid(ncid, 'lon_v', var_id) ) ! longitude at V-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonv) )
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, h) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, rmask) )
  call check( nf90_inq_varid(ncid, 'mask_u', var_id) ) 
  call check( nf90_get_var(ncid, var_id, umask) )
  call check( nf90_inq_varid(ncid, 'mask_v', var_id) ) 
  call check( nf90_get_var(ncid, var_id, vmask) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
  write(*,*) "CLOSE: ", trim( GRID_FILE )
      
  call set_scoord (  Nzr, Vtransform, Vstretching    &
        , THETA_S, THETA_B, TCLINE, DCRIT            &
!    output parameters
        , hc, sc_w, sc_r, Cs_w, Cs_r )

  call set_depth ( Nxr, Nyr, Nzr, h                  &
        , Vtransform, hc, sc_w, sc_r, Cs_w, Cs_r     &
!    output parameters
        , z_r, z_w)
  z_r = -z_r
  z_w = -z_w
  
  do i=1,L
    do j=0,M
      do k=1,N
        z_u(i,j,k) = (z_r(i-1,j,k)+z_r(i,j,k))*0.5d0
        dz_u(i,j,k) = ((z_w(i-1,j,k-1)+z_w(i,j,k-1))-(z_w(i-1,j,k)+z_w(i,j,k)))*0.5d0
      enddo
    enddo
  enddo
  do i=0,L
    do j=1,M
      do k=1,N
        z_v(i,j,k) = (z_r(i,j-1,k)+z_r(i,j,k))*0.5d0
        dz_v(i,j,k) = ((z_w(i,j-1,k-1)+z_w(i,j,k-1))-(z_w(i,j-1,k)+z_w(i,j,k)))*0.5d0
      enddo
    enddo
  enddo

  do i=1,L
    do j=0,M
      d_lat=latr(i,j)-latr(i-1,j)
      d_lon=lonr(i,j)-lonr(i-1,j)
      d_lon=d_lon*cos(latu(i,j)/180.0d0*PI)
      cosAu(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAu(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)

      hu(i,j) = (h(i-1,j)+h(i,j))*0.5d0
      if(abs( hu(i,j) ) < 1.0d-2) hu(i,j) = 0.01d0
    enddo
  enddo
  do i=0,L
    do j=1,M
      d_lat=latr(i,j)-latr(i,j-1)
      d_lon=lonr(i,j-1)-lonr(i,j)
      d_lon=d_lon*cos(latv(i,j)/180.0d0*PI)
      cosAv(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAv(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)

      hv(i,j) = (h(i,j-1)+h(i,j))*0.5d0
      if(abs( hv(i,j) ) < 1.0d-2 ) hv(i,j) = 0.01d0
    enddo         
  enddo


!======= Set initial time and Ending time ========================
  call jd(INIyear, INImonth, INIday, jdate_Start)
  d_jdate_Start = dble(jdate_Start) + dble(INIhour)/24.0d0
  write(*,*) d_jdate_Start

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref
!======= Read grid coorinate of Ocean model ========================

#if defined ROMS_MODEL
  write(*,*) "************ ROMS MODEL ************"

!-Read ROMS parent grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( parent_grid )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( parent_grid ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr_dg)
  call get_dimension(ncid, 'eta_rho', Nyr_dg)
  Nxu_dg = Nxr_dg-1
  Nyu_dg = Nyr_dg
  Nxv_dg = Nxr_dg
  Nyv_dg = Nyr_dg-1
  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1
       
  allocate( latr_dg(0:Ldg, 0:Mdg) )
  allocate( lonr_dg(0:Ldg, 0:Mdg) )
  allocate( latu_dg(1:Ldg, 0:Mdg) )
  allocate( lonu_dg(1:Ldg, 0:Mdg) )
  allocate( latv_dg(0:Ldg, 1:Mdg) )
  allocate( lonv_dg(0:Ldg, 1:Mdg) )
  allocate( rmask_dg(0:Ldg, 0:Mdg) )
  allocate( umask_dg(1:Ldg, 0:Mdg) )
  allocate( vmask_dg(0:Ldg, 1:Mdg) )
  allocate( cosAu_dg(1:Ldg, 0:Mdg) )
  allocate( sinAu_dg(1:Ldg, 0:Mdg) )
  allocate( cosAv_dg(0:Ldg, 1:Mdg) )
  allocate( sinAv_dg(0:Ldg, 1:Mdg) )
  allocate( h_dg(0:Ldg, 0:Mdg) )
  
  ! Get variables
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latr_dg) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr_dg) )
  call check( nf90_inq_varid(ncid, 'lat_u', var_id) ) ! latitude at U-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latu_dg) )
  call check( nf90_inq_varid(ncid, 'lon_u', var_id) ) ! longitude at U-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonu_dg) )
  call check( nf90_inq_varid(ncid, 'lat_v', var_id) ) ! latitude at V-points (degree_east)
  call check( nf90_get_var(ncid, var_id, latv_dg) )
  call check( nf90_inq_varid(ncid, 'lon_v', var_id) ) ! longitude at V-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonv_dg) )
  call check( nf90_inq_varid(ncid, 'h', var_id) ) 
  call check( nf90_get_var(ncid, var_id, h_dg) )
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, rmask_dg) )
  call check( nf90_inq_varid(ncid, 'mask_u', var_id) ) 
  call check( nf90_get_var(ncid, var_id, umask_dg) )
  call check( nf90_inq_varid(ncid, 'mask_v', var_id) ) 
  call check( nf90_get_var(ncid, var_id, vmask_dg) )

  ! Close NetCDF file
  call check( nf90_close(ncid) )

  write(*,*) "CLOSE: ", trim( parent_grid )

!-Read ROMS ocean_his netCDF file --------------------------------

  ! Open NetCDF file
  write(*,*) "OPEN: ", trim( ROMS_HISFILE )
  call check( nf90_open(trim( ROMS_HISFILE ), nf90_nowrite, ncid) )
  call get_dimension(ncid, 'ocean_time', Nt)
  call get_dimension(ncid, 's_rho', Nzr_dg)

  allocate(time_all(Nt))
  Ndg = Nzr_dg
  allocate( sc_w_dg(0:Ndg ) )
  allocate( sc_r_dg(1:Ndg ) )
  allocate( Cs_w_dg(0:Ndg ) )
  allocate( Cs_r_dg(1:Ndg ) )
  allocate( z_r_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
  allocate( z_w_dg(0:Ldg, 0:Mdg, 0:Ndg ) )
  allocate( z_u_dg(1:Ldg, 0:Mdg, 1:Ndg ) )
  allocate( z_v_dg(0:Ldg, 1:Mdg, 1:Ndg ) )

  call check( nf90_inq_varid(ncid, 'ocean_time', var_id) )
  call check( nf90_get_var(ncid, var_id, time_all) )
  call check( nf90_inq_varid(ncid, 's_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, sc_r_dg ) )
  call check( nf90_inq_varid(ncid, 's_w', var_id) )
  call check( nf90_get_var(ncid, var_id, sc_w_dg ) )
  call check( nf90_inq_varid(ncid, 'Cs_r', var_id) )
  call check( nf90_get_var(ncid, var_id, Cs_r_dg ) )
  call check( nf90_inq_varid(ncid, 'Cs_w', var_id) )
  call check( nf90_get_var(ncid, var_id, Cs_w_dg ) )
  call check( nf90_inq_varid(ncid, 'Vtransform', var_id) )
  call check( nf90_get_var(ncid, var_id, Vtransform_dg ) )
  call check( nf90_inq_varid(ncid, 'hc', var_id) )
  call check( nf90_get_var(ncid, var_id, hc_dg ) )
  
  if (itime <= 0 ) then
  
    itime = 1
    do i=2,Nt
      d_jdate = d_jdate_Ref + time_all(i)/86400.0d0
      if( d_jdate > dble(jdate_Start) ) then
        itime = i-1
        exit
      endif
    enddo
  endif
  write(*,*) 'NetCDF time index: ', itime
  
# if defined WET_DRY
  start3D = (/ 1,      1,      itime /)
  count3D = (/ Nxr_dg, Nyr_dg, 1     /)
  call check( nf90_inq_varid(ncid, 'wetdry_mask_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, rmask_dg, start3D, count3D)  )
  start3D = (/ 1,      1,      itime /)
  count3D = (/ Nxu_dg, Nyu_dg, 1     /)
  call check( nf90_inq_varid(ncid, 'wetdry_mask_u', var_id) )
  call check( nf90_get_var(ncid, var_id, umask_dg, start3D, count3D)  )
  start3D = (/ 1,      1,      itime /)
  count3D = (/ Nxv_dg, Nyv_dg, 1     /)
  call check( nf90_inq_varid(ncid, 'wetdry_mask_v', var_id) )
  call check( nf90_get_var(ncid, var_id, vmask_dg, start3D, count3D)  )
# endif
 
  call set_depth ( Nxr_dg, Nyr_dg, Nzr_dg, h_dg                         &
        , Vtransform_dg, hc_dg, sc_w_dg, sc_r_dg, Cs_w_dg, Cs_r_dg   &
!    output parameters
        , z_r_dg, z_w_dg)

  z_r_dg = -z_r_dg
  z_w_dg = -z_w_dg
  
  do i=1,Ldg
    do j=0,Mdg
      do k=1,Ndg
        z_u_dg(i,j,k)= (z_r_dg(i-1,j,k)+z_r_dg(i,j,k))*0.5d0
      enddo
    enddo
  enddo
  do i=0,Ldg
    do j=1,Mdg
      do k=1,Ndg
        z_v_dg(i,j,k)= (z_r_dg(i,j-1,k)+z_r_dg(i,j,k))*0.5d0
      enddo
    enddo
  enddo

  do i=1,Ldg
    do j=0,Mdg
      d_lat=latr_dg(i,j)-latr_dg(i-1,j)
      d_lon=lonr_dg(i,j)-lonr_dg(i-1,j)
      d_lon=d_lon*cos(latu_dg(i,j)/180.0d0*PI)
      cosAu_dg(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAu_dg(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo
  enddo
  do i=0,Ldg
    do j=1,Mdg
      d_lat=latr_dg(i,j)-latr_dg(i,j-1)
      d_lon=lonr_dg(i,j-1)-lonr_dg(i,j)
      d_lon=d_lon*cos(latv_dg(i,j)/180.0d0*PI)
      cosAv_dg(i,j)= d_lat/sqrt(d_lat*d_lat+d_lon*d_lon)
      sinAv_dg(i,j)= d_lon/sqrt(d_lat*d_lat+d_lon*d_lon)
    enddo         
  enddo

  write(*,*) "*************************************" 

#elif defined HYCOM_MODEL
!---- Read HYCOM netCDF file --------------------------------

  write(*,*) "************ HYCOM MODEL ************"
  call jd(2000, 1, 1, jdate_20000101)
  d_jdate_20000101 = dble(jdate_20000101)
  write(*,*) d_jdate_20000101


# if defined GOFS_31
#  if defined ANALYSIS_Y
  open(50, file=trim(HYCOM_TIME_DIR)//'/time_HYCOM_GOF_31_analysisY.dat')
#  elif defined ANALYSIS
  open(50, file=trim(HYCOM_TIME_DIR)//'/time_HYCOM_GOF_31_analysis.dat')
#  elif defined REANALYSIS
  open(50, file=trim(HYCOM_TIME_DIR)//'/time_HYCOM_GOF_31_reanalysis.dat')
#  endif
# elif defined GOFS_30
#  if defined ANALYSIS
  open(50, file=trim(HYCOM_TIME_DIR)//'/time_HYCOM_GOF_30_analysis.dat')
#  elif defined REANALYSIS
  open(50, file=trim(HYCOM_TIME_DIR)//'/time_HYCOM_GOF_30_reanalysis.dat')
#  endif
# endif
# if defined SKIP_CHECK_TIME
  write(*,*) 'READ: Time'
  do iNC=1, NCnum
    read(50,*) NC(iNC)%Nt
    write(*,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    read(50,*) NC(iNC)%time_all
  end do
# else
  do iNC=1, NCnum
    write(*,*) 'CHECK: Time'
    ! Open NetCDF file
    write(*,*) "OPEN: ", HYCOM_FILE(iNC)
    call try_nf_open(HYCOM_FILE(iNC), nf90_nowrite, ncid)
    call get_dimension(ncid, 'time', NC(iNC)%Nt)
    write(*,*) NC(iNC)%Nt
    write(50,*) NC(iNC)%Nt
    allocate( NC(iNC)%time_all(NC(iNC)%Nt) )
    allocate( time2(NC(iNC)%Nt) )
    call readNetCDF_1d(ncid, 'time', NC(iNC)%Nt, time2)
    call check( nf90_close(ncid) )
    write(*,*) "CLOSE: ", HYCOM_FILE(iNC)
    NC(iNC)%time_all = time2
    write(50,*) NC(iNC)%time_all
    deallocate(time2)
  end do
#endif
  close(50)   

  write(*,*) NC(:)%Nt
  
  do iNC=1, NCnum-1
    do i=2,NC(iNC)%Nt
      if( NC(iNC+1)%time_all(1) <= NC(iNC)%time_all(i) ) then
        NC(iNC)%Nt = i-1
        exit
      end if
    end do
  end do
  
  write(*,*) NC(:)%Nt
  
  write(*,*) "******************************************************************"
    
  
  NC(:)%ItE = -1
  NC(:)%ItS = -1
  
!  do iNC=NCnum,1,-1
!    do i=NC(iNC)%Nt-1,1,-1
!      d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24.0d0
!      if(d_jdate < dble(jdate_End)) then
!        write(*,*) '*** FOUND: Ending point @ HYCOM_FILE',iNC
!        NC(iNC)%ItE=i+1
!        exit
!      endif
!    end do
!  end do
!  write(*,*) NC(:)%ItE 
  
  LOOP1 : do iNC=1,NCnum
    do i=2,NC(iNC)%Nt
      d_jdate = d_jdate_20000101 + NC(iNC)%time_all(i)/24.0d0
      if(d_jdate>d_jdate_Start) then
        write(*,*) '*** FOUND: Starting point @ HYCOM_FILE',iNC
        NC(iNC)%ItS=i-1
        itime = NC(iNC)%ItS
        exit LOOP1
      endif
    end do
  end do LOOP1

  write(*,*) NC(:)%ItS 
  write(*,*) 'HYCOM NetCDF #: ', iNC
  write(*,*) 'HYCOM NetCDF time index: ', itime
    
!---- Read HYCOM netCDF grid cooredinate --------------------------------
  ! Open NetCDF file
  write(*,*) "OPEN: ", trim( HYCOM_FILE(iNC) )
  call check( nf90_open(trim( HYCOM_FILE(iNC) ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'lat', Nyr_dg)
  call get_dimension(ncid, 'lon', Nxr_dg)
  call get_dimension(ncid, 'depth', Nzr_dg)
  write(*,*) Nxr_dg, Nyr_dg, Nzr_dg
  ! Allocate variable
  allocate(lat(Nyr_dg))
  allocate(lon(Nxr_dg))
  allocate(depth(Nzr_dg))

  Nxu_dg = Nxr_dg
  Nyu_dg = Nyr_dg
  Nxv_dg = Nxr_dg
  Nyv_dg = Nyr_dg
  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1

  allocate( latr_dg(0:Ldg, 0:Mdg) )
  allocate( lonr_dg(0:Ldg, 0:Mdg) )
  allocate( latu_dg(0:Ldg, 0:Mdg) )
  allocate( lonu_dg(0:Ldg, 0:Mdg) )
  allocate( latv_dg(0:Ldg, 0:Mdg) )
  allocate( lonv_dg(0:Ldg, 0:Mdg) )
  allocate( rmask_dg(0:Ldg, 0:Mdg) )
  allocate( umask_dg(0:Ldg, 0:Mdg) )
  allocate( vmask_dg(0:Ldg, 0:Mdg) )
  allocate( cosAu_dg(0:Ldg, 0:Mdg) )
  allocate( sinAu_dg(0:Ldg, 0:Mdg) )
  allocate( cosAv_dg(0:Ldg, 0:Mdg) )
  allocate( sinAv_dg(0:Ldg, 0:Mdg) )
  allocate( h_dg(0:Ldg, 0:Mdg) )
  Ndg = Nzr_dg
  allocate( z_r_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
!  allocate( z_w_dg(0:Ldg, 0:Mdg, 0:Ndg ) )
  allocate( z_u_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
  allocate( z_v_dg(0:Ldg, 0:Mdg, 1:Ndg ) )

  
!  call check( nf90_inq_varid(ncid, 'lat', var_id) )
!  call check( nf90_get_var(ncid, var_id, lat) )
!  call check( nf90_inq_varid(ncid, 'lon', var_id) )
!  call check( nf90_get_var(ncid, var_id, lon) )
!  call check( nf90_inq_varid(ncid, 'depth', var_id) )
!  call check( nf90_get_var(ncid, var_id, depth) )
  call readNetCDF_1d(ncid, 'lat', Nyr_dg, lat)
  call readNetCDF_1d(ncid, 'lon', Nxr_dg, lon)
  call readNetCDF_1d(ncid, 'depth', Nzr_dg, depth)
  ! Close NetCDF file
!  call check( nf90_close(ncid) )
!  write(*,*) "CLOSE: ", trim( HYCOM_FILE(iNC) )

  do i=0,Ldg
    do j=0,Mdg
      latr_dg(i,j) = lat(j+1)
      lonr_dg(i,j) = lon(i+1)
    enddo
  enddo
  latu_dg(:,:) = latr_dg(:,:)
  lonu_dg(:,:) = lonr_dg(:,:)
  latv_dg(:,:) = latr_dg(:,:)
  lonv_dg(:,:) = lonr_dg(:,:)
  do i=1,Ldg
    do j=0,Mdg
      do k=1,Ndg
        z_r_dg(i,j,k)= depth(Ndg-k+1)
!        z_r_dg(i,j,k)= depth(k)
      enddo
    enddo
  enddo
  z_u_dg(:,:,:) = z_r_dg(:,:,:)
  z_v_dg(:,:,:) = z_r_dg(:,:,:)

  rmask_dg(:,:) = 1.0d0

  cosAu_dg(:,:) = 1.0d0
  sinAu_dg(:,:) = 0.0d0

  cosAv_dg(:,:) = 1.0d0
  sinAv_dg(:,:) = 0.0d0

  write(*,*) "*************************************" 
 
#elif defined JCOPE_MODEL
  write(*,*) "************ JCOPE MODEL ************"

  write(*,*) "Read domain description file: "

  call read_jcope_info( JCOPE_info_dir, Nxr_dg, Nyr_dg, Nzr_dg )

  allocate( lon( Nxr_dg ) )
  allocate( lat( Nyr_dg ) )
  allocate( z( Nxr_dg, Nyr_dg, Nzr_dg ) )
  allocate( zz( Nxr_dg, Nyr_dg, Nzr_dg ) )
  allocate( dz( Nxr_dg, Nyr_dg, Nzr_dg ) )

  Nxu_dg = Nxr_dg-1
  Nyu_dg = Nyr_dg
  Nxv_dg = Nxr_dg
  Nyv_dg = Nyr_dg-1
  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1
  Ndg = Nzr_dg-1
       
  allocate( latr_dg(0:Ldg, 0:Mdg) )
  allocate( lonr_dg(0:Ldg, 0:Mdg) )
  allocate( latu_dg(1:Ldg, 0:Mdg) )
  allocate( lonu_dg(1:Ldg, 0:Mdg) )
  allocate( latv_dg(0:Ldg, 1:Mdg) )
  allocate( lonv_dg(0:Ldg, 1:Mdg) )
  allocate( rmask_dg(0:Ldg, 0:Mdg) )
  allocate( umask_dg(1:Ldg, 0:Mdg) )
  allocate( vmask_dg(0:Ldg, 1:Mdg) )
  allocate( pmask_dg(1:Ldg, 1:Mdg) )
  allocate( cosAu_dg(1:Ldg, 0:Mdg) )
  allocate( sinAu_dg(1:Ldg, 0:Mdg) )
  allocate( cosAv_dg(0:Ldg, 1:Mdg) )
  allocate( sinAv_dg(0:Ldg, 1:Mdg) )
  allocate( h_dg(0:Ldg, 0:Mdg) )
  allocate( z_r_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
!  allocate( z_w_dg(0:Ldg, 0:Mdg, 0:Ndg ) )
  allocate( z_u_dg(1:Ldg, 0:Mdg, 1:Ndg ) )
  allocate( z_v_dg(0:Ldg, 1:Mdg, 1:Ndg ) )

  call read_jcope_latlon( JCOPE_info_dir, Nxr_dg, Nyr_dg, 1, lon, lat )

  do i=0,Ldg
    do j=0,Mdg
      latr_dg(i,j) = lat(j+1)
      lonr_dg(i,j) = lon(i+1)
    enddo
  enddo
  write(*, *) " lonr_dg range : ", lonr_dg(0,0), lonr_dg(Ldg,Mdg)
  write(*, *) " latr_dg range : ", latr_dg(0,0), latr_dg(Ldg,Mdg)
  
  call read_jcope_latlon( JCOPE_info_dir, Nxr_dg, Nyr_dg, 2, lon, lat )

  do i=1,Ldg
    do j=0,Mdg
      latu_dg(i,j) = lat(j+1)
      lonu_dg(i,j) = lon(i+1)
    enddo
  enddo
  write(*, *) " lonu_dg range : ", lonu_dg(1,0), lonu_dg(Ldg,Mdg)
  write(*, *) " latu_dg range : ", latu_dg(1,0), latu_dg(Ldg,Mdg)

  call read_jcope_latlon( JCOPE_info_dir, Nxr_dg, Nyr_dg, 3, lon, lat )

  do i=0,Ldg
    do j=1,Mdg
      latv_dg(i,j) = lat(j+1)
      lonv_dg(i,j) = lon(i+1)
    enddo
  enddo
  write(*, *) " lonv_dg range : ", lonv_dg(0,1), lonv_dg(Ldg,Mdg)
  write(*, *) " latv_dg range : ", latv_dg(0,1), latv_dg(Ldg,Mdg)

  call read_jcope_depth( JCOPE_info_dir, Nxr_dg, Nyr_dg, Nzr_dg, z, zz, dz )

  do k=1,Ndg
    z_r_dg(0:Ldg, 0:Mdg, k ) = -zz(1:Nxr_dg, 1:Nyr_dg, Nzr_dg-k )
  enddo
  h_dg(0:Ldg, 0:Mdg) = -z(1:Nxr_dg, 1:Nyr_dg, Nzr_dg )

  do i=1,Ldg
    do j=0,Mdg
      do k=1,Ndg
        z_u_dg(i,j,k)= (z_r_dg(i-1,j,k)+z_r_dg(i,j,k))*0.5d0
      enddo
    enddo
  enddo
  do i=0,Ldg
    do j=1,Mdg
      do k=1,Ndg
        z_v_dg(i,j,k)= (z_r_dg(i,j-1,k)+z_r_dg(i,j,k))*0.5d0
      enddo
    enddo
  enddo

! Land masking
  rmask_dg(:,:) = 1.0d0
  do i=0,Ldg
    do j=0,Mdg
      if( h_dg(i,j) < 1.001d0 ) then  !! Land: Depth < 1.0 m grids
        rmask_dg(i,j) = 0.0d0
      endif
    enddo
  enddo

  call uvp_masks( Nxr_dg, Nyr_dg, rmask_dg, umask_dg, vmask_dg, pmask_dg )


  cosAu_dg(:,:) = 1.0d0
  sinAu_dg(:,:) = 0.0d0

  cosAv_dg(:,:) = 1.0d0
  sinAv_dg(:,:) = 0.0d0  

#else
  write(*,*) "************ Please choose a OCEAN MODEL ************"
  Stop
#endif
  
! ----------------------------------------------------
  write(*,*) "*** Calculate weight factor ***"

  allocate( ID_cnt2Dr(4, Nxr*Nyr) )
  allocate( w_cnt2Dr (4, Nxr*Nyr) )
  allocate( ID_cnt2Du(4, Nxu*Nyu) )
  allocate( w_cnt2Du (4, Nxu*Nyu) )
  allocate( ID_cnt2Dv(4, Nxv*Nyv) )
  allocate( w_cnt2Dv (4, Nxv*Nyv) )

  allocate( ID_cnt3Dr(6, Nxr*Nyr*Nzr) )
  allocate( w_cnt3Dr (8, Nxr*Nyr*Nzr) )
  allocate( ID_cnt3Du(6, Nxu*Nyu*Nzr) )
  allocate( w_cnt3Du (8, Nxu*Nyu*Nzr) )
  allocate( ID_cnt3Dv(6, Nxv*Nyv*Nzr) )
  allocate( w_cnt3Dv (8, Nxv*Nyv*Nzr) )

  write(*,*) "Seek rho point donor IJ range"
  call seek_IJrange(                                   &
          0, Ldg, 0, Mdg, lonr_dg, latr_dg             & 
        , 0, L,   0, M,   lonr,    latr                &
        , Irdg_min, Irdg_max, Jrdg_min, Jrdg_max)

  Irdg_min = max(Irdg_min-2, 0 )
  Irdg_max = min(Irdg_max+2,Ldg)
  Jrdg_min = max(Jrdg_min-2, 0 )
  Jrdg_max = min(Jrdg_max+2,Mdg)
  Iudg_min = Irdg_min
  Iudg_max = Irdg_max
  Judg_min = Jrdg_min
  Judg_max = Jrdg_max
  Ivdg_min = Irdg_min
  Ivdg_max = Irdg_max
  Jvdg_min = Jrdg_min
  Jvdg_max = Jrdg_max
#if defined ARAKAWA_C_GRID
  Iudg_min = Irdg_min+1
  Jvdg_min = Jrdg_min+1
#endif
  Nxr_dg =Irdg_max-Irdg_min+1
  Nyr_dg =Jrdg_max-Jrdg_min+1
  Nxu_dg =Iudg_max-Iudg_min+1
  Nyu_dg =Judg_max-Judg_min+1
  Nxv_dg =Ivdg_max-Ivdg_min+1
  Nyv_dg =Jvdg_max-Jvdg_min+1

  write(*,*) Irdg_min, Irdg_max, Jrdg_min, Jrdg_max

  allocate( zeta_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, 1) )
  allocate( ubar_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1) )
  allocate( vbar_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1) )
  allocate( u_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
  allocate( v_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
  allocate( ull_dg (Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
  allocate( ullu_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
  allocate( vllu_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
  allocate( vll_dg (Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
  allocate( vllv_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
  allocate( ullv_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
  allocate( t_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, 1:Nzr_dg, 1) )
#if defined WET_DRY
  allocate( rmask_wet_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, 1) )
  allocate( umask_wet_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1) )
  allocate( vmask_wet_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1) )
#endif

#if defined HYCOM_MODEL
  write(*,*) 'Read surf_el and create land mask'
  start3D = (/ Irdg_min+1, Jrdg_min+1, itime /)
  count3D = (/ Nxr_dg,     Nyr_dg,     1     /)
  call readNetCDF_3d_hycom( ncid, 'surf_el'      &
                    , Nxr_dg, Nyr_dg, 1, start3D, count3D, zeta_dg )
  do i=Irdg_min,Irdg_max
    do j=Jrdg_min,Jrdg_max
      if (zeta_dg(i,j,1)<-10.0d0) then
        rmask_dg(i,j) = 0.0d0
      endif
    enddo
  enddo
  umask_dg(:,:) = rmask_dg(:,:)
  vmask_dg(:,:) = rmask_dg(:,:)
#endif

  write(*,*) "Calculate 2D rho point weight factor"
  call weight2D_grid3_2(                                &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , 0, L,   0, M,   lonr,    latr                 &
        , ID_cnt2Dr, w_cnt2Dr )

  write(*,*) "Calculate 2D u point weight factor"
  call weight2D_grid3_2(                                &
          Iudg_min, Iudg_max, Judg_min, Judg_max        &
        , lonu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
        , latu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
        , umask_dg(Iudg_min:Iudg_max,Judg_min:Judg_max) & 
        , 1, L,   0, M,   lonu,    latu                 &
        , ID_cnt2Du, w_cnt2Du )

  write(*,*) "Calculate 2D v point weight factor"
  call weight2D_grid3_2(                                &
          Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max        &
        , lonv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
        , latv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
        , vmask_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max) & 
        , 0, L,   1, M,   lonv,    latv                 &
        , ID_cnt2Dv, w_cnt2Dv )

  write(*,*) "Calculate 3D rho point weight factor"
  call weight3D_grid3_2(                                    &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, 1, Ndg    &
        , z_r_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max,1:Ndg) &
        , 0, L,   0, M,   1, N,   z_r                       &
        , ID_cnt2Dr, w_cnt2Dr                               &
        , ID_cnt3Dr, w_cnt3Dr )

  write(*,*) "Calculate 3D u point weight factor"
  call weight3D_grid3_2(                                    &
          Iudg_min, Iudg_max, Judg_min, Judg_max, 1, Ndg    &
        , z_u_dg(Iudg_min:Iudg_max,Judg_min:Judg_max,1:Ndg) &
        , 1, L,   0, M,   1, N,   z_u                       &
        , ID_cnt2Du, w_cnt2Du                               &
        , ID_cnt3Du, w_cnt3Du )

  write(*,*) "Calculate 3D v point weight factor"
  call weight3D_grid3_2(                                    &
          Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, 1, Ndg    &
        , z_v_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max,1:Ndg) &
        , 0, L,   1, M,   1, N,   z_v                       &
        , ID_cnt2Dv, w_cnt2Dv                               &
        , ID_cnt3Dv, w_cnt3Dv )

  write(*,*) "*********************************************"

#if defined ROMS_MODEL
!-Read ROMS ocean_his netCDF file --------------------------------

  ocean_time(1) = time_all(itime)
!-Initial netcdf file name -------------------------------------------------
  call oceantime2cdate(ocean_time(1), Ryear, Rmonth, Rday, YYYYMMDDpHH)
  
  INI_suffix(2:12)=YYYYMMDDpHH
  INI_FILE = trim( INI_prefix )//INI_suffix

  !---- Create the ROMS initial conditions netCDF file --------------------------------
        
  call createNetCDFini2(  trim( ROMS_HISFILE ),  trim( INI_FILE )   &
        , TIME_ATT , Nxr, Nyr, Nzr, 1 ,romsvar  )
#else
# if defined HYCOM_MODEL
  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  ocean_time(1) =(  NC(iNC)%time_all( NC(iNC)%ItS )           &
                  + ( d_jdate_Ref - d_jdate_20000101 )*24.0d0 )*3600.0d0

# elif defined JCOPE_MODEL

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  call jd(INIyear, INImonth, INIday, jdate_Start)
  d_jdate_Start = dble(jdate_Start) + dble(INIhour)/24.0d0    ! day

  ocean_time(1) =( d_jdate_Start - d_jdate_Ref )*86400.0d0  ! sec

!---- Preparation for bin file name string --------------------------------

  write (YYYY, "(I4.4)") INIyear
  write (MM, "(I2.2)") INImonth
  write (DD, "(I2.2)") INIday
  write (hh, "(I2.2)") INIhour
#  if defined JCOPE_T
  JCOPE_sufix = YYYY//MM//DD//hh//"00.bin"
#  else
  JCOPE_sufix = YYYY//MM//DD
#  endif
! Check data file existance
  JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(1) )//trim( JCOPE_sufix )
  status = access( trim( JCOPE_data_file ), ' ' )
  if ( status /= 0 ) then
    write(*,*) status, 'File not found: ', trim( JCOPE_data_file )
    write(*,*) '*** Please set different initial year/month/day/hour'
    stop
  endif

# endif

  call oceantime2cdate(ocean_time(1), Ryear, Rmonth, Rday, YYYYMMDDpHH)

  INI_suffix(2:12)=YYYYMMDDpHH
  INI_FILE = trim( INI_prefix )//INI_suffix

!---- Create the ROMS initial conditions netCDF file --------------------------------
          
  call createNetCDFini(  trim( INI_FILE ), TIME_ATT, Nxr, Nyr, Nzr, 1 ) 



#endif

!-Write ROMS initial conditions netCDF file --------------------------------
 
  call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
  call check( nf90_inq_varid(ncid2, 'spherical', var_id) )
  call check( nf90_put_var(ncid2, var_id, spherical) )
  call check( nf90_inq_varid(ncid2, 'Vtransform', var_id) )
  call check( nf90_put_var(ncid2, var_id, Vtransform ) )
  call check( nf90_inq_varid(ncid2, 'Vstretching', var_id) )
  call check( nf90_put_var(ncid2, var_id, Vstretching ) )
  call check( nf90_inq_varid(ncid2, 'theta_s', var_id) )
  call check( nf90_put_var(ncid2, var_id, THETA_S ) )
  call check( nf90_inq_varid(ncid2, 'theta_b', var_id) )
  call check( nf90_put_var(ncid2, var_id, THETA_B) )
  call check( nf90_inq_varid(ncid2, 'Tcline', var_id) )
  call check( nf90_put_var(ncid2, var_id, TCLINE) )
  call check( nf90_inq_varid(ncid2, 'hc', var_id) )
  call check( nf90_put_var(ncid2, var_id, hc ) )
  start1D = (/ 1 /)
  count1D = (/ Nzr /)
  call check( nf90_inq_varid(ncid2, 's_rho', var_id) )
  call check( nf90_put_var(ncid2, var_id, sc_r, start1D, count1D) )
  call check( nf90_inq_varid(ncid2, 'Cs_r', var_id) )
  call check( nf90_put_var(ncid2, var_id, Cs_r, start1D, count1D) )
  start1D = (/ 1 /)
  count1D = (/ Nzr+1 /)
  call check( nf90_inq_varid(ncid2, 's_w', var_id) )
  call check( nf90_put_var(ncid2, var_id, sc_w, start1D,  count1D) )
  call check( nf90_inq_varid(ncid2, 'Cs_w', var_id) )
  call check( nf90_put_var(ncid2, var_id, Cs_w, start1D, count1D) )

!-Write ROMS initial conditions netCDF file --------------------------------

  start1D = (/ 1 /)
  count1D = (/ 1 /)
  write(*,*)  'Write: ocean_time', ncid2, var_id2, ocean_time
  call check( nf90_inq_varid(ncid2, 'ocean_time', var_id2) )
  call check( nf90_put_var(ncid2, var_id2, ocean_time, start = start1D, count = count1D) )
  call check( nf90_close(ncid2) )

  if( romsvar(1)==1 ) then
  !- zeta --------------------------------
    write(*,*) 'Read: zeta'
#if defined ROMS_MODEL
    start3D = (/ Irdg_min+1, Jrdg_min+1, itime /)
    count3D = (/ Nxr_dg,     Nyr_dg,     1     /)
    call check( nf90_inq_varid(ncid, 'zeta', var_id) )
    call check( nf90_get_var(ncid, var_id, zeta_dg, start3D, count3D)  )
#elif defined HYCOM_MODEL
!    start3D = (/ Irdg_min+1, Jrdg_min+1, itime /)
!    count3D = (/ Nxr_dg,     Nyr_dg,     1     /)
!    call readNetCDF_3d_hycom( ncid, 'surf_el'      &
!                      , Nxr_dg, Nyr_dg, 1, start3D, count3D, zeta_dg )
#elif defined JCOPE_MODEL
    JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(1) )//trim( JCOPE_sufix )
    write(*,*) 'Read: ', trim( JCOPE_data_file )
    call read_jcope_data2D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
          , Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, zeta_dg )
#endif

    write(*,*) 'Linear Interporation: zeta'
    call interp2D_grid3_2(                                              &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, zeta_dg             &
          , 0, L, 0, M                                                  &
          , Id_cnt2Dr, w_cnt2Dr                                         &
          , zeta  ) 
          
#if defined NAOTIDE || defined NAOTIDEJ
!---- NAOTIDE calculation --------------------------------
    call mjdymd( Smjd, Ryear, Rmonth, Rday , 0    &
               , 0,    0,     1                      )  ! from naotidej.f
        
    jtime = dble(Smjd) + ocean_time(1)/86400.0d0
!    write(*,*) jtime, Smjd

    do i=0,L
      do j=0,M
# if defined NAOTIDE
        call naotide ( lonr(i,j), latr(i,j), jtime  , itmode, lpmode       &
                     , zeta_tide(i,j), hsp   , hlp   , Ldata          )
# elif defined NAOTIDEJ
        call naotidej( lonr(i,j), latr(i,j), jtime  , itmode, lpmode       &
                     , zeta_tide(i,j), hsp   , hlp   , Ldata          )
# endif          
        zeta_tide(i,j)=zeta_tide(i,j)*0.01d0

      enddo
    enddo

    call FillExtraPoints2D(Nxr, Nyr, zeta_tide, -5.0d0, 5.0d0)

    zeta(:,:,1)=zeta(:,:,1)+zeta_tide(:,:)
#endif

    start3D = (/ 1,  1,  1 /)
    count3D = (/ Nxr, Nyr, 1 /)
    write(*,*)  'Write: zeta'
    call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
    call check( nf90_inq_varid(ncid2, 'zeta', var_id2) )
    call check( nf90_put_var(ncid2, var_id2, zeta, start = start3D, count = count3D) )
    call check( nf90_close(ncid2) )
  endif


  

  if( romsvar(2)==1 .or. romsvar(3)==1 .or.         &
      romsvar(4)==1 .or. romsvar(5)==1      ) then
  !- u --------------------------------
  
    write(*,*) 'Read: u'
#if defined ROMS_MODEL
    start4D = (/ Iudg_min, Judg_min+1, 1,      itime /)
    count4D = (/ Nxu_dg,     Nyu_dg,   Nzr_dg, 1     /)
    call check( nf90_inq_varid(ncid, 'u', var_id) )
    call check( nf90_get_var(ncid, var_id, u_dg, start4D, count4D)  )
#elif defined HYCOM_MODEL
    start4D = (/ Iudg_min+1, Judg_min+1,   1,      itime /)
    count4D = (/ Nxu_dg,     Nyu_dg,   Nzr_dg, 1     /)
    call readNetCDF_4d_hycom( ncid, 'water_u'     &
          , Nxu_dg, Nyu_dg, Nzr_dg, 1, start4D, count4D, u_dg )
#elif defined JCOPE_MODEL
    JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(2) )//trim( JCOPE_sufix )
    write(*,*) 'Read: ', trim( JCOPE_data_file )
    call read_jcope_data3D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
          , Iudg_min, Iudg_max, Judg_min, Judg_max, u_dg(:,:,:,1) )
#endif

#if defined WET_DRY
    start3D = (/ Iudg_min, Judg_min+1, itime /)
    count3D = (/ Nxu_dg,     Nyu_dg,     1     /)
    call check( nf90_inq_varid(ncid, 'wetdry_mask_u', var_id) )
    call check( nf90_get_var(ncid, var_id, umask_wet_dg, start3D, count3D)  )
    
    do i=Iudg_min,Iudg_max
      do j=Judg_min,Judg_max
        u_dg(i,j,:,1)= u_dg(i,j,:,1)*umask_wet_dg(i,j,1)
      enddo
    enddo
#endif

  !- v --------------------------------
  
#if defined ROMS_MODEL
    write(*,*) 'Read: v'
    start4D = (/ Ivdg_min+1, Jvdg_min, 1,      itime /)
    count4D = (/ Nxv_dg,     Nyv_dg,   Nzr_dg, 1     /)
    call check( nf90_inq_varid(ncid, 'v', var_id) )
    call check( nf90_get_var(ncid, var_id, v_dg, start4D, count4D)  )
#elif defined HYCOM_MODEL
    write(*,*) 'Read: water_v'
    start4D = (/ Ivdg_min+1, Jvdg_min+1, 1,      itime /)
    count4D = (/ Nxv_dg,     Nyv_dg,   Nzr_dg, 1     /)
    call readNetCDF_4d_hycom( ncid, 'water_v'     &
          , Nxv_dg, Nyv_dg, Nzr_dg, 1, start4D, count4D, v_dg )
#elif defined JCOPE_MODEL
    JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(3) )//trim( JCOPE_sufix )
    write(*,*) 'Read: ', trim( JCOPE_data_file )
    call read_jcope_data3D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
          , Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, v_dg(:,:,:,1) )
#endif

#if defined WET_DRY
    start3D = (/ Ivdg_min+1, Jvdg_min, itime /)
    count3D = (/ Nxv_dg,     Nyv_dg,     1     /)
    call check( nf90_inq_varid(ncid, 'wetdry_mask_v', var_id) )
    call check( nf90_get_var(ncid, var_id, vmask_wet_dg, start3D, count3D)  )
    
    do i=Ivdg_min,Ivdg_max
      do j=Jvdg_min,Jvdg_max
        v_dg(i,j,:,1)= v_dg(i,j,:,1)*vmask_wet_dg(i,j,1)
      enddo
    enddo
#endif

  !- convert ROMS donor coordinate to lat lon coordinate --------------------------------
  
    do i=Iudg_min,Iudg_max
      do j=Judg_min,Judg_max
        ullu_dg(i,j,:,1) = u_dg(i,j,:,1)*cosAu_dg(i,j)
        vllu_dg(i,j,:,1) = u_dg(i,j,:,1)*sinAu_dg(i,j)
      enddo
    enddo
    do i=Ivdg_min,Ivdg_max
      do j=Jvdg_min,Jvdg_max
        ullv_dg(i,j,:,1) = -v_dg(i,j,:,1)*sinAv_dg(i,j)
        vllv_dg(i,j,:,1) =  v_dg(i,j,:,1)*cosAv_dg(i,j)
      enddo
    enddo
#if defined ARAKAWA_C_GRID
    do i=Iudg_min,Iudg_max
      do j=Judg_min,Judg_max-1
        ull_dg(i,j,:,1) = ullu_dg(i,j,:,1)+ullv_dg(i-1,j+1,:,1)
      enddo
      ull_dg(i,Judg_max,:,1) = ullu_dg(i,Judg_max,:,1)+ullv_dg(i-1,Judg_max,:,1)
    enddo
    do j=Jvdg_min,Jvdg_max
      do i=Ivdg_min,Ivdg_max-1
        vll_dg(i,j,:,1) = vllv_dg(i,j,:,1)+vllu_dg(i+1,j-1,:,1)
      enddo
      vll_dg(Ivdg_max,j,:,1) = vllv_dg(Ivdg_max,j,:,1)+vllu_dg(Ivdg_max,j-1,:,1)
    enddo
#else
  do i=Iudg_min,Iudg_max
    do j=Judg_min,Judg_max
      ull_dg(i,j,:,1) = ullu_dg(i,j,:,1)+ullv_dg(i,j,:,1)
    enddo
  enddo
  do j=Jvdg_min,Jvdg_max
    do i=Ivdg_min,Ivdg_max
      vll_dg(i,j,:,1) = vllv_dg(i,j,:,1)+vllu_dg(i,j,:,1)
    enddo
  enddo
#endif  
  
    write(*,*) 'Linear Interporation: u'
    call interp3D_grid3_2(                                  &
            Iudg_min, Iudg_max, Judg_min, Judg_max, 1, Ndg  &
          , ull_dg                                          &
          , 1, L, 0, M , 1, N                               &
          , Id_cnt3Du, w_cnt3Du                             &
          , ull ) 
  
    write(*,*) 'Linear Interporation: v'
  
    call interp3D_grid3_2(                                  &
            Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, 1, Ndg  &
          , vll_dg                                          &
          , 0, L, 1, M , 1, N                               &
          , Id_cnt3Dv, w_cnt3Dv                             &
          , vll  ) 
    
  !- convert lat lon coordinate to ROMS refine coordinate --------------------------------
  
    do i=1,L
      do j=0,M
        uu(i,j,:,1) =  ull(i,j,:,1)*cosAu(i,j)
        vu(i,j,:,1) = -ull(i,j,:,1)*sinAu(i,j)
      enddo
    enddo
    do i=0,L
      do j=1,M
        uv(i,j,:,1) = vll(i,j,:,1)*sinAv(i,j)
        vv(i,j,:,1) = vll(i,j,:,1)*cosAv(i,j)
      enddo
    enddo
    do i=1,L
      do j=0,M-1
        u(i,j,:,1) = uu(i,j,:,1)+uv(i-1,j+1,:,1)
      enddo
      u(i,M,:,1) = uu(i,M,:,1)+uv(i-1,M,:,1)
    enddo
    do j=1,M
      do i=0,L-1
        v(i,j,:,1) = vv(i,j,:,1)+vu(i+1,j-1,:,1)
      enddo
      v(L,j,:,1) = vv(L,j,:,1)+vu(L,j-1,:,1)
    enddo
  
    if( romsvar(2)==1 ) then
      !- Write u v --------------------------------
      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ Nxu, Nyu, Nzr, 1 /)
      write(*,*)  'Write: u'
      call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
      call check( nf90_inq_varid(ncid2, 'u', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, u, start = start4D, count = count4D) )
      call check( nf90_close(ncid2) )
    endif
    if( romsvar(3)==1 ) then
      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ Nxv, Nyv, Nzr, 1 /)
      write(*,*)  'Write: v'
      call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
      call check( nf90_inq_varid(ncid2, 'v', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, v, start = start4D, count = count4D) )
      call check( nf90_close(ncid2) )
    endif
  
  !-Depth averaged velocity calculation --------------------------------
  
    if( romsvar(4)==1 ) then
      ubar(:,:,1)=0.0d0
      do i=1,L
        do j=0,M
          do k=1,N
            ubar(i,j,1) = ubar(i,j,1)+u(i,j,k,1)*dz_u(i,j,k)/hu(i,j)*umask(i,j)
          enddo
        enddo
      enddo
      start3D = (/ 1,  1,  1 /)
      count3D = (/ Nxu, Nyu, 1 /)
      write(*,*)  'Write: ubar'
      call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
      call check( nf90_inq_varid(ncid2, 'ubar', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, ubar, start = start3D, count = count3D) )
      call check( nf90_close(ncid2) )
    endif

    if( romsvar(5)==1 ) then   
      vbar(:,:,1)=0.0d0
      do i=0,L
        do j=1,M
          do k=1,N
            vbar(i,j,1) = vbar(i,j,1)+v(i,j,k,1)*dz_v(i,j,k)/hv(i,j)*vmask(i,j)
          enddo
        enddo
      enddo
      start3D = (/ 1,  1,  1 /)
      count3D = (/ Nxv, Nyv, 1 /)
      write(*,*)  'Write: vbar'
      call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
      call check( nf90_inq_varid(ncid2, 'vbar', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, vbar, start = start3D, count = count3D) )
      call check( nf90_close(ncid2) )
    endif
  
  endif
  
!- Tracer (temp, salt, etc.) --------------------------------
#if defined ROMS_MODEL
  i = 6
  j = 1
  do while ( i<= N_var )         
    if( romsvar(i)==0 ) then
      i=i+1
      cycle
    endif
#else
  do i=6, 7
    if( romsvar(i)==0 ) cycle
#endif
    varname = trim( VAR_NAME(i) )

#if defined ROMS_MODEL
    start4D = (/ Irdg_min+1, Jrdg_min+1, 1,      itime /)
    count4D = (/ Nxr_dg,     Nyr_dg,     Nzr_dg, 1     /)
    if( i==8  .or. i==9  .or. i==13 .or. i==14 .or. i==15 .or. i==16 .or.  &
        i==17 .or. i==21 .or. i==22 .or. i==23 .or. i==24 .or. i==26 .or.  &
        i==27 .or. i==28 .or. i==29 .or. i==30 .or. i==33 .or. i==34 .or.  &
        i==35 .or. i==36  ) then  ! mud_, sand_, etc...
      write(varnum,'(I2.2)') j
      varname = trim( VAR_NAME(i) )//varnum
      status = nf90_inq_varid(ncid, trim( varname ), var_id)
      if (status /= nf90_noerr) then
        i=i+1
        j = 1
        cycle
      endif
      write(*,*) 'Read: ', trim( varname )
      call check( nf90_get_var(ncid, var_id, t_dg, start4D, count4D)  )
      j=j+1
    else
      write(*,*) 'Read: ', trim( varname )
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_get_var(ncid, var_id, t_dg, start4D, count4D)  )
      i=i+1
    endif

#elif defined HYCOM_MODEL
    if( i == 6 ) then
      HY_VAR_NAME = 'water_temp'
    elseif( i == 7 ) then
      HY_VAR_NAME = 'salinity'
    else
      cycle
    endif
    start4D = (/ Irdg_min+1, Jrdg_min+1, 1,      itime /)
    count4D = (/ Nxr_dg,     Nyr_dg,     Nzr_dg, 1     /)
    call readNetCDF_4d_hycom( ncid, trim( HY_VAR_NAME )  &
          , Nxr_dg, Nyr_dg, Nzr_dg, 1, start4D, count4D, t_dg )

#elif defined JCOPE_MODEL
    JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(i) )//trim( JCOPE_sufix )
    write(*,*) 'Read: ', trim( JCOPE_data_file )
    call read_jcope_data3D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
          , Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, t_dg(:,:,:,1) )
# if !defined JCOPE_T
    if( i == 6 ) then
      t_dg(:,:,:,1) = t_dg(:,:,:,1) + 10.0d0
    elseif( i == 7 ) then
      t_dg(:,:,:,1) = t_dg(:,:,:,1) + 35.0d0
    endif
# endif
#endif

    write(*,*) 'Linear Interporation: ', trim( varname )
    call interp3D_grid3_2(                                &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, 1, Ndg  &
        , t_dg                                            &
        , 0, L, 0, M , 1, N                               &
        , Id_cnt3Dr, w_cnt3Dr                             &
        , t  ) 

    start4D = (/ 1,   1,   1,   1 /)
    count4D = (/ Nxr, Nyr, Nzr, 1 /)
    write(*,*)  'Write: ', trim( varname )
    call check( nf90_open(trim( INI_FILE ), NF90_WRITE, ncid2) )
    call check( nf90_inq_varid(ncid2, trim( varname ), var_id2) )
    call check( nf90_put_var(ncid2, var_id2, t, start = start4D, count = count4D) )
    call check( nf90_close(ncid2) )
  
  enddo

#if defined ROMS_MODEL
  call check( nf90_close(ncid) )
#endif
!  call check( nf90_close(ncid2) )

  write(*,*) 'FINISH!!'
      
      
END PROGRAM iniOCN2ROMS
      
