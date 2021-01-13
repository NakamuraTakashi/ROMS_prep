
!!!=== Copyright (c) 2014-2021 Takashi NAKAMURA  ===== 

PROGRAM bryOCN2ROMS
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
  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  integer :: Ryear, Rmonth, Rday
  integer :: itime,ibry
  integer :: mode
  character(256) :: GRID_FILE
  character(256) :: parent_grid
  integer :: parent_Imin, parent_Imax
  integer :: parent_Jmin, parent_Jmax
  integer :: refine_factor
  character(256) :: ROMS_HISFILE
#if defined ROMS_MODEL
  integer :: romsvar(N_var)
#else
  integer :: romsvar(7)
#endif 
  character(256) :: BRY_prefix
  integer :: SNWE(4)
  integer :: N_s_rho
  integer :: Nzr
  integer :: spherical
  integer :: Vtransform, Vstretching
  real(8) :: THETA_S , THETA_B, TCLINE, DCRIT
  integer :: Nzr_dg
  integer :: Vtransform_dg, Vstretching_dg
  integer :: LBri, UBri, LBrj, UBrj
  integer :: LBui, UBui, LBuj, UBuj
  integer :: LBvi, UBvi, LBvj, UBvj
 
  integer :: Nrbry,Nubry,Nvbry

  character(33) :: TIME_ATT  = "seconds since 2000-01-01 00:00:00"

  character(15) :: BRY_suffix   = "_20000101.00.nc"
  character(256) :: BRY_FILE
  
  real(8), allocatable :: time_all(:)
  real(8), allocatable :: bry_time(:)
  integer, allocatable :: iNCt(:)
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
  
  real(8), allocatable :: zeta(:,:,:)    ! free-surface (meter)
  real(8), allocatable :: ubar(:,:,:)    ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar(:,:,:)    ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u(:,:,:,:)     ! u-momentum component (meter second-1)
  real(8), allocatable :: v(:,:,:,:)     ! v-momentum component (meter second-1)
  real(8), allocatable :: t(:,:,:,:)     ! tracer 
  real(8), allocatable :: ull(:,:,:,:)   ! u-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: vll(:,:,:,:)   ! v-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: uu(:,:,:,:)
  real(8), allocatable :: uv(:,:,:,:)
  real(8), allocatable :: vu(:,:,:,:)
  real(8), allocatable :: vv(:,:,:,:)
  real(8), allocatable :: zeta_bry(:,:)  ! free-surface (meter)
  real(8), allocatable :: ubar_bry(:,:)  ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar_bry(:,:)  ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u_bry(:,:,:)   ! u-momentum component (meter second-1)
  real(8), allocatable :: v_bry(:,:,:)   ! v-momentum component (meter second-1)
  real(8), allocatable :: t_bry(:,:,:)   ! tracer 

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

  real(8) :: ocean_time(1)                ! Ocean time (sec)
  real(8), allocatable :: zeta_dg(:,:,:)  ! free-surface (meter)
  real(8), allocatable :: ubar_dg(:,:,:)  ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar_dg(:,:,:)  ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u_dg(:,:,:,:)   ! u-momentum component (meter second-1)
  real(8), allocatable :: v_dg(:,:,:,:)   ! v-momentum component (meter second-1)
  real(8), allocatable :: t_dg(:,:,:,:)   ! tracer 
  real(8), allocatable :: ull_dg(:,:,:,:) ! u-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: vll_dg(:,:,:,:) ! v-momentum component on lat lon coordinate (meter second-1)
  real(8), allocatable :: ullu_dg(:,:,:,:) 
  real(8), allocatable :: ullv_dg(:,:,:,:) 
  real(8), allocatable :: vllu_dg(:,:,:,:) 
  real(8), allocatable :: vllv_dg(:,:,:,:) 
#if defined WET_DRY
  real(8), allocatable :: umask_wet_dg(:,:,:)
  real(8), allocatable :: vmask_wet_dg(:,:,:)
#endif
  integer, allocatable :: ID_cnt2Dr(:,:)
  real(8), allocatable :: w_cnt2Dr (:,:)
  integer, allocatable :: ID_cnt2Du(:,:)
  real(8), allocatable :: w_cnt2Du (:,:)
  integer, allocatable :: ID_cnt2Dv(:,:)
  real(8), allocatable :: w_cnt2Dv (:,:)
  integer, allocatable :: ID_cnt3Dr(:,:)
  real(8), allocatable :: w_cnt3Dr (:,:)
  integer, allocatable :: ID_cnt3Du(:,:)
  real(8), allocatable :: w_cnt3Du (:,:)
  integer, allocatable :: ID_cnt3Dv(:,:)
  real(8), allocatable :: w_cnt3Dv (:,:)

  integer :: i,j,k
  real(8) :: d_lat,d_lon
  real(8) :: cff
  
  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh
  character(11) :: YYYYMMDDpHH
  character(12) :: YYYYMMDDhhmm
  character(8) :: YYYYMMDD

  real(8) :: d_jdate
  integer :: N_days
  integer :: jdate_Start, jdate_End, jdate_Ref
  real(8) :: d_jdate_Start, d_jdate_Ref
  
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
  integer, allocatable :: idt(:)

  integer :: ncid,var_id
  integer :: ncid2,var_id2
  integer :: iNC, iNCm
  integer :: ItS, ItE
  character(256) :: varname
  character(256) :: bryname
  character(2) :: varnum
  integer :: status, access

  !
#if defined HYCOM_MODEL

  TYPE T_NC
    real(8), pointer :: time_all(:)
    real(8), pointer :: lat(:), lon(:), depth(:)
    integer :: Nxr_dg, Nyr_dg, Nzr_dg
    integer :: Nt
    integer :: ItS, ItE
  END TYPE T_NC
  TYPE (T_NC) :: NC(NCnum)

  real(8), allocatable :: time2(:)

  character(256) :: HY_VAR_NAME
  integer :: iNCs, iNCe
  integer :: jdate_20000101
  real(8) :: d_jdate_20000101
#endif
#if defined JCOPE_MODEL
  character(256) :: JCOPE_info_dir, JCOPE_data_dir
  character(256) :: JCOPE_data_file
  real(8), allocatable :: jcope_data(:,:,:)
  real(8), allocatable :: z(:,:,:), zz(:,:,:), dz(:,:,:)
  real(8), allocatable :: lat(:), lon(:)
  character(20) :: JCOPE_sufix(10000)
  character(20) :: JCOPE_sufix_tmp
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: idays,ihours,ijdate,Sjdate
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
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/bry/BRY_prefix, SNWE
  namelist/hcoord/spherical
  namelist/zcoord/N_s_rho
  namelist/zcoord/Vtransform, Vstretching
  namelist/zcoord/THETA_S, THETA_B, TCLINE, DCRIT

  ! Read parameters in namelist file
  
  read (*, nml=grd)
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
  read (*, nml=sdate)
  rewind(5)
  read (*, nml=edate)
  rewind(5)
  read (*, nml=refdate)
  rewind(5)
  read (*, nml=bry)
  rewind(5)
  read (*, nml=hcoord)
  rewind(5)
  read (*, nml=zcoord)
           
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
        z_u(i,j,k)= (z_r(i-1,j,k)+z_r(i,j,k))*0.5d0
        dz_u(i,j,k) = ((z_w(i-1,j,k-1)+z_w(i,j,k-1))-(z_w(i-1,j,k)+z_w(i,j,k)))*0.5d0
      enddo
    enddo
  enddo
  do i=0,L
    do j=1,M
      do k=1,N
        z_v(i,j,k)= (z_r(i,j-1,k)+z_r(i,j,k))*0.5d0
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

!======= Set starting time and Ending time ========================

  call ndays(Emonth, Eday, Eyear, Smonth, Sday, Syear, N_days)
  call jd(Syear, Smonth, Sday, jdate_Start)

  jdate_End = jdate_Start + N_days
  
  write(*,*) jdate_Start, jdate_End, N_days

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

  do i=0,Ldg
    do j=0,Mdg
      if ( h_dg(i,j)<0.0d0 ) then
        rmask_dg(i,j) = 0.0d0
      end if
    end do
  end do
    
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

! Check time

  ItS = 1
  ItE = Nt
  do i=2,Nt
    d_jdate = d_jdate_Ref + time_all(i)/86400.0d0
    if( d_jdate > dble(jdate_Start) ) then
      ItS = i-1
      exit
    endif
  enddo
  do i=ItS,Nt
    d_jdate = d_jdate_Ref + time_all(i)/86400.0d0
    if( d_jdate >= dble(jdate_End) ) then
      ItE = i-1
      exit
    endif
  enddo

  Nt = ItE - ItS + 1
  
  allocate(bry_time(Nt))
  allocate(idt(Nt))
  allocate(iNCt(Nt))

  bry_time(:) = time_all(ItS:ItE)
  
  do i=1, Nt
    idt(i) = ItS + i-1
    iNCt(i) = 1
  enddo

! Extract grid data

  call set_depth ( Nxr_dg, Nyr_dg, Nzr_dg, h_dg                      &
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
  
  do iNC=NCnum,1,-1
    do i=NC(iNC)%Nt,1,-1
      d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24.0d0
      if(d_jdate < dble(jdate_End)) then
        write(*,*) '*** FOUND: Ending point @ HYCOM_FILE',iNC
        NC(iNC)%ItE=i
        exit
      endif
    end do
  end do
  write(*,*) NC(:)%ItE 
  
  do iNC=1,NCnum
    do i=2,NC(iNC)%ItE
      d_jdate=d_jdate_20000101+NC(iNC)%time_all(i)/24.0d0
      if(d_jdate>dble(jdate_Start)) then
        write(*,*) '*** FOUND: Starting point @ HYCOM_FILE',iNC
        NC(iNC)%ItS=i-1
        exit
      endif
    end do
  end do
  write(*,*) NC(:)%ItS 

  Nt = 0
  iNCs = NCnum
  iNCe = 1

  do iNC=1,NCnum
    if(NC(iNC)%ItS==-1) then
      cycle
    end if
    Nt = Nt + NC(iNC)%ItE - NC(iNC)%ItS + 1
    iNCs = min(iNCs,iNC)
    iNCe = max(iNCe,iNC)
  enddo

  allocate(bry_time(Nt))
  allocate(idt(Nt))
  allocate(iNCt(Nt))

  i=1
  do iNC=iNCs,iNCe
    do k=0,NC(iNC)%ItE-NC(iNC)%ItS
      iNCt(i+k) = iNC
      idt(i+k)  = NC(iNC)%ItS + k
    enddo
    j=i+NC(iNC)%ItE-NC(iNC)%ItS
    bry_time(i:j) = NC(iNC)%time_all( NC(iNC)%ItS : NC(iNC)%ItE )
    i=i+j
  enddo

  bry_time(:) = ( bry_time(:) + ( d_jdate_Ref - d_jdate_20000101 )*24.0d0)*3600.0d0 ! hour -> sec


  !---- Read HYCOM netCDF grid cooredinate --------------------------------
  do iNC=iNCs,iNCe
    write(*,*) "OPEN: ", HYCOM_FILE(iNC)
    call try_nf_open(HYCOM_FILE(iNC), nf90_nowrite, ncid)
    ! Get dimension data
    call get_dimension(ncid, 'lat', Nyr_dg)
    call get_dimension(ncid, 'lon', Nxr_dg)
    call get_dimension(ncid, 'depth', Nzr_dg)
    write(*,*) Nxr_dg, Nyr_dg, Nzr_dg
    ! Allocate variable
    NC(iNC)%Nyr_dg = Nyr_dg
    NC(iNC)%Nxr_dg = Nxr_dg
    NC(iNC)%Nzr_dg = Nzr_dg

    allocate( NC(iNC)%lat(Nyr_dg) )
    allocate( NC(iNC)%lon(Nxr_dg) )
    allocate( NC(iNC)%depth(Nzr_dg) )

    call readNetCDF_1d(ncid, 'lat', Nyr_dg, NC(iNC)%lat)
    call readNetCDF_1d(ncid, 'lon', Nxr_dg, NC(iNC)%lon)
    call readNetCDF_1d(ncid, 'depth', Nzr_dg, NC(iNC)%depth)
    call check( nf90_close(ncid) )
    write(*,*) "CLOSE: ", HYCOM_FILE(iNC)
  enddo

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

! ---------------- bry_time check ----------------------------------

  allocate( bry_time(10000) ) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Nt = 0
  itime = 1
  ihours = 0
  iyear = Syear
  imonth = Smonth
  iday = Sday
  call jd(iyear, imonth, iday, Sjdate)

  do
    ihour = mod(ihours,24)
    ijdate = Sjdate + int(ihours/24)
    call cdate( ijdate, iyear, imonth, iday )
    ! Check end date
    if(iyear==Eyear .and. imonth==Emonth .and. iday==Eday) then
      write(*,*) "Total time step: ", Nt
      write(*,*) "Start file: ", trim( JCOPE_data_dir )//trim( JCOPE_prefix(4) )//trim( JCOPE_sufix(1) )
      write(*,*) "End file:   ", trim( JCOPE_data_dir )//trim( JCOPE_prefix(4) )//trim( JCOPE_sufix(Nt) )
      write(*,*) bry_time(1), bry_time(Nt)
      exit
    endif

    write (YYYY, "(I4.4)") iyear
    write (MM, "(I2.2)") imonth
    write (DD, "(I2.2)") iday
    write (hh, "(I2.2)") ihour
    call jd(iyear, imonth, iday, ijdate)

#  if defined JCOPE_T
    ihours = ihours + 1  !!! Files exist 1 hourly interval
    JCOPE_sufix_tmp = YYYY//MM//DD//hh//"00.bin"
#  else
    ihours = ihours + 24  !!! Files exist daily interval
    JCOPE_sufix_tmp = YYYY//MM//DD
#  endif

    JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(1) )//trim( JCOPE_sufix_tmp )
    status = access( trim( JCOPE_data_file ), ' ' )
    if ( status /= 0 ) cycle 
    Nt = Nt + 1
    JCOPE_sufix(Nt) = JCOPE_sufix_tmp  
    bry_time(Nt) = ( dble( ijdate ) + dble( ihour )/24.0d0 - d_jdate_Ref ) * 86400.0d0  ! sec
  enddo

!  allocate(idt(Nt))
  allocate(iNCt(Nt))
  
  do i=1, Nt
!    idt(i) = ItS + i-1
    iNCt(i) = 1
  enddo

#else
  write(*,*) "************ Please choose a OCEAN MODEL ************"
  Stop
#endif

! ----------------------------------------------------
!-Read ROMS ocean_his netCDF file --------------------------------
  ocean_time(1) = bry_time(1)
!-Boundary condition netcdf file name -------------------------------------------------
  call oceantime2cdate(ocean_time(1), Ryear, Rmonth, Rday, YYYYMMDDpHH)
  BRY_suffix(2:12)=YYYYMMDDpHH
  BRY_FILE = trim( BRY_prefix )//BRY_suffix

!-Create the ROMS initial conditions netCDF file --------------------------------
#if defined ROMS_MODEL     
  call createNetCDFbry2(   trim( ROMS_HISFILE ), trim( BRY_FILE )   &
        , TIME_ATT , Nxr, Nyr, Nzr, 1 ,romsvar ,SNWE   )
#else
  call createNetCDFbry(  trim( BRY_FILE ), TIME_ATT       &
        , Nxr, Nyr, Nzr, 1 ,romsvar ,SNWE   )
#endif
!-Write ROMS boundary conditions netCDF file --------------------------------
 
  call check( nf90_open(trim( BRY_FILE ), NF90_WRITE, ncid2) )
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

  start1D = (/ 1 /)
  count1D = (/ Nt /)
  call check( nf90_inq_varid(ncid2, 'bry_time', var_id) )
  call check( nf90_put_var(ncid2, var_id, bry_time, start = start1D, count = count1D) )
  call check( nf90_close(ncid2) )      


!===== LOOP START ==========================================================

!-Write ROMS initial conditions netCDF file --------------------------------
  do ibry=1,4
    if( SNWE(ibry)==0 ) cycle

    write(*,*) '***************** Start extracting '//trim(BRY_NAME(ibry))//' boundary *****************'

    if (ibry == 1) then     ! south
      LBri = 0
      UBri = L
      LBrj = 0
      UBrj = 0
      LBui = 1
      UBui = L
      LBuj = 0
      UBuj = 1
      LBvi = 0
      UBvi = L
      LBvj = 1
      UBvj = 1
      Nrbry = UBri-LBri+1
      Nubry = UBui-LBui+1
      Nvbry = UBvi-LBvi+1
    elseif (ibry == 2) then ! North
      LBri = 0
      UBri = L
      LBrj = M
      UBrj = M
      LBui = 1
      UBui = L
      LBuj = M-1
      UBuj = M
      LBvi = 0
      UBvi = L
      LBvj = M
      UBvj = M
      Nrbry = UBri-LBri+1
      Nubry = UBui-LBui+1
      Nvbry = UBvi-LBvi+1
    elseif (ibry == 3) then ! West
      LBri = 0
      UBri = 0
      LBrj = 0
      UBrj = M
      LBui = 1
      UBui = 1
      LBuj = 0
      UBuj = M
      LBvi = 0
      UBvi = 1
      LBvj = 1
      UBvj = M
      Nrbry = UBrj-LBrj+1
      Nubry = UBuj-LBuj+1
      Nvbry = UBvj-LBvj+1
    else                    ! East
      LBri = L
      UBri = L
      LBrj = 0
      UBrj = M
      LBui = L
      UBui = L
      LBuj = 0
      UBuj = M
      LBvi = L-1
      UBvi = L
      LBvj = 1
      UBvj = M
      Nrbry = UBrj-LBrj+1
      Nubry = UBuj-LBuj+1
      Nvbry = UBvj-LBvj+1
    endif
    Nxr = UBri-LBri+1
    Nxu = UBui-LBui+1
    Nxv = UBvi-LBvi+1
    Nyr = UBrj-LBrj+1
    Nyu = UBuj-LBuj+1
    Nyv = UBvj-LBvj+1

    allocate( zeta(LBri:UBri, LBrj:UBrj, 1) )
    allocate( ubar(LBui:UBui, LBuj:UBuj, 1) )
    allocate( vbar(LBvi:UBvi, LBvj:UBvj, 1) )
    allocate( u(LBui:UBui, LBuj:UBuj, 1:N, 1) )
    allocate( v(LBvi:UBvi, LBvj:UBvj, 1:N, 1) )
    allocate( ull(LBui:UBui, LBuj:UBuj, 1:N, 1) )
    allocate( uu (LBui:UBui, LBuj:UBuj, 1:N, 1) )
    allocate( vu (LBui:UBui, LBuj:UBuj, 1:N, 1) )
    allocate( vll(LBvi:UBvi, LBvj:UBvj, 1:N, 1) )
    allocate( uv (LBvi:UBvi, LBvj:UBvj, 1:N, 1) )
    allocate( vv (LBvi:UBvi, LBvj:UBvj, 1:N, 1) )
    allocate( t(LBri:UBri, LBrj:UBrj, 1:N, 1) )
    allocate( zeta_bry(Nrbry, 1) )
    allocate( ubar_bry(Nubry, 1) )
    allocate( vbar_bry(Nvbry, 1) )
    allocate( u_bry(Nubry, 1:N, 1) )
    allocate( v_bry(Nvbry, 1:N, 1) )
    allocate( t_bry(Nrbry, 1:N, 1) )
    
    write(*,*) "*** Calculate weight factor ***"

    allocate( ID_cnt2Dr(4, Nxr*Nyr) )
    allocate( w_cnt2Dr (4, Nxr*Nyr) )
    allocate( ID_cnt3Dr(6, Nxr*Nyr*Nzr) )
    allocate( w_cnt3Dr (8, Nxr*Nyr*Nzr) )
    allocate( ID_cnt2Du(4, Nxu*Nyu) )
    allocate( w_cnt2Du (4, Nxu*Nyu) )
    allocate( ID_cnt2Dv(4, Nxv*Nyv) )
    allocate( w_cnt2Dv (4, Nxv*Nyv) )
    allocate( ID_cnt3Du(6, Nxu*Nyu*Nzr) )
    allocate( w_cnt3Du (8, Nxu*Nyu*Nzr) )
    allocate( ID_cnt3Dv(6, Nxv*Nyv*Nzr) )
    allocate( w_cnt3Dv (8, Nxv*Nyv*Nzr) )

    iNC = 0

    do itime=1,Nt
      iNCm = iNC
      iNC  = iNCt(itime)

      if(iNCm < iNC) then
#if defined HYCOM_MODEL
  !---- Load HYCOM netCDF grid cooredinate --------------------------------

        Nyr_dg = NC(iNC)%Nyr_dg
        Nxr_dg = NC(iNC)%Nxr_dg
        Nzr_dg = NC(iNC)%Nzr_dg
      
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
        allocate( z_u_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
        allocate( z_v_dg(0:Ldg, 0:Mdg, 1:Ndg ) )
      
        do i=0,Ldg
          do j=0,Mdg
            latr_dg(i,j) = NC(iNC)%lat(j+1)
            lonr_dg(i,j) = NC(iNC)%lon(i+1)
          enddo
        enddo
        latu_dg(:,:) = latr_dg(:,:)
        lonu_dg(:,:) = lonr_dg(:,:)
        latv_dg(:,:) = latr_dg(:,:)
        lonv_dg(:,:) = lonr_dg(:,:)
        do i=1,Ldg
          do j=0,Mdg
            do k=1,Ndg
              z_r_dg(i,j,k)= NC(iNC)%depth(Ndg-k+1)
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
  
        write(*,*) "OPEN: ", HYCOM_FILE(iNC)
        call try_nf_open(HYCOM_FILE(iNC), nf90_nowrite, ncid)      
#endif            
  
        write(*,*) "Seek rho point donor IJ range"
        call seek_IJrange(                                 &
              0, Ldg, 0, Mdg, lonr_dg, latr_dg             & 
            , LBri, UBri, LBrj, UBrj                       &
            , lonr(LBri:UBri,LBrj:UBrj)                    &
            , latr(LBri:UBri,LBrj:UBrj)                    &
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
        allocate( t_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, 1:Ndg, 1) )
        allocate( ull_dg (Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
        allocate( ullu_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
        allocate( vllu_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1:Ndg, 1) )
        allocate( vll_dg (Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
        allocate( vllv_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
        allocate( ullv_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1:Ndg, 1) )
#if defined WET_DRY
        allocate( umask_wet_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, 1) )
        allocate( vmask_wet_dg(Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, 1) )
#endif

#if defined HYCOM_MODEL
        write(*,*) 'Read surf_el, and create land mask'
        start3D = (/ Irdg_min+1, Jrdg_min+1, NC(iNC)%ItS /)
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
        call weight2D_grid3_2(                              &
              Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
            , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
            , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
            , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
            , LBri, UBri, LBrj, UBrj                        &
            , lonr(LBri:UBri,LBrj:UBrj)                     &
            , latr(LBri:UBri,LBrj:UBrj)                     &
            , ID_cnt2Dr, w_cnt2Dr )
     
        write(*,*) "Calculate 2D u point weight factor"
        call weight2D_grid3_2(                              &
              Iudg_min, Iudg_max, Judg_min, Judg_max        &
            , lonu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
            , latu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
            , umask_dg(Iudg_min:Iudg_max,Judg_min:Judg_max) & 
            , LBui, UBui, LBuj, UBuj                        &
            , lonu(LBui:UBui,LBuj:UBuj)                     &
            , latu(LBui:UBui,LBuj:UBuj)                     &
            , ID_cnt2Du, w_cnt2Du )
    
        write(*,*) "Calculate 2D v point weight factor"
        call weight2D_grid3_2(                              &
              Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max        &
            , lonv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
            , latv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
            , vmask_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max) & 
            , LBvi, UBvi, LBvj, UBvj                        &
            , lonv(LBvi:UBvi,LBvj:UBvj)                     &
            , latv(LBvi:UBvi,LBvj:UBvj)                     &
            , ID_cnt2Dv, w_cnt2Dv )
    
        write(*,*) "Calculate 3D rho point weight factor"
        call weight3D_grid3_2(                                  &
              Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, 1, Ndg    &
            , z_r_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max,1:Ndg) &
            , LBri, UBri, LBrj, UBrj, 1, N                      &
            , z_r(LBri:UBri,LBrj:UBrj, 1:N)                     &
            , ID_cnt2Dr, w_cnt2Dr                               &
            , ID_cnt3Dr, w_cnt3Dr )
    
        write(*,*) "Calculate 3D u point weight factor"
        call weight3D_grid3_2(                                  &
              Iudg_min, Iudg_max, Judg_min, Judg_max, 1, Ndg    &
            , z_u_dg(Iudg_min:Iudg_max,Judg_min:Judg_max,1:Ndg) &
            , LBui, UBui, LBuj, UBuj, 1, N                      &
            , z_u(LBui:UBui,LBuj:UBuj, 1:N)                     &
            , ID_cnt2Du, w_cnt2Du                               &
            , ID_cnt3Du, w_cnt3Du )
    
        write(*,*) "Calculate 3D v point weight factor"
        call weight3D_grid3_2(                                  &
              Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, 1, Ndg    &
            , z_v_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max,1:Ndg) &
            , LBvi, UBvi, LBvj, UBvj, 1, N                      &
            , z_v(LBvi:UBvi,LBvj:UBvj, 1:N)                     &
            , ID_cnt2Dv, w_cnt2Dv                               &
            , ID_cnt3Dv, w_cnt3Dv )
  
        write(*,*) "***************************************************************"

      endif

      ocean_time(1) = bry_time(itime)

      call oceantime2cdate(ocean_time(1), Ryear, Rmonth, Rday, YYYYMMDDpHH)
      Write(*,*) 'Time = ',YYYYMMDDpHH, ',  bry = ',trim( BRY_NAME(ibry) )
!      Write(*,*) idt(itime) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if( romsvar(1)==1 ) then
        !- zeta --------------------------------
#if defined ROMS_MODEL
        write(*,*) 'Read: zeta'
        start3D = (/ Irdg_min+1, Jrdg_min+1, idt(itime) /)
        count3D = (/ Nxr_dg,     Nyr_dg,     1          /)
        call check( nf90_inq_varid(ncid, 'zeta', var_id) )
        call check( nf90_get_var(ncid, var_id, zeta_dg, start3D, count3D)  )
#elif defined HYCOM_MODEL
        write(*,*) 'Read: '//HYCOM_FILE(iNC)
        start3D = (/ Irdg_min+1, Jrdg_min+1, idt(itime) /)
        count3D = (/ Nxr_dg,     Nyr_dg,     1          /)
        call readNetCDF_3d_hycom( ncid, 'surf_el'      &
              , Nxr_dg, Nyr_dg, 1, start3D, count3D, zeta_dg )
#elif defined JCOPE_MODEL
       JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(1) )//trim( JCOPE_sufix(itime) )
       write(*,*) 'Read: ', trim( JCOPE_data_file )
       call read_jcope_data2D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
             , Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, zeta_dg )
#endif
      
        write(*,*) 'Linear Interporation: zeta'
        call interp2D_grid3_2(                                              &
                Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, zeta_dg             &
              , LBri, UBri, LBrj, UBrj                                      &
              , Id_cnt2Dr, w_cnt2Dr                                         &
              , zeta  ) 
        
        if(ibry == 1) then     ! South
          zeta_bry(:,1) = zeta(LBri:UBri, 0, 1)
        elseif(ibry == 2) then ! North
          zeta_bry(:,1) = zeta(LBri:UBri, M, 1)
        elseif(ibry == 3) then ! West
          zeta_bry(:,1) = zeta(0, LBrj:UBrj, 1)
        else                   ! East
          zeta_bry(:,1) = zeta(L, LBrj:UBrj, 1)
        endif
    
        start2D = (/ 1,     itime /)
        count2D = (/ Nrbry, 1     /)
        write(*,*)  'Write: ', 'zeta_'//trim(BRY_NAME(ibry))
        call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
        call check( nf90_inq_varid(ncid2, 'zeta_'//trim(BRY_NAME(ibry)), var_id2) )
        call check( nf90_put_var(ncid2, var_id2, zeta_bry, start = start2D, count = count2D) )
        call check( nf90_close(ncid2) )      
    
      endif

      if( romsvar(2)==1 .or. romsvar(3)==1 .or.         &
          romsvar(4)==1 .or. romsvar(5)==1      ) then
    !- u --------------------------------
    
#if defined ROMS_MODEL
        write(*,*) 'Read: u'
        start4D = (/ Iudg_min, Judg_min+1, 1,      idt(itime) /)
        count4D = (/ Nxu_dg,     Nyu_dg,   Nzr_dg, 1          /)
        call check( nf90_inq_varid(ncid, 'u', var_id) )
        call check( nf90_get_var(ncid, var_id, u_dg, start4D, count4D)  )
#elif defined HYCOM_MODEL
        write(*,*) 'Read: '//HYCOM_FILE(iNC)
        start4D = (/ Iudg_min+1, Judg_min+1, 1,      idt(itime) /)
        count4D = (/ Nxu_dg,     Nyu_dg,     Nzr_dg, 1          /)
# if defined FAST_READ
        call readNetCDF_4d_hycom_fast( ncid, 'water_u'             &
              , Nxu_dg, Nyu_dg, Nzr_dg, 1, start4D, count4D, u_dg )
# else
        call readNetCDF_4d_hycom( ncid, 'water_u'                  &
              , Nxu_dg, Nyu_dg, Nzr_dg, 1, start4D, count4D, u_dg )
# endif
#elif defined JCOPE_MODEL
        JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(2) )//trim( JCOPE_sufix(itime) )
        write(*,*) 'Read: ', trim( JCOPE_data_file )
        call read_jcope_data3D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
              , Iudg_min, Iudg_max, Judg_min, Judg_max, u_dg(:,:,:,1) )
#endif
    
#if defined WET_DRY
        start3D = (/ Iudg_min, Judg_min+1, idt(itime) /)
        count3D = (/ Nxu_dg,   Nyu_dg,     1          /)
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
        start4D = (/ Ivdg_min+1, Jvdg_min, 1,      idt(itime) /)
        count4D = (/ Nxv_dg,     Nyv_dg,   Nzr_dg, 1          /)
        call check( nf90_inq_varid(ncid, 'v', var_id) )
        call check( nf90_get_var(ncid, var_id, v_dg, start4D, count4D)  )
#elif defined HYCOM_MODEL
        write(*,*) 'Read: '//HYCOM_FILE(iNC)
        start4D = (/ Ivdg_min+1, Jvdg_min+1, 1,    idt(itime) /)
        count4D = (/ Nxv_dg,     Nyv_dg,   Nzr_dg, 1          /)
# if defined FAST_READ
        call readNetCDF_4d_hycom_fast( ncid, 'water_v'             &
              , Nxv_dg, Nyv_dg, Nzr_dg, 1, start4D, count4D, v_dg )
# else
        call readNetCDF_4d_hycom( ncid, 'water_v'                  &
              , Nxv_dg, Nyv_dg, Nzr_dg, 1, start4D, count4D, v_dg )
# endif
#elif defined JCOPE_MODEL
        JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(3) )//trim( JCOPE_sufix(itime) )
        write(*,*) 'Read: ', trim( JCOPE_data_file )
        call read_jcope_data3D( trim( JCOPE_data_file ), 0, Ldg, 0, Mdg     &
              , Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, v_dg(:,:,:,1) )
#endif    
#if defined WET_DRY
        start3D = (/ Ivdg_min+1, Jvdg_min, idt(itime) /)
        count3D = (/ Nxv_dg,     Nyv_dg,   1          /)
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
              , LBui, UBui, LBuj, UBuj, 1, N                    &
              , Id_cnt3Du, w_cnt3Du                             &
              , ull )
        
        write(*,*) 'Linear Interporation: v'
        call interp3D_grid3_2(                                  &
                Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, 1, Ndg  &
              , vll_dg                                          &
              , LBvi, UBvi, LBvj, UBvj, 1, N                    &
              , Id_cnt3Dv, w_cnt3Dv                             &
              , vll  ) 
        
        ! convert lat lon coordinate to ROMS refine coordinate --------------------------------
    
        do i=LBui,UBui
          do j=LBuj,UBuj
            uu(i,j,:,1) =  ull(i,j,:,1)*cosAu(i,j)
            vu(i,j,:,1) = -ull(i,j,:,1)*sinAu(i,j)
          enddo
        enddo
        do i=LBvi,UBvi
          do j=LBvj,UBvj
            uv(i,j,:,1) = vll(i,j,:,1)*sinAv(i,j)
            vv(i,j,:,1) = vll(i,j,:,1)*cosAv(i,j)
          enddo
        enddo
        do i=LBui,UBui
          do j=LBuj,UBuj-1
            u(i,j,:,1) = uu(i,j,:,1)+uv(i-1,j+1,:,1)
          enddo
!          write(*,*) 'DEBUG 1', i,j, UBuj
          u(i,UBuj,:,1) = uu(i,UBuj,:,1)+uv(i-1,UBuj,:,1)
        enddo
        do j=LBvj,UBvj
          do i=LBvi,UBvi-1
            v(i,j,:,1) = vv(i,j,:,1)+vu(i+1,j-1,:,1)
          enddo
          v(UBvi,j,:,1) = vv(UBvi,j,:,1)+vu(UBvi,j-1,:,1)
        enddo
  
        if( romsvar(2)==1 ) then
        ! Write u v --------------------------------
            
          if(ibry == 1) then     ! South
            u_bry(:,:,1) = u(LBui:UBui, 0, 1:N, 1)
          elseif(ibry == 2) then ! North
            u_bry(:,:,1) = u(LBui:UBui, M, 1:N, 1)
          elseif(ibry == 3) then ! West
            u_bry(:,:,1) = u(1, LBuj:UBuj, 1:N, 1)
          else                   ! East
            u_bry(:,:,1) = u(L, LBuj:UBuj, 1:N, 1)
          endif
          start3D = (/ 1,     1,   itime /)
          count3D = (/ Nubry, Nzr, 1     /)
          write(*,*)  'Write: ', 'u_'//trim(BRY_NAME(ibry))
          call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
          call check( nf90_inq_varid(ncid2, 'u_'//trim(BRY_NAME(ibry)), var_id2) )
          call check( nf90_put_var(ncid2, var_id2, u_bry, start = start3D, count = count3D) )
          call check( nf90_close(ncid2) )      
        endif
  
        if( romsvar(3)==1 ) then
          if(ibry == 1) then     ! South
            v_bry(:,:,1) = v(LBvi:UBvi, 1, 1:N, 1)
          elseif(ibry == 2) then ! North
            v_bry(:,:,1) = v(LBvi:UBvi, M, 1:N, 1)
          elseif(ibry == 3) then ! West
            v_bry(:,:,1) = v(0, LBvj:UBvj, 1:N, 1)
          else                   ! East
            v_bry(:,:,1) = v(L, LBvj:UBvj, 1:N, 1)
          endif
          start3D = (/ 1,    1,    itime /)
          count3D = (/ Nvbry, Nzr, 1     /)
          write(*,*)  'Write: ', 'v_'//trim(BRY_NAME(ibry))
          call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
          call check( nf90_inq_varid(ncid2, 'v_'//trim(BRY_NAME(ibry)), var_id2) )
          call check( nf90_put_var(ncid2, var_id2, v_bry, start = start3D, count = count3D) )
          call check( nf90_close(ncid2) )      
      
        endif
      
    ! Depth averaged velocity calculation --------------------------------
        if( romsvar(4)==1 ) then
          ubar(:,:,1)=0.0d0
          do i=LBui, UBui
            do j=LBuj, UBuj
              do k=1,N
                ubar(i,j,1) = ubar(i,j,1)+u(i,j,k,1)*dz_u(i,j,k)/hu(i,j)*umask(i,j)
              enddo
            enddo
          enddo
    
          if(ibry == 1) then     ! South
            ubar_bry(:,1) = ubar(LBui:UBui, 0, 1)
          elseif(ibry == 2) then ! North
            ubar_bry(:,1) = ubar(LBui:UBui, M, 1)
          elseif(ibry == 3) then ! West
            ubar_bry(:,1) = ubar(1, LBuj:UBuj, 1)
          else                   ! East
            ubar_bry(:,1) = ubar(L, LBuj:UBuj, 1)
          endif      
          start2D = (/ 1,     itime /)
          count2D = (/ Nubry, 1     /)
          write(*,*)  'Write: ', 'ubar_'//trim(BRY_NAME(ibry))
          call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
          call check( nf90_inq_varid(ncid2, 'ubar_'//trim(BRY_NAME(ibry)), var_id2) )
          call check( nf90_put_var(ncid2, var_id2, ubar_bry, start = start2D, count = count2D) )
          call check( nf90_close(ncid2) )      
          
        endif
          
        if( romsvar(5)==1 ) then 
          vbar(:,:,1)=0.0d0
          do i=LBvi, UBvi
            do j=LBvj, UBvj
              do k=1,N
                vbar(i,j,1) = vbar(i,j,1)+v(i,j,k,1)*dz_v(i,j,k)/hv(i,j)*vmask(i,j)
              enddo
            enddo
          enddo
    
          if(ibry == 1) then     ! South
            vbar_bry(:,1) = vbar(LBvi:UBvi, 1, 1)
          elseif(ibry == 2) then ! North
            vbar_bry(:,1) = vbar(LBvi:UBvi, M, 1)
          elseif(ibry == 3) then ! West
            vbar_bry(:,1) = vbar(0, LBvj:UBvj, 1)
          else                   ! East
            vbar_bry(:,1) = vbar(L, LBvj:UBvj, 1)
          endif      
          start2D = (/ 1,    itime /)
          count2D = (/ Nvbry, 1     /)
          write(*,*)  'Write: ', 'vbar_'//trim(BRY_NAME(ibry))
          call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
          call check( nf90_inq_varid(ncid2, 'vbar_'//trim(BRY_NAME(ibry)), var_id2) )
          call check( nf90_put_var(ncid2, var_id2, vbar_bry, start = start2D, count = count2D) )
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
        bryname = trim( VAR_NAME(i) )//'_'//trim(BRY_NAME(ibry))

#if defined ROMS_MODEL   
        start4D = (/ Irdg_min+1, Jrdg_min+1, 1,      idt(itime) /)
        count4D = (/ Nxr_dg,     Nyr_dg,     Nzr_dg, 1          /)
        if( i==8  .or. i==9  .or. i==13 .or. i==14 .or. i==15 .or. i==16 .or.  &
            i==17 .or. i==21 .or. i==22 .or. i==23 .or. i==24 .or. i==26 .or.  &
            i==27 .or. i==28 .or. i==29 .or. i==30 .or. i==33 .or. i==34 .or.  &
            i==35 .or. i==36  ) then  ! mud_, sand_, etc...
          write(varnum,'(I2.2)') j
          varname = trim( VAR_NAME(i) )//varnum
          bryname = trim( VAR_NAME(i) )//trim(BRY_NAME(ibry))//'_'//varnum
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
        write(*,*) 'Read: '//HYCOM_FILE(iNC)
        start4D = (/ Irdg_min+1, Jrdg_min+1, 1,      idt(itime) /)
        count4D = (/ Nxr_dg,     Nyr_dg,     Nzr_dg, 1          /)
# if defined FAST_READ
        call readNetCDF_4d_hycom_fast( ncid, trim( HY_VAR_NAME )   &
              , Nxr_dg, Nyr_dg, Nzr_dg, 1, start4D, count4D, t_dg )
# else
        call readNetCDF_4d_hycom( ncid, trim( HY_VAR_NAME )        &
              , Nxr_dg, Nyr_dg, Nzr_dg, 1, start4D, count4D, t_dg )
# endif
#elif defined JCOPE_MODEL
        JCOPE_data_file = trim( JCOPE_data_dir )//trim( JCOPE_prefix(i) )//trim( JCOPE_sufix(itime) )
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
            , LBri, UBri, LBrj, UBrj, 1, N                    &
            , Id_cnt3Dr, w_cnt3Dr                             &
            , t  ) 
  
        if(ibry == 1) then     ! South
          t_bry(:,:,1) = t(LBri:UBri, 0, 1:N, 1)
        elseif(ibry == 2) then ! North
          t_bry(:,:,1) = t(LBri:UBri, M, 1:N, 1)
        elseif(ibry == 3) then ! West
          t_bry(:,:,1) = t(0, LBrj:UBrj, 1:N, 1)
        else                   ! East
          t_bry(:,:,1) = t(L, LBrj:UBrj, 1:N, 1)
        endif
        start3D = (/ 1,    1,   itime /)
        count3D = (/ Nrbry, Nzr, 1     /)
        write(*,*)  'Write: ', trim( bryname )
        call check( nf90_open( trim( BRY_FILE ), NF90_WRITE, ncid2) )
        call check( nf90_inq_varid(ncid2, trim( bryname ), var_id2 ) )
        call check( nf90_put_var(ncid2, var_id2, t_bry, start = start3D, count = count3D) )
        call check( nf90_close(ncid2) )      
      enddo

      if ( itime < Nt ) then
        if( iNC == iNCt(itime+1) ) cycle
      endif
  
      deallocate( zeta_dg )
      deallocate( ubar_dg )
      deallocate( vbar_dg )
      deallocate( u_dg )
      deallocate( v_dg )
      deallocate( ull_dg  )
      deallocate( ullu_dg )
      deallocate( vllu_dg )
      deallocate( vll_dg  )
      deallocate( vllv_dg )
      deallocate( ullv_dg )
      deallocate( t_dg )
#if defined WET_DRY
      deallocate( umask_wet_dg )
      deallocate( vmask_wet_dg )
#endif

#if defined HYCOM_MODEL
      deallocate( latr_dg )
      deallocate( lonr_dg )
      deallocate( latu_dg )
      deallocate( lonu_dg )
      deallocate( latv_dg )
      deallocate( lonv_dg )
      deallocate( rmask_dg )
      deallocate( umask_dg )
      deallocate( vmask_dg )
      deallocate( cosAu_dg )
      deallocate( sinAu_dg )
      deallocate( cosAv_dg )
      deallocate( sinAv_dg )
      deallocate( h_dg )
      deallocate( z_r_dg )
      deallocate( z_u_dg )
      deallocate( z_v_dg )
    
      call check( nf90_close(ncid) )
      write(*,*) "CLOSE: ", HYCOM_FILE(iNC)
#endif

    enddo

    deallocate( zeta )
    deallocate( ubar )
    deallocate( vbar )
    deallocate( u )
    deallocate( v )
    deallocate( ull )
    deallocate( uu  )
    deallocate( vu  )
    deallocate( vll )
    deallocate( uv  )
    deallocate( vv  )
    deallocate( t )
    deallocate( zeta_bry )
    deallocate( ubar_bry )
    deallocate( vbar_bry )
    deallocate( u_bry )
    deallocate( v_bry )
    deallocate( t_bry )
  
    deallocate( ID_cnt2Dr )
    deallocate( w_cnt2Dr  )
    deallocate( ID_cnt2Du )
    deallocate( w_cnt2Du  )
    deallocate( ID_cnt2Dv )
    deallocate( w_cnt2Dv  )
    deallocate( ID_cnt3Dr )
    deallocate( w_cnt3Dr  )
    deallocate( ID_cnt3Du )
    deallocate( w_cnt3Du  )
    deallocate( ID_cnt3Dv )
    deallocate( w_cnt3Dv  )

    write(*,*) '***************** '//trim(BRY_NAME(ibry))//' boundary complete! *****************'
  enddo

#if defined ROMS_MODEL
  call check( nf90_close(ncid) )
#endif

  write(*,*) 'FINISH!!'
      
END PROGRAM bryOCN2ROMS
      
