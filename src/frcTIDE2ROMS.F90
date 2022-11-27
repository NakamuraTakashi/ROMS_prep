
!!!=== Copyright (c) 2022 Takashi NAKAMURA  ===== 

#if defined NAOTIDE
# undef  UV_TIDES
#elif defined NAOTIDEJ
# define UV_TIDES
#elif defined FES2014
# define UV_TIDES
#elif defined TPXO9_ATLAS_NC
# define UV_TIDES
#endif

PROGRAM frcTIDE2ROMS
  use netcdf
  use mod_utility
  use mod_roms_netcdf
  use mod_calendar
  use mod_interpolation
  use mod_tide
 
  implicit none
      
! ---------------------------------------------------------------------

  character(256) :: GRID_FILE

  character(256) :: TIDE_DATA_dir
  character(256) :: TIDE_prefix
  character(15) :: TIDE_suffix   = "_20000101.00.nc"
  character(256) :: TIDE_FILE

  character(256) :: TDATA_FILE  !NAO99b / NAO99Jb data file name
  
  integer :: start1D(1), count1D(1)
  integer :: start2D(2), count2D(2)
  integer :: start3D(3), count3D(3)
  integer :: start4D(4), count4D(4)
  
  real(8), allocatable :: h(:,:)        ! depth (meter)
  real(8), allocatable :: rmask(:,:)    ! land mask
  real(8), allocatable :: latr(:,:)
  real(8), allocatable :: lonr(:,:)
  real(8), allocatable :: SSH_Tamp(:,:,:)
  real(8), allocatable :: SSH_Tphase(:,:,:)
  real(8), allocatable :: SSH_Tcos(:,:,:)
  real(8), allocatable :: SSH_Tsin(:,:,:)

  real(8), allocatable :: latr_dg(:,:)
  real(8), allocatable :: lonr_dg(:,:)
  real(8), allocatable :: rmask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: SSH_Tamp_dg(:,:,:)
  real(8), allocatable :: SSH_Tphase_dg(:,:,:)
  real(8), allocatable :: SSH_Tcos_dg(:,:,:)
  real(8), allocatable :: SSH_Tsin_dg(:,:,:)
  real(8), allocatable :: amp(:,:), phs(:,:)

  integer, allocatable :: ID_cnt2Dr(:,:)
  real(8), allocatable :: w_cnt2Dr(:,:)

#if defined UV_TIDES
  real(8), allocatable :: UV_Tangle(:,:,:)
  real(8), allocatable :: UV_Tmajor(:,:,:)
  real(8), allocatable :: UV_Tminor(:,:,:)
  real(8), allocatable :: UV_Tphase(:,:,:)

  real(8), allocatable :: latu_dg(:,:)
  real(8), allocatable :: lonu_dg(:,:)
  real(8), allocatable :: latv_dg(:,:)
  real(8), allocatable :: lonv_dg(:,:)
  real(8), allocatable :: umask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: vmask_dg(:,:)  ! land mask of donor grid

  integer, allocatable :: ID_cnt2Du(:,:)
  real(8), allocatable :: w_cnt2Du(:,:)
  integer, allocatable :: ID_cnt2Dv(:,:)
  real(8), allocatable :: w_cnt2Dv(:,:)

  real(8), allocatable :: Ua_dg(:,:), Up_dg(:,:), Va_dg(:,:), Vp_dg(:,:)

  real(8), allocatable :: Uamp  (:,:)
  real(8), allocatable :: Uphase(:,:)
  real(8), allocatable :: Upcos(:,:)
  real(8), allocatable :: Upsin(:,:)
  real(8), allocatable :: Vamp  (:,:)
  real(8), allocatable :: Vphase(:,:)
  real(8), allocatable :: Vpcos(:,:)
  real(8), allocatable :: Vpsin(:,:)

  real(8), allocatable :: Uamp_dg  (:,:,:)
  real(8), allocatable :: Uphase_dg(:,:,:)
  real(8), allocatable :: Upcos_dg(:,:,:)
  real(8), allocatable :: Upsin_dg(:,:,:)
  real(8), allocatable :: Vamp_dg  (:,:,:)
  real(8), allocatable :: Vphase_dg(:,:,:)
  real(8), allocatable :: Vpcos_dg(:,:,:)
  real(8), allocatable :: Vpsin_dg(:,:,:)
#endif

  integer :: ncid,var_id
  integer :: Nxr, Nyr !, Nxu, Nyu, Nxv, Nyv
  integer :: L, M  
  integer :: Nxr_dg, Nyr_dg, Nxu_dg, Nyu_dg, Nxv_dg, Nyv_dg
  integer :: Irdg_min, Irdg_max, Jrdg_min, Jrdg_max  ! Minimum and maximum ID of donor grid index.
  integer :: Iudg_min, Iudg_max, Judg_min, Judg_max  ! Minimum and maximum ID of donor grid index.
  integer :: Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max  ! Minimum and maximum ID of donor grid index.
  integer :: Ldg, Mdg, Ndg
  integer :: dimids(3)
  character(256) :: varname

  integer :: irg, jrg, krg
  integer :: Idg, Jdg, Kdg
  integer :: i,j,iTC
  real(8) :: cff

#if defined NAOTIDE || defined NAOTIDEJ

# if defined NAOTIDE
  character(6) :: SOURCE_NAME = "NAO99b"
# elif defined NAOTIDEJ
  character(7) :: SOURCE_NAME = "NAO99Jb"
# endif

  integer, parameter :: NTC = 16
  character(3), parameter :: Tconsti(NTC) = (/    &
      'k2  '          &
    , 's2  '          &
    , 't2  '          &
    , 'l2  '          &
    , 'm2  '          &
    , 'nu2 '          &
    , 'n2  '          &
    , 'mu2 '          &
    , '2n2 '          &
    , 'oo1 '          &
    , 'j1  '          &
    , 'k1  '          &
    , 'p1  '          &
    , 'm1  '          &
    , 'o1  '          &
    , 'q1  '          &
    /)
  real(8), parameter :: Tperiod(NTC) = (/     &           
      11.96723479d0  &  ! k2 
    , 12.00000000d0  &  ! s2 
    , 12.01644920d0  &  ! t2 
    , 12.19162016d0  &  ! l2 
    , 12.42060122d0  &  ! m2 
    , 12.62600437d0  &  ! nu2
    , 12.65834824d0  &  ! n2 
    , 12.87175759d0  &  ! mu2
    , 12.90537448d0  &  ! 2n2
    , 22.30607420d0  &  ! oo1
    , 23.09847677d0  &  ! j1 
    , 23.93446966d0  &  ! k1 
    , 24.06589016d0  &  ! p1 
    , 24.83324836d0  &  ! m1 
    , 25.81934166d0  &  ! o1 
    , 26.86835667d0  &  ! q1 
    /)
  real(8), parameter :: Tphsini(NTC) = (/     &           
      199.9355894d0  &  ! k2 
    , 359.9902258d0  &  ! s2 
    ,   2.9548684d0  &  ! t2 
    ,  84.8961568d0  &  ! l2 
    , 136.4519732d0  &  ! m2 
    ,  41.3579043d0  &  ! nu2
    ,   8.0077897d0  &  ! n2 
    , 272.9137207d0  &  ! mu2
    , 239.5636061d0  &  ! 2n2
    ,  73.4514109d0  &  ! oo1
    , 138.4119783d0  &  ! j1 
    ,   9.9677947d0  &  ! k1 
    , 350.0224311d0  &  ! p1 
    , 241.5236112d0  &  ! m1 
    , 126.4841785d0  &  ! o1 
    , 358.0399950d0  &  ! q1 
    /)

  real(8) :: xmin, xmax    ! 
  real(8) :: ymin, ymax    ! 
  real(8) :: dx, dy        ! 
  integer :: mmax, nmax    ! 
  integer :: ideff         ! 
  character(6) :: fmt      !
  real(8) :: aunit         ! 
  real(8) :: punit         ! 
 
#elif defined FES2014
  character(7) :: SOURCE_NAME = "FES2014"
  integer, parameter :: NTC = 34
  character(4), parameter :: Tconsti(NTC) = (/    &
      '2n2 '          &
    , 'k2  '          &
    , 'm3  '          &
    , 'mf  '          &
    , 'ms4 '          &
    , 'mu2 '          &
    , 'o1  '          &
    , 's1  '          &
    , 'ssa '          &
    , 'eps2'          &
    , 'l2  '          &
    , 'm4  '          &
    , 'mks2'          &
    , 'msf '          &
    , 'n2  '          &
    , 'p1  '          &
    , 's2  '          &
    , 't2  '          &
    , 'j1  '          &
    , 'la2 '          &
    , 'm6  '          &
    , 'mm  '          &
    , 'msqm'          &
    , 'n4  '          &
    , 'q1  '          &
    , 's4  '          &
    , 'k1  '          &
    , 'm2  '          &
    , 'm8  '          &
    , 'mn4 '          &
    , 'mtm '          &
    , 'nu2 '          &
    , 'r2  '          &
    , 'sa  '          &
    /)
  real(8), parameter :: Tperiod(NTC) = (/     &           
      11.96723479d0  &  ! 2n2 
    , 12.00000000d0  &  ! k2  
    , 12.01644920d0  &  ! m3  
    , 12.19162016d0  &  ! mf  
    , 12.42060122d0  &  ! ms4 
    , 12.62600437d0  &  ! mu2 
    , 12.65834824d0  &  ! o1  
    , 12.87175759d0  &  ! s1  
    , 12.90537448d0  &  ! ssa 
    , 22.30607420d0  &  ! eps2
    , 23.09847677d0  &  ! l2  
    , 23.93446966d0  &  ! m4  
    , 24.06589016d0  &  ! mks2
    , 24.83324836d0  &  ! msf 
    , 25.81934166d0  &  ! n2  
    , 26.86835667d0  &  ! p1  
    , 12.00000000d0  &  ! s2  
    , 12.01644920d0  &  ! t2  
    , 12.19162016d0  &  ! j1  
    , 12.42060122d0  &  ! la2 
    , 12.62600437d0  &  ! m6  
    , 12.65834824d0  &  ! mm  
    , 12.87175759d0  &  ! msqm
    , 12.90537448d0  &  ! n4  
    , 22.30607420d0  &  ! q1  
    , 23.09847677d0  &  ! s4  
    , 23.93446966d0  &  ! k1  
    , 24.06589016d0  &  ! m2  
    , 24.83324836d0  &  ! m8  
    , 25.81934166d0  &  ! mn4 
    , 26.86835667d0  &  ! mtm 
    , 12.00000000d0  &  ! nu2 
    , 12.01644920d0  &  ! r2  
    , 12.19162016d0  &  ! sa  
    /)
  real(8), parameter :: Tphsini(NTC) = (/     &           
      199.9355894d0  &  ! 2n2 
    , 359.9902258d0  &  ! k2  
    ,   2.9548684d0  &  ! m3  
    ,  84.8961568d0  &  ! mf  
    , 136.4519732d0  &  ! ms4 
    ,  41.3579043d0  &  ! mu2 
    ,   8.0077897d0  &  ! o1  
    , 272.9137207d0  &  ! s1  
    , 239.5636061d0  &  ! ssa 
    ,  73.4514109d0  &  ! eps2
    , 138.4119783d0  &  ! l2  
    ,   9.9677947d0  &  ! m4  
    , 350.0224311d0  &  ! mks2
    , 241.5236112d0  &  ! msf 
    , 126.4841785d0  &  ! n2  
    , 358.0399950d0  &  ! p1  
    , 359.9902258d0  &  ! s2  
    ,   2.9548684d0  &  ! t2  
    ,  84.8961568d0  &  ! j1  
    , 136.4519732d0  &  ! la2 
    ,  41.3579043d0  &  ! m6  
    ,   8.0077897d0  &  ! mm  
    , 272.9137207d0  &  ! msqm
    , 239.5636061d0  &  ! n4  
    ,  73.4514109d0  &  ! q1  
    , 138.4119783d0  &  ! s4  
    ,   9.9677947d0  &  ! k1  
    , 136.4519732d0  &  ! m2  
    , 241.5236112d0  &  ! m8  
    , 126.4841785d0  &  ! mn4 
    , 358.0399950d0  &  ! mtm 
    , 359.9902258d0  &  ! nu2 
    ,   2.9548684d0  &  ! r2  
    ,  84.8961568d0  &  ! sa  
    /)

  real(8), allocatable :: lat_dg(:), lon_dg(:)
    
#endif

!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/tide/TIDE_DATA_dir
  namelist/tide/TIDE_prefix

  read (5, nml=grd)
  rewind(5)
  read (5, nml=tide)

!---- Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
  L = Nxr-1
  M = Nyr-1
  allocate( latr(0:L, 0:M) )
  allocate( lonr(0:L, 0:M) )

  allocate ( SSH_Tamp(0:L, 0:M, NTC) )
  allocate ( SSH_Tphase(0:L, 0:M, NTC) )
  allocate ( SSH_Tcos(0:L, 0:M, NTC) )
  allocate ( SSH_Tsin(0:L, 0:M, NTC) )
#if defined UV_TIDES
  allocate ( UV_Tangle(0:L, 0:M, NTC) )
  allocate ( UV_Tmajor(0:L, 0:M, NTC) )
  allocate ( UV_Tminor(0:L, 0:M, NTC) )
  allocate ( UV_Tphase(0:L, 0:M, NTC) )
  allocate ( Uamp(0:L, 0:M) )
  allocate ( Uphase(0:L, 0:M) )
  allocate ( Upcos(0:L, 0:M) )
  allocate ( Upsin(0:L, 0:M) )
  allocate ( Vamp(0:L, 0:M) )
  allocate ( Vphase(0:L, 0:M) )
  allocate ( Vpcos(0:L, 0:M) )
  allocate ( Vpsin(0:L, 0:M) )
#endif


  allocate(ID_cnt2Dr(4, Nxr*Nyr))
  allocate(w_cnt2Dr(4, Nxr*Nyr))
#if defined UV_TIDES
  allocate(ID_cnt2Du(4, Nxr*Nyr))
  allocate(w_cnt2Du(4, Nxr*Nyr))
  allocate(ID_cnt2Dv(4, Nxr*Nyr))
  allocate(w_cnt2Dv(4, Nxr*Nyr))
#endif

  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_north)
  call check( nf90_get_var(ncid, var_id, latr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, lonr) )
  
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
!---- Read NAO99/NAO99J grid --------------------------------
#if defined NAOTIDE || defined NAOTIDEJ
# if defined NAOTIDEJ
  TDATA_FILE = trim(TIDE_DATA_dir)//'omapj/'//trim(Tconsti(1))//'_j.nao'
# else
  TDATA_FILE = trim(TIDE_DATA_dir)//'omap/'//trim(Tconsti(1))//'.nao'
# endif
  write(*,*) "OPEN: "//trim(TDATA_FILE)

  open(50, file=trim(TDATA_FILE), status='old')
  call readNAO99_head( 50, xmin, xmax, ymin, ymax, dx, dy         &
                     , Nxr_dg, Nyr_dg, ideff, fmt, aunit, punit )

  Nxu_dg = Nxr_dg
  Nyu_dg = Nyr_dg
  Nxv_dg = Nxr_dg
  Nyv_dg = Nyr_dg
  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1

  allocate( amp(0:Ldg, 0:Mdg) )
  allocate( phs(0:Ldg, 0:Mdg) )

  call readNAO99_data( 50, Nxr_dg, Nyr_dg, fmt, aunit, punit  &
                           , amp, phs  )
  close(50)
!---------------------------------------

  allocate( latr_dg(0:Ldg, 0:Mdg) )
  allocate( lonr_dg(0:Ldg, 0:Mdg) )
  allocate( rmask_dg(0:Ldg, 0:Mdg) )

  do j=0,Mdg
    do i=0,Ldg
      latr_dg(i,j) = ymin + dy*dble(j)
      lonr_dg(i,j) = xmin + dx*dble(i)
    enddo
  enddo
  rmask_dg(:,:) = 1.0d0
  do j=0,Mdg
    do i=0,Ldg
      if (amp(i, j)>=99.0d0) rmask_dg(i,j) = 0.0d0 
    enddo
  enddo

# if defined UV_TIDES

  allocate( Ua_dg(0:Ldg, 0:Mdg) )
  allocate( Up_dg(0:Ldg, 0:Mdg) )
  allocate( Va_dg(0:Ldg, 0:Mdg) )
  allocate( Vp_dg(0:Ldg, 0:Mdg) )
  allocate( latu_dg(0:Ldg, 0:Mdg) )
  allocate( lonu_dg(0:Ldg, 0:Mdg) )
  allocate( latv_dg(0:Ldg, 0:Mdg) )
  allocate( lonv_dg(0:Ldg, 0:Mdg) )
  allocate( umask_dg(0:Ldg, 0:Mdg) )
  allocate( vmask_dg(0:Ldg, 0:Mdg) )

!---- Read NAO99J velocity grid --------------------------------
!! ----- Read velocity from nao99Jb original data ----------------------------
!  TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/vfield.'//trim(Tconsti(1))
!  write(*,*) "OPEN: "//trim(TDATA_FILE)
!
!  open(50, file=trim(TDATA_FILE), status='old')
!  call readNAO99_vel( 50, Nxr_dg, Nyr_dg, latr_dg, lonr_dg    &
!                    , Ua_dg, Up_dg, Va_dg, Vp_dg  )
!  close(50)

! ----- Read velocity of gridded nao99Jb data --------------------------------
  TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Au_'//trim(Tconsti(1))//'.dat'
  write(*,*) "OPEN: "//trim(TDATA_FILE)
  open(50, file=trim(TDATA_FILE), status='old')
  do j=0,Mdg
    read(50,*) (Ua_dg(i,j), i=0,Ldg)
  enddo
  close(50)

  TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Av_'//trim(Tconsti(1))//'.dat'
  write(*,*) "OPEN: "//trim(TDATA_FILE)
  open(50, file=trim(TDATA_FILE), status='old')
  do j=0,Mdg
    read(50,*) (Va_dg(i,j), i=0,Ldg)
  enddo
  close(50)
! ----------------------------------------------------------------------------- 

  do j=0,Mdg
    do i=0,Ldg
      latu_dg(i,j) = latr_dg(i,j)
      lonu_dg(i,j) = lonr_dg(i,j) - 0.5d0*dx
      latv_dg(i,j) = latr_dg(i,j) - 0.5d0*dy
      lonv_dg(i,j) = lonr_dg(i,j)
    enddo
  enddo

  umask_dg(:,:) = 1.0d0
  vmask_dg(:,:) = 1.0d0
  do j=0,Mdg
    do i=0,Ldg
      if (Ua_dg(i, j)>=9999.0d0) umask_dg(i,j) = 0.0d0
      if (Va_dg(i, j)>=9999.0d0) vmask_dg(i,j) = 0.0d0
    enddo
  enddo

# endif

!---- FES2014 data grid preparation --------------------------------
#elif defined FES2014

  TDATA_FILE = trim(TIDE_DATA_dir)//'ocean_tide/'//trim(Tconsti(1))//'.nc'
  write(*,*) "OPEN: "//trim(TDATA_FILE)

  ! Open NetCDF FES2014 onean tide file
  call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'lon', Nxr_dg)
  call get_dimension(ncid, 'lat', Nyr_dg)

  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1

  allocate( lat_dg(0:Mdg) )
  allocate( lon_dg(0:Ldg) )
  allocate( amp(0:Ldg, 0:Mdg) )
  allocate( phs(0:Ldg, 0:Mdg) )
  allocate( latr_dg(0:Ldg, 0:Mdg) )
  allocate( lonr_dg(0:Ldg, 0:Mdg) )
  allocate( rmask_dg(0:Ldg, 0:Mdg) )

  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat', var_id) ) ! latitude (degrees_north)
  call check( nf90_get_var(ncid, var_id, lat_dg) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) ) ! longitude (degrees_east)
  call check( nf90_get_var(ncid, var_id, lon_dg) )
  call check( nf90_inq_varid(ncid, 'amplitude', var_id) ) ! Sea surface height amplitude (cm)
  call check( nf90_get_var(ncid, var_id, amp) )
  call check( nf90_inq_varid(ncid, 'phase', var_id) ) ! Sea surface height phaselag (degree)
  call check( nf90_get_var(ncid, var_id, phs) )
  
  ! Close NetCDF file
  call check( nf90_close(ncid) )
!---------------------------------------

  do i=0,Ldg
      lonr_dg(i,:) = lon_dg(i)
  enddo
  do j=0,Mdg
      latr_dg(:,j) = lat_dg(j)
  enddo
  rmask_dg(:,:) = 1.0d0
  do j=0,Mdg
    do i=0,Ldg
      if (amp(i, j)>=99999.0d0) rmask_dg(i,j) = 0.0d0 
    enddo
  enddo
  deallocate( lat_dg, lon_dg )

# if defined UV_TIDES

! ----- Read velocity of gridded FES2014 data --------------------------------
  TDATA_FILE = trim(TIDE_DATA_dir)//'eastward_velocity/'//trim(Tconsti(1))//'.nc'
  write(*,*) "OPEN: "//trim(TDATA_FILE)
  open(50, file=trim(TDATA_FILE), status='old')
  ! Open NetCDF FES2014 onean tide file
  call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'lon', Nxu_dg)
  call get_dimension(ncid, 'lat', Nyu_dg)

  Ldg = Nxu_dg-1
  Mdg = Nyu_dg-1

  allocate( lat_dg(0:Mdg) )
  allocate( lon_dg(0:Ldg) )
  allocate( Ua_dg(0:Ldg, 0:Mdg) )
  allocate( Up_dg(0:Ldg, 0:Mdg) )
  allocate( latu_dg(0:Ldg, 0:Mdg) )
  allocate( lonu_dg(0:Ldg, 0:Mdg) )
  allocate( umask_dg(0:Ldg, 0:Mdg) )
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat', var_id) ) ! latitude (degrees_north)
  call check( nf90_get_var(ncid, var_id, lat_dg) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) ) ! longitude (degrees_east)
  call check( nf90_get_var(ncid, var_id, lon_dg) )
  call check( nf90_inq_varid(ncid, 'Ua', var_id) ) ! Eastward sea water velocity amplitude (cm/s)
  call check( nf90_get_var(ncid, var_id, Ua_dg) )
  call check( nf90_inq_varid(ncid, 'Ug', var_id) ) ! Eastward sea water velocity phaselag (degree)
  call check( nf90_get_var(ncid, var_id, Up_dg) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )
!---------------------------------------

  do i=0,Ldg
      lonu_dg(i,:) = lon_dg(i)
  enddo
  do j=0,Mdg
      latu_dg(:,j) = lat_dg(j)
  enddo
  umask_dg(:,:) = 1.0d0
  do j=0,Mdg
    do i=0,Ldg
      if (Ua_dg(i, j)>=99999.0d0) umask_dg(i,j) = 0.0d0 
    enddo
  enddo
  deallocate( lat_dg, lon_dg )

! ----------------------------------------------------------------------------- 
 
  TDATA_FILE = trim(TIDE_DATA_dir)//'northward_velocity/'//trim(Tconsti(1))//'.nc'
  write(*,*) "OPEN: "//trim(TDATA_FILE)
  ! Open NetCDF FES2014 onean tide file
  call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'lon', Nxr_dg)
  call get_dimension(ncid, 'lat', Nyr_dg)

  allocate( lat_dg(0:Mdg) )
  allocate( lon_dg(0:Ldg) )
  allocate( Va_dg(0:Ldg, 0:Mdg) )
  allocate( Vp_dg(0:Ldg, 0:Mdg) )
  allocate( latv_dg(0:Ldg, 0:Mdg) )
  allocate( lonv_dg(0:Ldg, 0:Mdg) )
  allocate( vmask_dg(0:Ldg, 0:Mdg) )
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'lat', var_id) ) ! latitude (degrees_north)
  call check( nf90_get_var(ncid, var_id, lat_dg) )
  call check( nf90_inq_varid(ncid, 'lon', var_id) ) ! longitude (degrees_east)
  call check( nf90_get_var(ncid, var_id, lon_dg) )
  call check( nf90_inq_varid(ncid, 'Va', var_id) ) ! Northward sea water velocity amplitude (cm/s)
  call check( nf90_get_var(ncid, var_id, Va_dg) )
  call check( nf90_inq_varid(ncid, 'Vg', var_id) ) ! Northward sea water velocity phaselag (degree)
  call check( nf90_get_var(ncid, var_id, Vp_dg) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )

  do i=0,Ldg
      lonv_dg(i,:) = lon_dg(i)
  enddo
  do j=0,Mdg
      latv_dg(:,j) = lat_dg(j)
  enddo
  vmask_dg(:,:) = 1.0d0
  do j=0,Mdg
    do i=0,Ldg
      if (Va_dg(i, j)>=99999.0d0) vmask_dg(i,j) = 0.0d0 
    enddo
  enddo
  deallocate( lat_dg, lon_dg )

# endif
#endif

!---- Seek data range --------------------------------

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

  Nxr_dg =Irdg_max-Irdg_min+1
  Nyr_dg =Jrdg_max-Jrdg_min+1
  Nxu_dg =Iudg_max-Iudg_min+1
  Nyu_dg =Judg_max-Judg_min+1
  Nxv_dg =Ivdg_max-Ivdg_min+1
  Nyv_dg =Jvdg_max-Jvdg_min+1

  write(*,*) Irdg_min, Irdg_max, Jrdg_min, Jrdg_max

!  ---- Calc. weight parameters for interpolation -------------- 

  write(*,*) "Calculate 2D rho point weight factor"
  call weight2D_grid3_2(                                &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , 0, L,   0, M,   lonr,    latr                 &
        , ID_cnt2Dr, w_cnt2Dr )

#if defined UV_TIDES
  write(*,*) "Calculate 2D u point weight factor"
  call weight2D_grid3_2(                                &
          Iudg_min, Iudg_max, Judg_min, Judg_max        &
        , lonu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
        , latu_dg(Iudg_min:Iudg_max,Judg_min:Judg_max)  &
        , umask_dg(Iudg_min:Iudg_max,Judg_min:Judg_max) & 
        , 0, L,   0, M,   lonr,    latr                 &
        , ID_cnt2Du, w_cnt2Du )

  write(*,*) "Calculate 2D v point weight factor"
  call weight2D_grid3_2(                                &
          Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max        &
        , lonv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
        , latv_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max)  &
        , vmask_dg(Ivdg_min:Ivdg_max,Jvdg_min:Jvdg_max) & 
        , 0, L,   0, M,   lonr,    latr                 &
        , ID_cnt2Dv, w_cnt2Dv )
#endif

  allocate( SSH_Tcos_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, NTC) )
  allocate( SSH_Tsin_dg(Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, NTC) )
  allocate( SSH_Tamp_dg  (Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, NTC) )


!=== Read Tide data ================================================

  DO iTC=1,NTC
!---- Read NAO99/NAO99Jb data --------------------------------
#if defined NAOTIDE || defined NAOTIDEJ
# if defined NAOTIDEJ
    TDATA_FILE = trim(TIDE_DATA_dir)//'omapj/'//trim(Tconsti(iTC))//'_j.nao'
# else
    TDATA_FILE = trim(TIDE_DATA_dir)//'omap/'//trim(Tconsti(iTC))//'.nao'
# endif
    write(*,*) "OPEN: "//trim(TDATA_FILE)
  
    open(50, file=trim(TDATA_FILE), status='old')
    call readNAO99_head( 50, xmin, xmax, ymin, ymax, dx, dy         &
                       , Nxr_dg, Nyr_dg, ideff, fmt, aunit, punit )
    call readNAO99_data( 50, Nxr_dg, Nyr_dg, fmt, aunit, punit  &
                             , amp, phs  )
    close(50)

!---- Read FES2014 data --------------------------------
#elif defined FES2014
    TDATA_FILE = trim(TIDE_DATA_dir)//'ocean_tide/'//trim(Tconsti(iTC))//'.nc'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
  
    ! Open NetCDF FES2014 onean tide file
    call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
    
!    start2D = (/ Irdg_min+1, Jrdg_min+1 /)
!    count2D = (/ Nxr_dg,     Nyr_dg     /)
    call check( nf90_inq_varid(ncid, 'amplitude', var_id) ) ! Sea surface height amplitude (cm)
    call check( nf90_get_var(ncid, var_id, amp ) )
    call check( nf90_inq_varid(ncid, 'phase', var_id) ) ! Sea surface height phaselag (degree)
    call check( nf90_get_var(ncid, var_id, phs ) )   
    ! Close NetCDF file
    call check( nf90_close(ncid) )

    amp(:,:)=amp(:,:)*0.01d0  ! cm -> m
  
#endif
  
    do j=Jrdg_min,Jrdg_max
      do i=Irdg_min,Irdg_max
        SSH_Tamp_dg(i,j,iTC)   = amp(i,j)
        cff = (phs(i,j)-Tphsini(iTC))*deg2rad
        SSH_Tcos_dg(i,j,iTC) = cos(cff)
        SSH_Tsin_dg(i,j,iTC) = sin(cff)
      enddo
    enddo

  END DO

!---- Read velocity data --------------------------------
#if defined UV_TIDES
  allocate( Uamp_dg  (Irdg_min:Irdg_max, Jrdg_min:Jrdg_max, NTC) )
  allocate( Uphase_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )
  allocate( Upcos_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )
  allocate( Upsin_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )
  allocate( Vamp_dg  (Ivdg_min:Ivdg_max, Jvdg_min:Jvdg_max, NTC) )
  allocate( Vphase_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )
  allocate( Vpcos_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )
  allocate( Vpsin_dg(Iudg_min:Iudg_max, Judg_min:Judg_max, NTC) )


  DO iTC=1,NTC
!---- Read NAO99/NAO99J data --------------------------------
# if defined NAOTIDEJ
!! ----- Read velocity from nao99Jb original data ----------------------------
!    
!    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/vfield.'//trim(Tconsti(iTC))
!    write(*,*) "OPEN: "//trim(TDATA_FILE)
!  
!    open(50, file=trim(TDATA_FILE), status='old')
!    call readNAO99_vel( 50, Nxr_dg, Nyr_dg, latr_dg, lonr_dg    &
!                      , Ua_dg, Up_dg, Va_dg, Vp_dg  )
!    close(50)
!
!! ----- Write velocity of gridded nao99Jb data ------------------------------
!
!    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Au_'//trim(Tconsti(iTC))//'.dat'
!    open(50, file=trim(TDATA_FILE), status='replace')
!    do j=0,Mdg
!      write(50,*) (Ua_dg(i,j), i=0,Ldg)
!    enddo
!    close(50)
!
!    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Av_'//trim(Tconsti(iTC))//'.dat'
!    open(50, file=trim(TDATA_FILE), status='replace')
!    do j=0,Mdg
!      write(50,*) (Va_dg(i,j), i=0,Ldg)
!    enddo
!    close(50)
!
!    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Pu_'//trim(Tconsti(iTC))//'.dat'
!    open(50, file=trim(TDATA_FILE), status='replace')
!    do j=0,Mdg
!      write(50,*) (Up_dg(i,j), i=0,Ldg)
!    enddo
!    close(50)
!
!    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Pv_'//trim(Tconsti(iTC))//'.dat'
!    open(50, file=trim(TDATA_FILE), status='replace')
!    do j=0,Mdg
!      write(50,*) (Vp_dg(i,j), i=0,Ldg)
!    enddo
!    close(50)
!

! ----- Read velocity of gridded nao99Jb data --------------------------------

    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Au_'//trim(Tconsti(iTC))//'.dat'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    open(50, file=trim(TDATA_FILE), status='old')
    do j=0,Mdg
      read(50,*) (Ua_dg(i,j), i=0,Ldg)
    enddo
    close(50)

    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Av_'//trim(Tconsti(iTC))//'.dat'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    open(50, file=trim(TDATA_FILE), status='old')
    do j=0,Mdg
      read(50,*) (Va_dg(i,j), i=0,Ldg)
    enddo
    close(50)

    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Pu_'//trim(Tconsti(iTC))//'.dat'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    open(50, file=trim(TDATA_FILE), status='old')
    do j=0,Mdg
      read(50,*) (Up_dg(i,j), i=0,Ldg)
    enddo
    close(50)

    TDATA_FILE = trim(TIDE_DATA_dir)//'nao99Jb_vel/Pv_'//trim(Tconsti(iTC))//'.dat'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    open(50, file=trim(TDATA_FILE), status='old')
    do j=0,Mdg
      read(50,*) (Vp_dg(i,j), i=0,Ldg)
    enddo
    close(50)

!---- Read FES2014 data --------------------------------
# elif defined FES2014
! ----- Read velocity of gridded FES2014 data --------------------------------
    TDATA_FILE = trim(TIDE_DATA_dir)//'eastward_velocity/'//trim(Tconsti(iTC))//'.nc'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    open(50, file=trim(TDATA_FILE), status='old')
    ! Open NetCDF FES2014 onean tide file
    call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
    ! Read data
    call check( nf90_inq_varid(ncid, 'Ua', var_id) ) ! Eastward sea water velocity amplitude (cm/s)
    call check( nf90_get_var(ncid, var_id, Ua_dg) )
    call check( nf90_inq_varid(ncid, 'Ug', var_id) ) ! Eastward sea water velocity phaselag (degree)
    call check( nf90_get_var(ncid, var_id, Up_dg) )
  ! Close NetCDF file
    call check( nf90_close(ncid) )

! ----- Read velocity of gridded FES2014 data --------------------------------
    TDATA_FILE = trim(TIDE_DATA_dir)//'northward_velocity/'//trim(Tconsti(iTC))//'.nc'
    write(*,*) "OPEN: "//trim(TDATA_FILE)
    ! Open NetCDF FES2014 onean tide file
    call check( nf90_open(trim( TDATA_FILE ), nf90_nowrite, ncid) )
    ! Read data
    call check( nf90_inq_varid(ncid, 'Va', var_id) ) ! Northward sea water velocity amplitude (cm/s)
    call check( nf90_get_var(ncid, var_id, Va_dg) )
    call check( nf90_inq_varid(ncid, 'Vg', var_id) ) ! Northward sea water velocity phaselag (degree)
    call check( nf90_get_var(ncid, var_id, Vp_dg) )
    ! Close NetCDF file
    call check( nf90_close(ncid) )

# endif
! ----------------------------------------------------------------------------- 
    
    do j=Judg_min,Judg_max
      do i=Iudg_min,Iudg_max
        Uamp_dg(i,j,iTC)  = Ua_dg(i,j)
!        cff = (Up_dg(i,j)-Tphsini(iTC))*deg2rad
        cff = Up_dg(i,j)*deg2rad
        Upcos_dg(i,j,iTC) = cos(cff)
        Upsin_dg(i,j,iTC) = sin(cff)
      enddo
    enddo
    do j=Jvdg_min,Jvdg_max
      do i=Ivdg_min,Ivdg_max
        Vamp_dg(i,j,iTC)  = Va_dg(i,j)
!        cff = (Vp_dg(i,j)-Tphsini(iTC))*deg2rad
        cff = Vp_dg(i,j)*deg2rad
        Vpcos_dg(i,j,iTC) = cos(cff)
        Vpsin_dg(i,j,iTC) = sin(cff)
      enddo
    enddo

  END DO

#endif



!---- Create the ROMS tidal forcing netCDF file --------------------------------
  
  TIDE_FILE = trim( TIDE_prefix )//TIDE_suffix

  call createNetCDFtide(  trim( TIDE_FILE ), trim( SOURCE_NAME ), Nxr, Nyr, NTC ) 

!==== LOOP START ====================================================================

  DO iTC=1,NTC

    write(*,*) 'Linear Interporation: SSH_Tamp'
    call interp2D_grid3_2(                                                 &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, SSH_Tamp_dg(:,:,iTC)   &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dr, w_cnt2Dr                                            &
          , SSH_Tamp(:,:,iTC)  ) 
    write(*,*) 'Linear Interporation: SSH_Tphase'
    call interp2D_grid3_2(                                                 &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, SSH_Tcos_dg(:,:,iTC) &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dr, w_cnt2Dr                                            &
          , SSH_Tcos(:,:,iTC)  )
    call interp2D_grid3_2(                                                 &
            Irdg_min, Irdg_max, Jrdg_min, Jrdg_max, SSH_Tsin_dg(:,:,iTC) &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dr, w_cnt2Dr                                            &
          , SSH_Tsin(:,:,iTC)  )
    
    do j=0,M
      do i=0,L
        cff = atan2( SSH_Tsin(i,j,iTC), SSH_Tcos(i,j,iTC) ) !!!!!!!!!!!!!!!!!
        SSH_Tphase(i,j,iTC) = mod(cff/deg2rad+360.0d0,360.0d0)
      enddo
    enddo

#if defined UV_TIDES
    write(*,*) 'Linear Interporation: Uamp'
    call interp2D_grid3_2(                                                 &
            Iudg_min, Iudg_max, Judg_min, Judg_max, Uamp_dg(:,:,iTC)       &
          , 0, L, 0, M                                                     &
          , Id_cnt2Du, w_cnt2Du                                            &
          , Uamp  ) 
    write(*,*) 'Linear Interporation: Upcos'
    call interp2D_grid3_2(                                                 &
            Iudg_min, Iudg_max, Judg_min, Judg_max, Upcos_dg(:,:,iTC)     &
          , 0, L, 0, M                                                     &
          , Id_cnt2Du, w_cnt2Du                                            &
          , Upcos  )
    write(*,*) 'Linear Interporation: Upsin'
    call interp2D_grid3_2(                                                 &
            Iudg_min, Iudg_max, Judg_min, Judg_max, Upsin_dg(:,:,iTC)     &
          , 0, L, 0, M                                                     &
          , Id_cnt2Du, w_cnt2Du                                            &
          , Upsin  )
    write(*,*) 'Linear Interporation: Vamp'
    call interp2D_grid3_2(                                                 &
            Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, Vamp_dg(:,:,iTC)       &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dv, w_cnt2Dv                                            &
          , Vamp )
    write(*,*) 'Linear Interporation: Vpcos'
    call interp2D_grid3_2(                                                 &
            Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, Vpcos_dg(:,:,iTC)     &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dv, w_cnt2Dv                                            &
          , Vpcos  )
    write(*,*) 'Linear Interporation: Vpsin'
    call interp2D_grid3_2(                                                 &
            Ivdg_min, Ivdg_max, Jvdg_min, Jvdg_max, Vpsin_dg(:,:,iTC)      &
          , 0, L, 0, M                                                     &
          , Id_cnt2Dv, w_cnt2Dv                                            &
          , Vpsin  )
    
    do j=0,M
      do i=0,L
        cff = atan2( Upsin(i,j), Upcos(i,j) )
        Uphase(i,j) = mod(cff/deg2rad+360.0d0,360.0d0)
        cff = atan2( Vpsin(i,j), Vpcos(i,j) )
        Vphase(i,j) = mod(cff/deg2rad+360.0d0,360.0d0)

        CALL Tvel_nao2roms( Uamp(i,j), Uphase(i,j), Vamp(i,j), Vphase(i,j) &
                          , Tphsini(iTC)                                   &
                          , UV_Tphase(i,j,iTC), UV_Tangle(i,j,iTC)         &
                          , UV_Tmajor(i,j,iTC), UV_Tminor(i,j,iTC)       )

!        UV_Tphase(i,j,iTC) = Uamp(i,j) !!!!!!!!!!!!!!!!!
!        UV_Tangle(i,j,iTC) = Uphase(i,j) !!!!!!!!!!!!!!!!!!
!        UV_Tminor(i,j,iTC) = Vamp(i,j) !!!!!!!!!!!!!!!!!!!
!        UV_Tmajor(i,j,iTC) = Vphase(i,j) !!!!!!!!!!!!!!!!!!
      enddo
    enddo  
#endif

  END DO

  call check( nf90_open(trim( TIDE_FILE ), NF90_WRITE, ncid) )

  varname = 'tide_period'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, Tperiod ) )
  varname = 'tidal_constituents'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, Tconsti ) )

  varname = 'tide_Ephase'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, SSH_Tphase ) )
  varname = 'tide_Eamp'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, SSH_Tamp ) )

#if defined UV_TIDES
  varname = 'tide_Cphase'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, UV_Tphase ) )
  varname = 'tide_Cangle'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, UV_Tangle ) )
  varname = 'tide_Cmin'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, UV_Tminor ) )
  varname = 'tide_Cmax'
  write(*,*)  'Write: ', trim( varname )
  call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
  call check( nf90_put_var(ncid, var_id, UV_Tmajor ) )
#endif

  call check( nf90_close(ncid) )


  write(*,*) 'FINISH!!'
      
      
END PROGRAM frcTIDE2ROMS

