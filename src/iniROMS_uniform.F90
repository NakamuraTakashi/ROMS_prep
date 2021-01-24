
!!!=== Copyright (c) 2020-2021 Takashi NAKAMURA  ===== 

PROGRAM iniROMS_uniform
  use netcdf
  use mod_utility
  use mod_roms_netcdf
  use mod_calendar
  use mod_interpolation
 
  implicit none
      
! ---------------------------------------------------------------------
  integer :: Syear, Smonth, Sday
  integer :: Ryear, Rmonth, Rday
  integer :: itime
  integer :: INIyear, INImonth, INIday, INIhour
  integer :: mode
  character(256) :: GRID_FILE
  character(256) :: ROMS_HISFILE
  integer :: romsvar(N_var)
  character(256) :: INI_prefix
  integer :: N_s_rho
  integer :: Nzr
  integer :: spherical
  integer :: Vtransform, Vstretching
  real(8) :: THETA_S , THETA_B, TCLINE, DCRIT
 
  character(33) :: TIME_ATT  = "seconds since 2000-01-01 00:00:00"
  character(15) :: INI_suffix   = "_20000101.00.nc"
  character(256) :: INI_FILE
  
  real(8), allocatable :: time_all(:), lat_all(:), lon_all(:)
  integer :: start1D(1), count1D(1)
  integer :: start2D(2), count2D(2)
  integer :: start3D(3), count3D(3)
  integer :: start4D(4), count4D(4)
  
  real(8) :: hc       
  real(8), allocatable :: sc_w(:)       
  real(8), allocatable :: sc_r(:)  
  real(8), allocatable :: Cs_w(:)       
  real(8), allocatable :: Cs_r(:)  
  real(8), allocatable :: z_r(:,:,:)
  real(8), allocatable :: z_w(:,:,:)
  real(8), allocatable :: z_u(:,:,:)
  real(8), allocatable :: z_v(:,:,:)
  
  real(8), allocatable :: zeta(:,:,:)   ! free-surface (meter)
  real(8), allocatable :: ubar(:,:,:)   ! vertically integrated u-momentum component (meter second-1)
  real(8), allocatable :: vbar(:,:,:)   ! vertically integrated v-momentum component (milibar=hPa)
  real(8), allocatable :: u(:,:,:,:)    ! u-momentum component (meter second-1)
  real(8), allocatable :: v(:,:,:,:)    ! v-momentum component (meter second-1)
  real(8), allocatable :: t(:,:,:,:)    ! tracer 

  real(8) :: ocean_time(1)                 ! Ocean time (sec)

  integer :: i,j,k
  integer :: idt,incdf

  integer :: jdate_Start, jdate_Ref
  real(8) :: d_jdate_Start, d_jdate_Ref
  
  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(11) :: YYYYMMDDpHH

  integer :: Nxr, Nyr, Nxu, Nyu, Nxv, Nyv
  integer :: L, M, N

  integer :: ncid,var_id
  integer :: ncid2,var_id2
  !
  namelist/grd/GRID_FILE
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/roms2roms/ROMS_HISFILE, romsvar
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
  read (5, nml=refdate)
  rewind(5)
  read (5, nml=roms2roms)
  rewind(5)
  read (5, nml=ini)
  rewind(5)
  read (5, nml=hcoord)
  rewind(5)
  read (5, nml=zcoord)
  rewind(5)
           
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
  
  allocate( sc_w(0:N) )
  allocate( sc_r(1:N) )
  allocate( Cs_w(0:N) )
  allocate( Cs_r(1:N) )
  allocate( z_r(0:L, 0:M, 1:N) )
  allocate( z_w(0:L, 0:M, 0:N) )
  allocate( z_u(1:L, 0:M, 1:N) )
  allocate( z_v(0:L, 1:M, 1:N) )
  
  allocate( zeta(0:L, 0:M, 1) )
  allocate( ubar(1:L, 0:M, 1) )
  allocate( vbar(0:L, 1:M, 1) )
  allocate( u(1:L, 0:M, 1:N, 1) )
  allocate( v(0:L, 1:M, 1:N, 1) )
  allocate( t(0:L, 0:M, 1:N, 1) )
!  allocate( salt(0:L, 0:M, 1:N, 1) )
  
  
  ! Get variables
  start2D = (/ 1,  1 /)
  count2D = (/ Nxr, Nyr /)
  ! Close NetCDF file
  call check( nf90_close(ncid) )
  
  write(*,*) "CLOSE: ", trim( GRID_FILE )
      
  call set_scoord (  Nzr, Vtransform, Vstretching    &
        , THETA_S, THETA_B, TCLINE, DCRIT            &
!    output parameters
        , hc, sc_w, sc_r, Cs_w, Cs_r )

  write(*,*) "*********************************************"
!-Read ROMS ocean_his netCDF file --------------------------------
  call jd(INIyear, INImonth, INIday, jdate_Start)
  d_jdate_Start = dble(jdate_Start) + dble(INIhour)/24.0d0
  write(*,*) d_jdate_Start

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
  write(*,*) d_jdate_Ref

  ocean_time(1) = (d_jdate_Start - d_jdate_Ref)*86400.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-Initial netcdf file name -------------------------------------------------
  call oceantime2cdate(ocean_time(1), Ryear, Rmonth, Rday, YYYYMMDDpHH)
  INI_suffix(2:12)=YYYYMMDDpHH
  INI_FILE = trim( INI_prefix )//INI_suffix

!-Create the ROMS initial conditions netCDF file --------------------------------
      
  call createNetCDFini2(  trim( ROMS_HISFILE ),  trim( INI_FILE )   &
        , TIME_ATT , Nxr, Nyr, Nzr, 1 ,romsvar  )

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
  write(*,*)  'Write: ocean_time'
  call check( nf90_inq_varid(ncid2, 'ocean_time', var_id2) )
  call check( nf90_put_var(ncid2, var_id2, ocean_time, start = start1D, count = count1D) )
      
  if( romsvar(1)==1 ) then
  !- zeta --------------------------------
    zeta = 1.484d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    start3D = (/ 1,  1,  1 /)
    count3D = (/ Nxr, Nyr, 1 /)
    write(*,*)  'Write: zeta'
    call check( nf90_inq_varid(ncid2, 'zeta', var_id2) )
    call check( nf90_put_var(ncid2, var_id2, zeta, start = start3D, count = count3D) )
  endif

  if( romsvar(2)==1 .or. romsvar(3)==1 .or.         &
      romsvar(4)==1 .or. romsvar(5)==1      ) then
   
    if( romsvar(2)==1 ) then
      !- Write u v --------------------------------
      u = 0.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ Nxu, Nyu, Nzr, 1 /)
      write(*,*)  'Write: u'
      call check( nf90_inq_varid(ncid2, 'u', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, u, start = start4D, count = count4D) )
    endif
    if( romsvar(3)==1 ) then
      v = 0.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      start4D = (/ 1,  1,  1,  1 /)
      count4D = (/ Nxv, Nyv, Nzr, 1 /)
      write(*,*)  'Write: v'
      call check( nf90_inq_varid(ncid2, 'v', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, v, start = start4D, count = count4D) )
    endif
  
  !-Depth averaged velocity calculation --------------------------------
  
    if( romsvar(4)==1 ) then
      ubar(:,:,1)=0.0d0
      do i=1,L
        do j=0,M
          do k=1,N
            ubar(i,j,1)= ubar(i,j,1)+u(i,j,k,1)*abs(Cs_w(k-1)-Cs_w(k))
          enddo
        enddo
      enddo
      start3D = (/ 1,  1,  1 /)
      count3D = (/ Nxu, Nyu, 1 /)
      write(*,*)  'Write: ubar'
      call check( nf90_inq_varid(ncid2, 'ubar', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, ubar, start = start3D, count = count3D) )
    endif

    if( romsvar(5)==1 ) then   
      vbar(:,:,1)=0.0d0
      do i=0,L
        do j=1,M
          do k=1,N
            vbar(i,j,1)= vbar(i,j,1)+v(i,j,k,1)*abs(Cs_w(k-1)-Cs_w(k))
          enddo
        enddo
      enddo
      start3D = (/ 1,  1,  1 /)
      count3D = (/ Nxv, Nyv, 1 /)
      write(*,*)  'Write: vbar'
      call check( nf90_inq_varid(ncid2, 'vbar', var_id2) )
      call check( nf90_put_var(ncid2, var_id2, vbar, start = start3D, count = count3D) )
    endif
  
  endif
  
!- Tracer (temp, salt, etc.) --------------------------------

  do i=6, N_var 
      
    if( romsvar(i)==0 ) cycle
    if( i==6 ) then
      t = 29.3d0 !!! temp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif(i==7) then
      t = 33.3d0 !!! salt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif

    start4D = (/ 1,  1,  1,  1 /)
    count4D = (/ Nxr, Nyr, Nzr, 1 /)
    write(*,*)  'Write: ', trim( VAR_NAME(i) )
    call check( nf90_inq_varid(ncid2, trim( VAR_NAME(i) ), var_id2) )
    call check( nf90_put_var(ncid2, var_id2, t, start = start4D, count = count4D) )
  
  enddo

  call check( nf90_close(ncid2) )

  write(*,*) 'FINISH!!'
      
      
END PROGRAM iniROMS_uniform
      
