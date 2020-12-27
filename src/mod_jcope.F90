
!!!=== Copyright (c) 2020-2021 Takashi NAKAMURA  =====

!!!**** JCOPE MODULE ************************************

MODULE mod_jcope

  implicit none  
#if defined JCOPE_T
  character(4), parameter :: JCOPE_prefix(7) = (/  &
    "EGT_", "UT_ ", "VT_ ", "*_  ", "*_  ", "TT_ ", "ST_ " /)
#else
  character(3), parameter :: JCOPE_prefix(7) = (/  &
    "EL_", "U_ ", "V_ ", "*_ ", "*_ ", "T_ ", "S_ " /)
#endif

  CONTAINS

!**** create initial conditions NetCDF file **********************************************

  SUBROUTINE read_jcope_info( jcope_info_dir, Nx, Ny, Nz )
    character(len=*), intent( in) :: jcope_info_dir
    integer, intent(out) :: Nx, Ny, Nz

    character(256) :: info_file
    character(5) :: c1
    character(6) :: c2
    real(8), allocatable :: lat(:), lon(:)

    integer :: i

#if defined JCOPE_T
    info_file = trim( jcope_info_dir )//"TT.ctl"
#else
    info_file = trim( jcope_info_dir )//"t.ctl"
#endif

    open(unit=20, file=trim( info_file ), action='read')
    ! Skip header lines
    do i=1, 4
      read(20, *)
    end do  
    ! Read Number of X grids
    read(20, '(a,i4,a)') c1, Nx, c2
    write(*, *) c1, Nx, c2

    do i=1, Nx/6+1
      read(20, *)
    enddo
    read(20, '(a,i4,a)') c1, Ny, c2
    write(*, *) c1, Ny, c2

!    do i=1, Ny/6+1
!      read(20, *)
!    enddo
    Nz = 47
    write(*, *) Nx, Ny, Nz

    close(20)
  END SUBROUTINE read_jcope_info

  SUBROUTINE read_jcope_latlon( jcope_info_dir, Nx, Ny, ruv, lon, lat )
    character(len=*), intent( in) :: jcope_info_dir
    integer, intent( in) :: Nx, Ny
    integer, intent( in) :: ruv ! 1 -> rho-point, 2 -> u-point, 3 -> v-point
    real(8), intent(out) :: lon(Nx), lat(Ny)

    character(256) :: info_file
    integer :: i,j

#if defined JCOPE_T
    if(ruv==1) then
      info_file = trim( jcope_info_dir )//"TT.ctl"
    elseif(ruv==2) then
      info_file = trim( jcope_info_dir )//"UT.ctl"
    elseif(ruv==3) then
      info_file = trim( jcope_info_dir )//"VT.ctl"
    else
      write(*,*) "ruv setting error in read_jcope_info."
      stop
    endif
#else
    if(ruv==1) then
      info_file = trim( jcope_info_dir )//"t.ctl"
    elseif(ruv==2) then
      info_file = trim( jcope_info_dir )//"u.ctl"
    elseif(ruv==3) then
      info_file = trim( jcope_info_dir )//"v.ctl"
    else
      write(*,*) "ruv setting error in read_jcope_info."
      stop
    endif
#endif

    open(unit=20, file=trim( info_file ), action='read')
    ! Skip header lines
    do i=1, 5
      read(20, *)
    end do  
    ! read lon data

    read(20, '(6(1x,f10.6))') (lon(i),i=1,Nx)

!    i = Nx/6
!    read(20, *) lon(i:i+Mod(Nx,6)-1)

    ! Skip header line
    read(20, *)
    ! read lat data

    read(20,'(6(1x,f10.6))') (lat(i),i=1,Ny)

!    i = Ny/6
!    read(20, *) lat(i:i+Mod(Ny,6)-1)

    write(*, *) " Lon range : ", lon(1), lon(Nx)
!    write(*, *) lon
    write(*, *) " Lat range : ", lat(1), lat(Ny)
!    write(*, *) lat

    close(20)
  END SUBROUTINE read_jcope_latlon

  SUBROUTINE read_jcope_depth( jcope_info_dir, Nx, Ny, Nz, z, zz, dz )
    character(len=*), intent( in) :: jcope_info_dir
    integer, intent( in) :: Nx, Ny, Nz
    real(8), intent(out) :: z(Nx, Ny, Nz), zz(Nx, Ny, Nz), dz(Nx, Ny, Nz)

    character(256) :: info_file
    character(4) :: fildsc
    real(4) :: z4(Nx, Ny, Nz), zz4(Nx, Ny, Nz), dz4(Nx, Ny, Nz)
    integer :: ichflg
    integer :: i,j
    integer :: ierr

#if defined JCOPE_T
    info_file = trim( jcope_info_dir )//"sbasic.dat"
#else
    info_file = trim( jcope_info_dir )//"basic.dat"
#endif

    open( 88, file=trim(info_file)                                          &
        , iostat=ierr, status='unknown', action='read'  &
        , access='sequential', form='unformatted', convert='big_endian' )

    read(88,iostat=ierr) fildsc, z4, zz4, dz4, ichflg

    write(*,*) fildsc, ierr, ichflg
    
    if(ierr /= 0 .or. ichflg .ne. 123456)then
      write(*,*)' Error reading unformated data from ', trim(info_file)
      stop
    endif

    z = dble(z4)
    zz = dble(zz4)
    dz = dble(dz4)


!    do j=111,255
!      write(99,*) z4(80:255,j,Nz)
!    enddo
!    do j=1,Ny
!      write(98,*) zz(:,j,Nz-1)
!    enddo
!    do j=1,Ny
!      write(97,*) dz(:,j,Nz-1)
!    enddo
!
    close(88)
  END SUBROUTINE read_jcope_depth  

  SUBROUTINE read_jcope_data2D( jcope_data_file        &
    , Imin, Imax, Jmin, Jmax, is, ie, js, je, data2d )
    character(len=*), intent( in) :: jcope_data_file
    integer, intent( in) :: Imin, Imax, Jmin, Jmax
    integer, intent( in) :: is, ie, js, je
    real(8), intent(out) :: data2d(is:ie, js:je)

    integer :: Nx, Ny
    real(4), allocatable :: data4(:,:)
    integer :: recsize, iol
    real(4) :: r4t
    integer :: i,j
    integer :: ierr
    character(4) :: fildscT
    integer :: iyf,imf,idf,mall

    Nx = Imax - Imin + 1
    Ny = Jmax - Jmin + 1
    
    allocate( data4(Imin:Imax, Jmin:Jmax) )

#if defined JCOPE_T

    inquire(iolength=iol) r4t
    recsize = Nx*Ny*iol

    open( 88, file=trim(jcope_data_file)                                  &
        , form=   'unformatted', access= 'direct', recl=   recsize        &
        , iostat= ierr, action= 'read', convert='big_endian' )
    write(*,*) ierr
#else

    open( 88, file=trim(jcope_data_file)                                  &
        , form=   'unformatted'                                           &
        , iostat= ierr, action= 'read', convert='big_endian' )
    write(*,*) ierr
    read(88,iostat=ierr) fildscT, iyf, imf, idf                           &
        , data4, mall
    write(*,*) ierr, fildscT, iyf, imf, idf, mall
#endif

    read(88,rec=1,iostat=ierr) data4   
    write(*,*) ierr

    data2d = dble( data4(is:ie, js:je) )

!    do j=js,je
!      write(98,*) data2d(:,j)
!    enddo

    close(88)
  END SUBROUTINE read_jcope_data2D

  SUBROUTINE read_jcope_data3D( jcope_data_file             &
    , Imin, Imax, Jmin, Jmax, is, ie, js, je, data3d )
    character(len=*), intent( in) :: jcope_data_file
    integer, intent( in) :: Imin, Imax, Jmin, Jmax
    integer, intent( in) :: is, ie, js, je
    real(8), intent(out) :: data3d(is:ie, js:je, 1:46)

    integer :: Nx, Ny, Nz
    real(4), allocatable :: data4(:,:,:)
    integer :: recsize, iol
    real(4) :: r4t
    integer :: i,j,k
    integer :: ierr
    character(4) :: fildscT
    integer :: iyf,imf,idf,mall

    Nx = Imax - Imin + 1
    Ny = Jmax - Jmin + 1
    Nz = 47
    
    allocate( data4(Imin:Imax, Jmin:Jmax, 1:47) )

#if defined JCOPE_T

    inquire(iolength=iol) r4t
    recsize = Nx*Ny*Nz*iol

    open( 88, file=trim(jcope_data_file)                                  &
        , form=   'unformatted', access= 'direct', recl=   recsize        &
        , iostat= ierr, action= 'read', convert='big_endian' )
    write(*,*) ierr
    read(88,rec=1,iostat=ierr) data4
    write(*,*) ierr
#else

    open( 88, file=trim(jcope_data_file)                                  &
        , form=   'unformatted'                                           &
        , iostat= ierr, action= 'read', convert='big_endian' )
    write(*,*) ierr
    read(88,iostat=ierr) fildscT, iyf, imf, idf                           &
        , (((data4(i,j,k),i=Imin,Imax),j=Jmin,Jmax),k=1,46), mall
    write(*,*) ierr, fildscT, iyf, imf, idf, mall
#endif

    do k=1, 46
      data3d(:,:,k) = dble( data4(is:ie, js:je, 47-k) )
    enddo

!    do j=js,je
!      write(98,*) data2d(:,j)
!    enddo

    close(88)
  END SUBROUTINE read_jcope_data3D
    
END MODULE mod_jcope
      
! -------------------------------------------------------------------------
