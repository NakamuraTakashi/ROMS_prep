
!!!=== Copyright (c) 2025 Takashi NAKAMURA  =====

!!!**** MOVE-JPN MODULE ************************************

MODULE mod_movejpn

  implicit none  

  CONTAINS

  SUBROUTINE read_movejpn_ctl( topo_dir, Nx, Ny, Nz, lon, lat, z )
    character(len=*), intent( in) :: topo_dir
    integer, intent(out) :: Nx, Ny, Nz
    real(8), intent(out), allocatable :: lon(:), lat(:), z(:)

    character(256) :: info_file
    character(4) :: c1
    character(7) :: c2
    integer :: i

    info_file = trim( topo_dir )//"topo.ctl"

    open(unit=20, file=trim( info_file ), action='read')
    ! Skip header lines
    do i=1, 5
      read(20, *)
    end do  
    ! Read Number of X grids
    read(20, '(a,i7,a)') c1, Nx, c2
    write(*, *) c1, Nx, c2
    allocate( lon(Nx) )
    ! read lon data
    read(20, '(5(1x,f11.5))') (lon(i),i=1,Nx)

    ! Read Number of y grids
    read(20, '(a,i7,a)') c1, Ny, c2
    write(*, *) c1, Ny, c2
    allocate( lat(Ny) )
    ! read lat data
    read(20, '(5(1x,f11.5))') (lat(i),i=1,Ny)

    ! Read Number of z grids
    read(20, '(a,i7,a)') c1, Nz, c2
    write(*, *) c1, Nz, c2
    allocate( z(Nz) )
    ! read lat data
    read(20, '(5(1x,f11.5))') (z(i),i=1,Nz)

    write(*, *) "Lon range: ", lon(1), lon(Nx)
    write(*, *) "Lat range: ", lat(1), lat(Ny)
    write(*, *) "Z range  : ", z(1), z(Nz)

    close(20)
  END SUBROUTINE read_movejpn_ctl

  SUBROUTINE read_movejpn_depth( topo_dir, Nx, Ny, dp, lv )
    character(len=*), intent( in) :: topo_dir
    integer, intent( in) :: Nx, Ny
    real(8), intent(out), allocatable :: dp(:,:)
    integer, intent(out), allocatable :: lv(:,:)

    character(256) :: info_file
!    real(4) :: dp4(Nx, Ny), lv4(Nx, Ny)
    integer(4) :: i4dp(Nx, Ny), i4lv(Nx, Ny)
    character(1) :: header_dummy(4)
    integer :: ierr

    info_file = trim( topo_dir )//"topo.d"

    allocate( dp(Nx, Ny), lv(Nx, Ny) )

    open( 88, file=trim(info_file)                                          &
        , iostat=ierr, status='old', action='read'  &
        , access='stream', form='unformatted', convert='big_endian' )
    if(ierr /= 0)then
      write(*,*) 'Error: ',ierr, trim(info_file)
      stop
    endif

    read(88) header_dummy
!    write(*,*) header_dummy
    
    read(88) i4dp
    read(88) i4lv

    close(88)

    dp = dble(i4dp)*1.0d-2  ! convert from cm to m
    lv = int(i4lv)

  END SUBROUTINE read_movejpn_depth  
    
END MODULE mod_movejpn
      
! -------------------------------------------------------------------------
